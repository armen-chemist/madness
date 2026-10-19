/// \file magneticatom.cc
/// \brief atoms and small molecules in a uniform magnetic field B (along z), 3D: one electron,
///        closed-shell RHF (2 electrons) or spin-polarized UHF (all spins anti-parallel to B).
///        Every orbital is updated with the magnetic Green's function (magnetic_green.h).
///
///   H = 1/2 (p + A)^2 + Vnuc = -1/2 nabla^2 + (B/2) Lz + (B^2/8)(x^2 + y^2) + Vnuc,   gauge origin r = 0
///   Vnuc = MADNESS smoothed nuclear potential (smoothing radius from eprec, as in moldft)
///   psi = -G_B(E) Veff psi,  G_B(E) = (1/2 (p + A)^2 - E)^{-1},  E < |B|/2
///
/// Spin Zeeman term: not in the one-electron and RHF energies; included in E_phys of UHF.
///
/// Usage:  magneticatom B=.. [atoms=Z,x,y,z;Z,x,y,z...] [L=20] [k=8] [thresh=1e-6] [eprec=1e-6]
///                      [E0=-0.5] [maxiter=60] [tol=1e-7] [ref=...] [m=0]
///                      [qextra=10] [rcut=..] [rtaper=..] [dphi=3] [h=0.25] [tmin=1e-8] [tmax_factor=45] [qbatch=64]
///   nelec=2: closed-shell RHF;  nelec=N state=highspin (or triplet): spin-polarized UHF with start
///   orbitals from mlist=m1,..,mN (axial=1, E0list=e1,..,eN), mloc=m1,..,mk (axial=0: k orbitals with
///   Lz = m_j about every nucleus), or one 1s per nucleus (axial=0, N nuclei);
///   ropt=1: bond-length (breathing) optimization with Hellmann-Feynman forces
/// Examples: hydrogen atom (default atoms=1,0,0,0); H2+ along B: atoms=1,0,0,-1;1,0,0,1;
///           H2+ perpendicular to B: atoms=1,-1,0,0;1,1,0,0

#include <madness/mra/mra.h>
#include <madness/mra/magnetic_green.h>
#include <madness/mra/nonlinsol.h>
#include <madness/mra/vmra.h>
#include <madness/chem/molecule.h>
#include <madness/chem/potentialmanager.h>
#include <madness/chem/molecular_functors.h>
#include <chrono>
#include <cstdio>
#include <map>
#include <memory>
#include <sstream>
#include <string>

using namespace madness;

namespace {

std::map<std::string, std::string> parse(int argc, char** argv) {
    std::map<std::string, std::string> a;
    for (int i = 1; i < argc; ++i) {
        std::string s(argv[i]);
        auto eq = s.find('=');
        if (eq != std::string::npos) a[s.substr(0, eq)] = s.substr(eq + 1);
    }
    return a;
}

double get(const std::map<std::string, std::string>& a, const std::string& key, double def) {
    auto it = a.find(key);
    return it == a.end() ? def : std::stod(it->second);
}

/// first-derivative stencil of Derivative<>: ABGV (default), BLE or smoothed B-spline
enum class Stencil { abgv, ble, bspline };

const char* stencil_name(Stencil t) {
    return t == Stencil::abgv ? "abgv" : (t == Stencil::ble ? "ble" : "bspline");
}

/// free-space d/dx_axis with the chosen stencil
Derivative<double_complex, 3> derivative(World& world, int axis, Stencil t) {
    auto D = free_space_derivative<double_complex, 3>(world, axis);
    if (t == Stencil::ble) D.set_ble1();
    if (t == Stencil::bspline) D.set_bspline1();
    return D;
}

/// Lz f = -i (x df/dy - y df/dx), angular momentum about the z-axis through r = 0
complex_function_3d Lz(World& world, const complex_function_3d& f, Stencil t) {
    real_function_3d x = real_factory_3d(world).functor([](const coord_3d& r) { return r[0]; });
    real_function_3d y = real_factory_3d(world).functor([](const coord_3d& r) { return r[1]; });
    complex_function_3d r = x * derivative(world, 1, t)(f) - y * derivative(world, 0, t)(f);
    r.scale(double_complex(0.0, -1.0));
    r.truncate();
    return r;
}

/// one orbital in the fixed potential V: psi~ = -G_B(E) V psi, Rayleigh-quotient update
/// E <- E + <psi~|V|psi~ - psi> / <psi~|psi~>, KAIN step on the residual psi - psi~, normalize.
/// On return psi is the normalized orbital; returns E.
double solve(World& world, const MagneticGreen& G, const real_function_3d& V, complex_function_3d& psi,
             double E, int maxiter, double tol, bool verbose, int* iterations = nullptr) {
    psi.scale(1.0 / psi.norm2());
    auto Vf = [&](const complex_function_3d& f) {
        complex_function_3d r = V * f;
        r.truncate();
        return r;
    };
    const double top = G.bottom();
    XNonlinearSolver<std::vector<complex_function_3d>, double_complex, vector_function_allocator<double_complex, 3>>
        solver(vector_function_allocator<double_complex, 3>(world, 1));
    solver.set_maxsub(8);
    int it = 0;
    for (; it < maxiter; ++it) {
        auto t0 = std::chrono::steady_clock::now();
        if (E >= top) E = top - 1e-3 * std::max(top, 1e-3);      // G_B(E) needs E < |B|/2
        complex_function_3d pt = G.green(Vf(psi), E);
        pt.scale(-1.0);
        const double dE = std::real(inner(pt, Vf(pt - psi))) / std::real(inner(pt, pt));
        complex_function_3d res = psi - pt;                      // zero at the fixed point
        const double resid = res.norm2();
        std::vector<complex_function_3d> unew =
            solver.update(std::vector<complex_function_3d>{psi}, std::vector<complex_function_3d>{res});
        psi = unew[0];
        psi.scale(1.0 / psi.norm2());
        psi.truncate();
        E += dE;
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (verbose && world.rank() == 0) {
            std::printf("  iter %2d  E = %.12f  dE = %+.3e  |psi - psi~| = %.3e  (%.1f s)\n", it, E, dE, resid, sec);
            std::fflush(stdout);
        }
        if (std::abs(dE) < tol && resid < std::sqrt(tol)) break;
    }
    if (iterations) *iterations = it + 1;
    return E;
}

}  // namespace

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    int status = 0;
    try {
        startup(world, argc, argv);
        const auto args = parse(argc, argv);
        const double B = get(args, "B", 1.0);
        const double L = get(args, "L", 20.0);
        const int k = int(get(args, "k", 8));
        const double thresh = get(args, "thresh", 1e-6);
        const double eprec = get(args, "eprec", 1e-6);
        const double ref = get(args, "ref", std::nan(""));

        FunctionDefaults<3>::set_cubic_cell(-L, L);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_autorefine(true);
        FunctionDefaults<3>::set_truncate_mode(1);
        FunctionDefaults<3>::set_initial_level(5);

        Molecule mol;
        std::string atoms = args.count("atoms") ? args.at("atoms") : "1,0,0,0";
        std::stringstream ss(atoms);
        std::string a;
        while (std::getline(ss, a, ';')) {
            double Z, x, y, z;
            char c;
            std::stringstream as(a);
            as >> Z >> c >> x >> c >> y >> c >> z;
            mol.add_atom(x, y, z, Z, int(Z));
        }
        mol.update_rcut_with_eprec(eprec);

        real_function_3d V = real_factory_3d(world).functor(MolecularPotentialFunctor(mol)).truncate_on_project();
        V.truncate();

        const bool axial = get(args, "axial", 1) != 0;   // UHF start orbitals: about the B axis (1) or per nucleus (0)
        const int maxiter = int(get(args, "maxiter", 60));
        const double tol = get(args, "tol", 1e-7);

        // <psi|H|psi> for normalized p with derivative stencil t:
        // kin = 1/2 sum_d <D_d p|D_d p>, para = (B/2) <p|Lz p>, dia = (B^2/8) <|p|^2 rho^2>, pot = <|p|^2 V>
        real_function_3d rho2 = real_factory_3d(world).functor(
            [](const coord_3d& r) { return r[0] * r[0] + r[1] * r[1]; });
        struct Hparts { double kin, para, dia, pot, E; };
        auto expect_H = [&](const complex_function_3d& p, Stencil t) {
            Hparts h{0.0, 0.0, 0.0, 0.0, 0.0};
            for (int ax = 0; ax < 3; ++ax) {
                auto D = derivative(world, ax, t);
                complex_function_3d dp = D(p);
                h.kin += 0.5 * std::real(inner(dp, dp));
            }
            h.para = 0.5 * B * std::real(inner(p, Lz(world, p, t)));
            real_function_3d dens = abs_square(p);
            h.dia = B * B / 8.0 * inner(dens, rho2);
            h.pot = inner(dens, V);
            h.E = h.kin + h.para + h.dia + h.pot;
            return h;
        };

        // guess: sum of 1s-like functions at the nuclei, squeezed transversely to the magnetic
        // length (lowest Landau orbital e^{-|B| rho^2/4}), times (x -+ i y)^|m| for m != 0
        const auto coords = mol.get_all_coords_vec();
        const int m = int(get(args, "m", 0));
        const double Bq = 0.25 * std::abs(B);
        complex_function_3d psi = complex_factory_3d(world).functor([coords, m, Bq](const coord_3d& r) {
            // each nucleus off the z-axis gets the magnetic-translation phase e^{-i (B/2)(X y - Y x)}:
            // A(r) = A(r - R) + grad chi, chi = (B/2)(X y - Y x), and pi = p + A for the electron
            double_complex s = 0;
            for (const auto& R : coords) {
                const double dx = r[0] - R[0], dy = r[1] - R[1];
                const double chi = 2 * Bq * (R[0] * r[1] - R[1] * r[0]);
                s += std::exp(double_complex(-(r - R).normf() - Bq * (dx * dx + dy * dy), -chi));
            }
            double_complex ang = 1.0;
            for (int i = 0; i < std::abs(m); ++i) ang *= double_complex(r[0], m < 0 ? -r[1] : r[1]);
            return ang * s;
        });

        MagneticGreenParams gp;
        gp.verbose = get(args, "lverbose", 0) != 0;
        gp.qextra = get(args, "qextra", 10.0);
        gp.rcut = get(args, "rcut", 0.0);
        gp.rtaper = get(args, "rtaper", 0.0);
        gp.dphi = get(args, "dphi", 3.0);
        gp.h = get(args, "h", 0.25);                       // proper-time step in u = ln(tau)
        gp.tmin = get(args, "tmin", 1e-8);                 // smallest proper-time node
        gp.tmax_factor = get(args, "tmax_factor", 45.0);   // tau_max = tmax_factor / (|B|/2 - E)
        gp.qbatch = long(get(args, "qbatch", 64));         // q-modes per vector operation (0: all at once)
        const MagneticGreen G(world, B, gp);
        if (world.rank() == 0) {
            print("magneticatom: B =", B, " atoms =", atoms, " k =", k, " thresh =", thresh, " eprec =", eprec,
                  " cell [", -L, ",", L, "]^3");
            print("Fourier approach: modes =", G.qgrid().size());
            print("            Vnuc at (0.5,0,0) =", V(coord_3d{0.5, 0.0, 0.0}));
        }
        const int nelec = int(get(args, "nelec", 1));
        MADNESS_CHECK(nelec >= 1);
        // state=triplet (nelec = 2) or state=highspin (any nelec): all spins anti-parallel to B
        const bool triplet = args.count("state") && (args.at("state") == "triplet" || args.at("state") == "highspin");
        MADNESS_CHECK(nelec <= 2 || triplet);
        if (nelec >= 2 && triplet) {
            // spin-polarized UHF, all N spins anti-parallel to B. Same-spin Fock operator for phi_i:
            //   F_i = h + sum_{j != i} (J_j - K_j),  J_j = Coul(|phi_j|^2),  K_j f = phi_j Coul(phi_j^* f)
            //   phi_i = -G_B(eps_i) W_i phi_i,  W_i = Vnuc + sum_{j != i} (J_j - K_j)   (Gauss-Seidel over i,
            //   phi_i Gram-Schmidt orthogonalized to phi_0..phi_{i-1})
            //   E_el = sum_i <h>_i + sum_{i<j} (J_ij - K_ij);  E_phys = E_el + E_nn - N |B|/2 (spin Zeeman)
            // start orbitals: axial = 1, Lz = m_i about the B axis (e.g. H2 along B: m = 0, -1, Lai, RMP 73,
            //   629, Sec. IV); axial = 0: sigma_g / sigma_u (N = 2), one 1s per nucleus (N nuclei, e.g.
            //   quartet H3: Bischoff, PRA 101, 053413 (2020)) or atomic orbitals on every nucleus (mloc)
            const auto coul = CoulombOperator(world, get(args, "lo", 1e-6), 0.1 * thresh);
            auto coulc = [&](const complex_function_3d& f) {   // Coul(f) for complex f: real and imaginary parts
                real_function_3d r = apply(coul, real(f)), i = apply(coul, imag(f));
                complex_function_3d c = convert<double, double_complex>(r) +
                                        convert<double, double_complex>(i) * double_complex(0.0, 1.0);
                c.truncate();
                return c;
            };
            const int ma = int(get(args, "ma", 0)), mb = int(get(args, "mb", -1));
            // comma-separated lists: mlist = Lz of the axial start orbitals (any N, e.g. He2: 0,-1,-2,-3),
            // E0list = their starting orbital energies (each below |B|/2)
            auto getlist = [&](const std::string& key) {
                std::vector<double> v;
                if (!args.count(key)) return v;
                std::stringstream ls(args.at(key));
                std::string t;
                while (std::getline(ls, t, ',')) v.push_back(std::stod(t));
                return v;
            };
            const std::vector<double> mlist = getlist("mlist"), e0list = getlist("E0list"), mloc = getlist("mloc");
            // guesses: Landau-squeezed 1s sum at the nuclei times (x -+ i y)^|m|
            auto guess = [&](int mm) {
                return complex_factory_3d(world).functor([coords, mm, Bq](const coord_3d& r) {
                    double s = 0;
                    for (const auto& R : coords)
                        s += std::exp(-(r - R).normf() - Bq * ((r[0] - R[0]) * (r[0] - R[0]) + (r[1] - R[1]) * (r[1] - R[1])));
                    double_complex ang = 1.0;
                    for (int i = 0; i < std::abs(mm); ++i) ang *= double_complex(r[0], mm < 0 ? -r[1] : r[1]);
                    return ang * s;
                });
            };
            // general scheme: sigma_g / sigma_u from 1s-like functions at the nuclei, each with the
            // magnetic-translation phase e^{-i (B/2)(X y - Y x)} of the main guess
            auto guess_su = [&](double sgn) {
                return complex_factory_3d(world).functor([coords, Bq, sgn](const coord_3d& r) {
                    double_complex s = 0;
                    double c = 1.0;
                    for (const auto& R : coords) {
                        const double dx = r[0] - R[0], dy = r[1] - R[1];
                        const double chi = 2 * Bq * (R[0] * r[1] - R[1] * r[0]);
                        s += c * std::exp(double_complex(-(r - R).normf() - Bq * (dx * dx + dy * dy), -chi));
                        c *= sgn;
                    }
                    return s;
                });
            };
            // N > 2: one phased 1s per nucleus
            auto guess_at = [&](const coord_3d& R) {
                return complex_factory_3d(world).functor([R, Bq](const coord_3d& r) {
                    const double dx = r[0] - R[0], dy = r[1] - R[1];
                    const double chi = 2 * Bq * (R[0] * r[1] - R[1] * r[0]);
                    return std::exp(double_complex(-(r - R).normf() - Bq * (dx * dx + dy * dy), -chi));
                });
            };
            // atomic start (axial=0, mloc=m1,..,mk): k orbitals per nucleus A, Lz = m_j about A,
            //   (dx -+ i dy)^|m| e^{-|r - R_A| - (B/4) rho_A^2} e^{-i (B/2)(X y - Y x)}, N = k x (number of nuclei)
            auto guess_loc = [&](const coord_3d& R, int mm) {
                return complex_factory_3d(world).functor([R, Bq, mm](const coord_3d& r) {
                    const double dx = r[0] - R[0], dy = r[1] - R[1];
                    const double chi = 2 * Bq * (R[0] * r[1] - R[1] * r[0]);
                    double_complex ang = 1.0;
                    for (int i = 0; i < std::abs(mm); ++i) ang *= double_complex(dx, mm < 0 ? -dy : dy);
                    return ang * std::exp(double_complex(-(r - R).normf() - Bq * (dx * dx + dy * dy), -chi));
                });
            };
            const int N = nelec;
            // start orbitals: mlist given (axial=1): guess(m_i), any N and any nuclei; mloc (axial=0): atomic start;
            // N = 2 without mlist: guess(ma), guess(mb) (axial=1) or sigma_g, sigma_u (axial=0);
            // N > 2 without mlist: one phased 1s per nucleus (axial=0, N nuclei)
            std::vector<int> mm;
            if (!mlist.empty()) {
                MADNESS_CHECK(axial && int(mlist.size()) == N);
                for (double x : mlist) mm.push_back(int(x));
            } else if (N == 2) {
                mm = {ma, mb};
            }
            if (!mloc.empty()) {
                MADNESS_CHECK(mlist.empty() && !axial && int(mloc.size() * coords.size()) == N);
                mm.clear();
                for (std::size_t A = 0; A < coords.size(); ++A)
                    for (double x : mloc) mm.push_back(int(x));
            }
            MADNESS_CHECK(!mm.empty() || (!axial && int(coords.size()) == N));
            std::vector<complex_function_3d> phi;
            if (!mloc.empty()) {
                for (const auto& R : coords)
                    for (double x : mloc) phi.push_back(guess_loc(R, int(x)));
            } else if (!mlist.empty()) for (int x : mm) phi.push_back(guess(x));
            else if (N == 2) phi = axial
                ? std::vector<complex_function_3d>{guess(ma), guess(mb)}
                : std::vector<complex_function_3d>{guess_su(1.0), guess_su(-1.0)};
            else
                for (const auto& R : coords) phi.push_back(guess_at(R));
            if (mm.empty()) mm.assign(N, 0);
            for (int i = 0; i < N; ++i) {   // Gram-Schmidt (localized or off-axis guesses overlap)
                for (int j = 0; j < i; ++j) phi[i] = phi[i] - phi[j] * inner(phi[j], phi[i]);
                phi[i].scale(1.0 / phi[i].norm2());
            }
            std::vector<double> eps(N, get(args, "E0", -0.5));
            if (N == 2) eps[1] = get(args, "E0b", eps[0]);
            if (!e0list.empty()) {
                MADNESS_CHECK(int(e0list.size()) == N);
                eps = e0list;
            }
            for (double e : eps) MADNESS_CHECK(e < 0.5 * std::abs(B));
            if (world.rank() == 0) {
                std::printf("geometry (bohr):");
                for (const auto& R : coords) std::printf("  (%.6f, %.6f, %.6f)", R[0], R[1], R[2]);
                std::printf("\nnuclear repulsion = %.10f   start Lz =", mol.nuclear_repulsion_energy());
                for (int x : mm) std::printf(" %d", x);
                std::printf("   start eps =");
                for (double e : eps) std::printf(" %.4f", e);
                std::printf("\n");
                std::fflush(stdout);
            }
            const int inner_it = int(get(args, "inner", 3));
            const int maxscf = int(get(args, "maxscf", 40));
            const double tolscf = get(args, "tolscf", 1e-6);
            // geometry (>= 2 nuclei): Hellmann-Feynman forces after each converged SCF; ropt = 1 scales the
            // structure about its centroid (breathing coordinate R = |R_1 - R_0|, secant on dE/dR) and restarts
            // the SCF from the current orbitals, until |dE/dR| < gtol or the step < rtol
            const bool ropt = get(args, "ropt", 0) != 0;
            const int maxgeo = ropt ? int(get(args, "maxgeo", 12)) : 1;
            const double gtol = get(args, "gtol", 1e-4), rtol = get(args, "rtol", 1e-4);
            const double dr0 = get(args, "dr0", 0.02), maxstep = get(args, "maxstep", 0.05);
            double Rprev = 0.0, gprev = 0.0;
            for (int geo = 0; geo < maxgeo; ++geo) {
                double Eold = 0.0, dElast = 1.0;
                auto t0 = std::chrono::steady_clock::now();
                if (world.rank() == 0) print("UHF high-spin, N =", N, ":", axial ? "axial, m_a =" : "general (sigma_g, sigma_u), m_a =", ma, " m_b =", mb, " inner iterations =", inner_it);
                for (int scf = 0; scf < maxscf; ++scf) {
                    auto ts = std::chrono::steady_clock::now();
                    for (int i = 0; i < N; ++i) {      // Gauss-Seidel: phi_i sees the updated phi_0..phi_{i-1}
                        real_function_3d rhoj = real_factory_3d(world);
                        for (int j = 0; j < N; ++j)
                            if (j != i) rhoj = rhoj + abs_square(phi[j]);
                        real_function_3d Veff = V + apply(coul, rhoj);
                        Veff.truncate();
                        auto W = [&](const complex_function_3d& f) {   // (Vnuc + sum_j (J_j - K_j)) f
                            complex_function_3d r = Veff * f;
                            for (int j = 0; j < N; ++j)
                                if (j != i) r = r - phi[j] * coulc(conj(phi[j]) * f);
                            r.truncate();
                            return r;
                        };
                        // phi_i <- -G_B(eps_i) W phi_i, eps_i by the Rayleigh-quotient update, Gram-Schmidt
                        const double top = G.bottom();
                        for (int it = 0; it < inner_it; ++it) {
                            double Em = eps[i];
                            if (Em >= top) Em = top - 1e-3 * std::max(top, 1e-3);
                            complex_function_3d pt = G.green(W(phi[i]), Em);
                            pt.scale(-1.0);
                            const double dE = std::real(inner(pt, W(pt - phi[i]))) / std::real(inner(pt, pt));
                            phi[i] = pt;
                            for (int j = 0; j < i; ++j) phi[i] = phi[i] - phi[j] * inner(phi[j], phi[i]);
                            phi[i].scale(1.0 / phi[i].norm2());
                            phi[i].truncate();
                            eps[i] = Em + dE;
                        }
                    }
                    double Jab = 0.0, Kab = 0.0, hsum = 0.0, ovl = 0.0;
                    for (int a = 0; a < N; ++a) {
                        hsum += expect_H(phi[a], Stencil::abgv).E;
                        const real_function_3d ja = apply(coul, abs_square(phi[a]));
                        for (int b = a + 1; b < N; ++b) {
                            Jab += inner(abs_square(phi[b]), ja);
                            Kab += std::real(inner(phi[b], phi[a] * coulc(conj(phi[a]) * phi[b])));
                            ovl = std::max(ovl, std::abs(inner(phi[a], phi[b])));
                        }
                    }
                    const double Etot = hsum + Jab - Kab + mol.nuclear_repulsion_energy();
                    const double dE = Etot - Eold;
                    Eold = Etot;
                    dElast = dE;
                    const double ezee = 0.5 * N * std::abs(B);   // spin Zeeman, all spins down
                    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
                    if (world.rank() == 0) {
                        std::printf("  scf %2d  eps =", scf);
                        for (int a = 0; a < N; ++a) std::printf(" %.8f", eps[a]);
                        std::printf("  J = %.8f  K = %.8f  max<a|b> = %.1e  E_phys = %.10f  dE = %+.3e  (%.1f s)\n",
                                    Jab, Kab, ovl, Etot - ezee, dE, sec);
                        std::fflush(stdout);
                    }
                    if (scf > 0 && std::abs(dE) < tolscf) break;
                }
                const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                if (world.rank() == 0) {
                    std::printf("final UHF high-spin N = %d  E_tot = %.10f   E_phys = E_tot - N B/2 (spin Zeeman incl.) = %.10f   wall = %.1f s   converged = %s (last dE %+.1e)\n",
                                N, Eold, Eold - 0.5 * N * std::abs(B), sec, std::abs(dElast) < tolscf ? "yes" : "NO", dElast);
                    if (!std::isnan(ref)) std::printf("reference = %.10f   E - ref = %+.3e\n", ref, Eold - 0.5 * N * std::abs(B) - ref);
                    std::fflush(stdout);
                }
                // Hellmann-Feynman (fixed gauge origin: only Vnuc and E_nn depend on the nuclei; MRA has no
                // nucleus-centred basis, so no Pulay terms):
                //   g_A = dE/dR_A = <rho| dVnuc/dR_A> + dE_nn/dR_A,  rho = sum_i |phi_i|^2
                // Breathing coordinate (any number of nuclei, shape and orientation fixed): R_A = c + (R/R0) d_A,
                //   d_A = R_A - c (c = centroid), R = |R_1 - R_0| (2 nuclei: the bond length);
                //   dE/dR = sum_A g_A . d_A / R   (2 nuclei: (g_2 - g_1).u / 2).
                // |F_other| = |G - (G.D/D.D) D| (3N-vectors G = (g_A), D = (d_A)): forces that would change
                //   the shape or rotate/translate the molecule (e.g. a torque turning a tilted bond).
                if (mol.natom() < 2) break;
                const int nat = int(mol.natom());
                real_function_3d rho = real_factory_3d(world);
                for (const auto& p : phi) rho = rho + abs_square(p);
                std::vector<coord_3d> g(nat), d(nat);
                coord_3d c(0.0);
                for (int A = 0; A < nat; ++A) c += mol.get_atom(A).get_coords() * (1.0 / nat);
                for (int A = 0; A < nat; ++A) {
                    d[A] = mol.get_atom(A).get_coords() - c;
                    for (int ax = 0; ax < 3; ++ax) {
                        const real_function_3d dv = real_factory_3d(world).functor(madchem::MolecularDerivativeFunctor(mol, A, ax)).truncate_on_project();
                        g[A][ax] = inner(rho, dv) + mol.nuclear_repulsion_derivative(A, ax);
                    }
                }
                const double R = (mol.get_atom(1).get_coords() - mol.get_atom(0).get_coords()).normf();
                double GD = 0.0, DD = 0.0;
                for (int A = 0; A < nat; ++A)
                    for (int ax = 0; ax < 3; ++ax) { GD += g[A][ax] * d[A][ax]; DD += d[A][ax] * d[A][ax]; }
                const double gR = GD / R;
                double other = 0.0;
                for (int A = 0; A < nat; ++A)
                    for (int ax = 0; ax < 3; ++ax) other += std::pow(g[A][ax] - GD / DD * d[A][ax], 2);
                other = std::sqrt(other);
                if (world.rank() == 0) {
                    std::printf("geo %2d  R = %.6f  E_phys = %.10f  dE/dR (HF) = %+.6e  |F_other| = %.2e  forces:",
                                geo, R, Eold - 0.5 * N * std::abs(B), gR, other);
                    for (int A = 0; A < nat; ++A) std::printf("  (%+.2e %+.2e %+.2e)", g[A][0], g[A][1], g[A][2]);
                    std::printf("\n");
                    std::fflush(stdout);
                }
                if (!ropt || std::abs(gR) < gtol) break;
                double Rn = R - (gR > 0 ? dr0 : -dr0);                 // first step: downhill by dr0
                if (geo > 0) {
                    const double kR = (gR - gprev) / (R - Rprev);      // secant curvature
                    Rn = kR > 0 ? R - gR / kR : R - (gR > 0 ? 2 * dr0 : -2 * dr0);
                }
                Rn = std::max(R - maxstep, std::min(R + maxstep, Rn));
                Rn = std::max(Rn, 0.02);
                Rprev = R;
                gprev = gR;
                if (std::abs(Rn - R) < rtol) break;
                // the nuclei must stay inside the region where the magnetic phase is kept (rcut - 0.70 off the axis)
                const double rcut_geo = get(args, "rcut", 0.0);
                double offmax = 0.0;
                for (int A = 0; A < nat; ++A) {
                    const coord_3d n = c + d[A] * (Rn / R);
                    offmax = std::max(offmax, std::sqrt(n[0] * n[0] + n[1] * n[1]));
                }
                if (rcut_geo > 0 && offmax > rcut_geo - 0.70) {
                    if (world.rank() == 0) std::printf("geo: STOP, R = %.6f would put a nucleus %.3f bohr off the B axis > rcut - 0.70\n", Rn, offmax);
                    break;
                }
                for (int A = 0; A < nat; ++A) {
                    const coord_3d n = c + d[A] * (Rn / R);
                    mol.set_atom_coords(A, n[0], n[1], n[2]);
                }
                V = real_factory_3d(world).functor(MolecularPotentialFunctor(mol)).truncate_on_project();
                V.truncate();
                if (world.rank() == 0) print("geo: new R =", Rn);
            }
        } else if (nelec == 2) {
            // closed-shell RHF, one orbital phi doubly occupied (singlet, no spin Zeeman term):
            //   F = h + v_phi,  v_phi = int |phi(r')|^2 / |r - r'| dr'   (2J - K = J for one orbital)
            //   eps = <h> + J11,  J11 = <phi|v_phi|phi>,  E_el = 2 <h> + J11
            // SCF: v_phi from the current phi, then a few Green-function iterations with V = Vnuc + v_phi
            const auto coul = CoulombOperator(world, get(args, "lo", 1e-6), 0.1 * thresh);
            const int inner_it = int(get(args, "inner", 3));
            const bool verbose = get(args, "verbose", 0) != 0;
            const int maxscf = int(get(args, "maxscf", 40));
            const double tolscf = get(args, "tolscf", 1e-6);
            double eps = get(args, "E0", -0.5);
            double Eold = 0.0;
            psi.scale(1.0 / psi.norm2());   // v_phi needs the normalized orbital from the first step
            auto t0 = std::chrono::steady_clock::now();
            if (world.rank() == 0) print("RHF: 2 electrons in one orbital, inner iterations =", inner_it);
            for (int scf = 0; scf < maxscf; ++scf) {
                auto ts = std::chrono::steady_clock::now();
                real_function_3d vphi = apply(coul, abs_square(psi));
                vphi.truncate();
                real_function_3d Veff = V + vphi;
                Veff.truncate();
                eps = solve(world, G, Veff, psi, eps, inner_it, tol, verbose);
                real_function_3d dens = abs_square(psi);
                vphi = apply(coul, dens);
                const double J11 = inner(dens, vphi);
                const double h1 = expect_H(psi, Stencil::abgv).E;
                const double Eel = 2.0 * h1 + J11;
                const double Etot = Eel + mol.nuclear_repulsion_energy();
                const double dE = Etot - Eold;
                Eold = Etot;
                const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
                if (world.rank() == 0) {
                    std::printf("  scf %2d  eps = %.10f  <h> = %.10f  J11 = %.10f  E_tot = %.10f  dE = %+.3e  (%.1f s)\n",
                                scf, eps, h1, J11, Etot, dE, sec);
                    std::fflush(stdout);
                }
                if (scf > 0 && std::abs(dE) < tolscf) break;
            }
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (world.rank() == 0) {
                std::printf("final RHF E_tot = %.10f   E_tot - B (2 x Landau zero point) = %.10f   wall = %.1f s\n",
                            Eold, Eold - std::abs(B), sec);
                if (!std::isnan(ref)) std::printf("reference = %.10f   E - ref = %+.3e\n", ref, Eold - ref);
                std::fflush(stdout);
            }
        } else {
        auto t0 = std::chrono::steady_clock::now();
        int iters = 0;
        const double E = solve(world, G, V, psi, get(args, "E0", -0.5), maxiter, tol, true, &iters);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (world.rank() == 0) {
            std::printf("final E = %.10f   binding B/2 - E = %.10f   iterations = %d   wall = %.1f s   converged = %s\n",
                        E, 0.5 * std::abs(B) - E, iters, sec, iters <= maxiter ? "yes" : "NO (maxiter reached)");
            if (mol.natom() > 1)
                std::printf("nuclear repulsion = %.10f   total E = %.10f   total binding B/2 - E_tot = %.10f\n",
                            mol.nuclear_repulsion_energy(), E + mol.nuclear_repulsion_energy(),
                            0.5 * std::abs(B) - E - mol.nuclear_repulsion_energy());
            if (!std::isnan(ref)) std::printf("reference = %.10f   E - ref = %+.3e\n", ref, E - ref);
            std::fflush(stdout);
        }
        // independent check: <psi|H|psi> of the converged psi with each derivative stencil
        for (Stencil t : {Stencil::abgv, Stencil::ble, Stencil::bspline}) {
            const Hparts h = expect_H(psi, t);
            const double Et = h.E;
            if (world.rank() == 0) {
                std::printf("<H> %-7s = %.10f  (kin %.8f para %.8f dia %.8f V %.8f)  <H> - E = %+.3e",
                            stencil_name(t), Et, h.kin, h.para, h.dia, h.pot, Et - E);
                if (!std::isnan(ref)) std::printf("  <H> - ref = %+.3e", Et - ref);
                std::printf("\n");
                std::fflush(stdout);
            }
        }
        }   // nelec == 1
    }
    catch (const madness::MadnessException& e) { print(e); status = 1; }
    catch (const std::exception& e) { print(e.what()); status = 1; }
    world.gop.fence();
    finalize();
    return status;
}
