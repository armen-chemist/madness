/// \file magnetic_green.h
/// \brief Green's function of the kinetic operator of an electron in a uniform magnetic field
///        (B along z, symmetric gauge), applied by the Fourier approach.
///
///   T_B = 1/2 (p + A)^2 = -1/2 nabla^2 + (B/2) Lz + (B^2/8)(x^2 + y^2),   A = (B/2)(-y, x, 0)
///
/// inf spec T_B = |B|/2 (lowest Landau level), and for E < |B|/2
///
///   G_B(E) = (T_B - E)^{-1} = int_0^inf dtau e^{E tau} e^{-tau T_B}.
///
/// T_B = T_perp(x, y) - 1/2 d^2/dz^2 with commuting terms, so e^{-tau T_B} = U_perp(tau) K_z(tau),
/// K_z(d) = (2 pi tau)^{-1/2} exp(-d^2/(2 tau)).  U_perp is not a convolution (its phase couples
/// output and input coordinates), but it separates exactly as
///
///   U_perp = e^{iBxy/2} (1/2pi) int dq e^{iq(y-y')} g(x-X_q) N(x-x') g(x'-X_q) e^{-iBx'y'/2},
///
///   X_q = -q/B,  g(u) = exp(-m u^2),  N(d) = sqrt(c/pi) exp(-c d^2),
///   m = (B/2) tanh(B tau/2),  c = B / (2 sinh(B tau)),
///
/// with the mask written as m (x - X_q)^2 = m x^2 + b q x + d q^2, b = tanh(B tau/2),
/// d = tanh(B tau/2)/(2B) (-> tau/4 as B -> 0), so that every factor is finite for any B.
/// Applied to f:
///
///   1. ft = e^{-iBxy/2} f                                      (3D multiply)
///   2. F_q(x, z) = int dy e^{-iqy} ft(x, y, z)                 (project_out, per q)
///   3. Psi_q = g_q(x) * [ (N_x K_z) * (g_q(x) F_q) ]           (x-mask, 2D separable conv., x-mask)
///   4. result = e^{iBxy/2} (1/P) sum_q e^{iqy} Psi_q(x, z)     (hartree_product)
///
/// The q-integral is a sum over q_j = 2 pi j / P (y treated as periodic with period P).  Steps 1, 2
/// and 4 are done once; only step 3 is summed over the proper-time nodes, a trapezoid rule in
/// u = ln tau from tmin to tmax_factor / (|B|/2 - E), with [0, tmin] added as tmin f.
/// MADNESS's project_out integrates over the last block of dimensions, so steps 2 and 4 work in
/// the order (x, z, y): mapdim with {0, 2, 1} (a swap, its own inverse) before and after.
///
/// G_B(E) needs only E < |B|/2, for ANY orientation of a molecule; no (B/2) Lz term is left in the
/// potential.  The cost is the number of q-modes, ~ (|B| L/2 + qextra) P / pi, times the nodes.
///
/// Requirements: FunctionDefaults<3> holds a cubic cell; the constructor sets FunctionDefaults<1>
/// and <2> (cell, k, thresh) to match it.

#ifndef MADNESS_MRA_MAGNETIC_GREEN_H__INCLUDED
#define MADNESS_MRA_MAGNETIC_GREEN_H__INCLUDED

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <cmath>
#include <memory>
#include <vector>

namespace madness {

struct MagneticGreenParams {
    double h = 0.25;            ///< proper-time step in u = ln(tau) ...
    double tmin = 1e-8;         ///<   ... from tmin ([0, tmin] is added as tmin * f) ...
    double tmax_factor = 45.0;  ///<   ... to tmax_factor / (|B|/2 - E)
    double period = 0.0;        ///< y-period P; 0 -> min(2 x cell width, 4 rext)
    double qextra = 10.0;       ///< |q| <= |B| L/2 + qextra (wavenumbers of f itself)
    double rcut = 0.0;          ///< phase e^{+-iBxy/2} cut off (erfc taper) at rho > rcut; 0 -> 1.5 sqrt(8 ln(1/thresh)/|B|)
    double rtaper = 0.0;        ///< width of that erfc taper; 0 -> min(0.5, rcut/4)
    double dphi = 3.0;          ///< refine f until the phase B x y/2 varies by <= dphi per box (before the multiply)
    double fq_tol = 0.0;        ///< drop modes with ||F_q|| < fq_tol*||f||; 0 -> 0.01*thresh
    long qbatch = 64;           ///< q-modes processed together (one fence per step and batch); 0 -> all
    bool verbose = false;
};

class MagneticGreen {
public:
    using Params = MagneticGreenParams;
    using operator2dT = SeparatedConvolution<double, 2>;

    MagneticGreen(World& world, double B, const Params& p = Params())
        : world_(world), B_(B), p_(p) {
        const Tensor<double> cell = FunctionDefaults<3>::get_cell();
        lo_ = cell(0, 0);
        hi_ = cell(0, 1);
        for (int d = 1; d < 3; ++d)
            MADNESS_CHECK_THROW(cell(d, 0) == lo_ && cell(d, 1) == hi_,
                                "MagneticGreen: 3D cell must be cubic");
        FunctionDefaults<1>::set_cubic_cell(lo_, hi_);
        FunctionDefaults<2>::set_cubic_cell(lo_, hi_);
        FunctionDefaults<1>::set_k(FunctionDefaults<3>::get_k());
        FunctionDefaults<2>::set_k(FunctionDefaults<3>::get_k());
        FunctionDefaults<1>::set_thresh(FunctionDefaults<3>::get_thresh());
        FunctionDefaults<2>::set_thresh(FunctionDefaults<3>::get_thresh());

        // the phase is stored only on the cylinder rho < rcut (+ erfc taper of width 0.5): f and G f
        // decay as e^{-|B| rho^2/4} across B, and a full-cell 3D phase refines everywhere (memory)
        const double thr = FunctionDefaults<3>::get_thresh();
        // (the default rcut assumes the nuclei near the z-axis: pass rcut for off-axis nuclei)
        rcut_ = p_.rcut > 0 ? p_.rcut
              : (B_ == 0.0 ? 1e10 : 1.5 * std::sqrt(8.0 * std::log(1.0 / thr) / std::abs(B_)));
        rtaper_ = p_.rtaper > 0 ? p_.rtaper : std::min(0.5, rcut_ / 4);
        // ft = e^{-iBxy/2} f lives on |x|, |y| < rext: y-wavenumbers up to |B| rext/2 + those of f,
        // and the period P = 4 rext keeps the y-images of ft and of G ft apart (P <= 2 x cell width)
        const double rext = rcut_ + 3 * rtaper_;
        const double L = std::min(std::max(std::abs(lo_), std::abs(hi_)), rext);
        period_ = p_.period > 0 ? p_.period : std::min(2.0 * (hi_ - lo_), 4.0 * rext);
        const double qmax = std::abs(B_) * L / 2 + p_.qextra;
        const long J = long(std::ceil(qmax * period_ / constants::pi / 2));
        for (long j = -J; j <= J; ++j) q_.push_back(2 * constants::pi * j / period_);

        // step 4 uses e^{iqy} windowed to |y| < rcut + 9 rtaper (1 to 1e-10 wherever the final taper
        // is non-negligible, 0 well before the y-images at |y| ~ P - rext): every term of the q-sum is
        // then local in y, as the sum is (otherwise each term fills the cell at the q-wavelength)
        const double yw = rcut_ + 9 * rtaper_, rt = rtaper_;
        for (double q : q_) {
            // e^{-iqy} of step 2, kept redundant (the state project_out reads) for the whole run
            complex_function_1d pw = complex_factory_1d(world_).functor(
                [q](const coord_1d& r) { return std::exp(double_complex(0.0, q * r[0])); });
            planewave_conj_.push_back(conj(pw));
            planewave_conj_.back().change_tree_state(redundant);
            planewave_win_.push_back(complex_factory_1d(world_).functor([q, yw, rt](const coord_1d& r) {
                return 0.5 * std::erfc((std::abs(r[0]) - yw) / rt) * std::exp(double_complex(0.0, q * r[0]));
            }));
        }
        if (p_.verbose && world_.rank() == 0)
            print("MagneticGreen: B =", B_, " cell [", lo_, ",", hi_, "]  period =", period_,
                  " modes =", q_.size(), " qmax =", q_.back(), " rcut =", rcut_, " rtaper =", rtaper_);
    }

    /// inf spec T_B = |B|/2
    double bottom() const { return 0.5 * std::abs(B_); }

    /// G_B(E) f = (T_B - E)^{-1} f,  E < |B|/2
    complex_function_3d green(const complex_function_3d& f, double E) const {
        MADNESS_CHECK_THROW(E < bottom(), "MagneticGreen::green: require E < |B|/2");
        const double gap = bottom() - E;
        const double umin = std::log(p_.tmin), h = p_.h;
        const double umax = std::log(p_.tmax_factor / gap);
        const long npt = long(std::floor((umax - umin) / h + 0.5)) + 1;
        std::vector<Node> nodes;
        for (long i = 0; i < npt; ++i) {
            const double tau = std::exp(umin + i * h);
            Node nd = node(tau);
            // weight of the NORMALIZED kernel sqrt(c/pi) sqrt(1/(2 pi tau)), in log form
            nd.log_w = std::log(h * tau) + E * tau + 0.5 * (nd.log_c - std::log(constants::pi))
                     + 0.5 * (-std::log(2 * tau) - std::log(constants::pi));
            if (i == 0) nd.log_w += std::log(0.5);       // trapezoid endpoint
            if (nd.log_w > -60.0) nodes.push_back(nd);
        }
        if (p_.verbose && world_.rank() == 0)
            print("MagneticGreen::green: E =", E, " nodes =", nodes.size(), "of", npt);
        return apply(f, nodes, std::exp(umin));
    }

    const std::vector<double>& qgrid() const { return q_; }

private:
    struct Node {
        double tau, m, b, d, log_c, log_w;
    };

    static double tanh_over_x(double x) { return x == 0.0 ? 1.0 : std::tanh(x) / x; }

    /// mask exponents m, b, d and log c of node tau (finite at B = 0)
    Node node(double tau) const {
        const double x = 0.5 * B_ * tau, bt = std::abs(B_) * tau;
        Node nd;
        nd.tau = tau;
        nd.m = 0.5 * B_ * std::tanh(x);
        nd.b = std::tanh(x);
        nd.d = 0.25 * tau * tanh_over_x(x);
        nd.log_c = -std::log(2 * tau)
                 + (bt == 0.0 ? 0.0 : std::log(2 * bt) - bt - std::log(-std::expm1(-2.0 * bt)));
        nd.log_w = 0.0;
        return nd;
    }

    /// unnormalized exp(-c dx^2 - dz^2/(2 tau)) on the (x, z) plane; an exponent below the floor
    /// (c underflows at large B tau) is raised to it: the kernel is then 1 on the cell to ~cfloor
    std::shared_ptr<operator2dT> op2d(const Node& nd) const {
        const Tensor<double>& width = FunctionDefaults<2>::get_cell_width();
        const int k = FunctionDefaults<2>::get_k();
        const auto lattice = FunctionDefaults<2>::get_bc().lattice_range();
        const double la[2] = {nd.log_c, -std::log(2 * nd.tau)};
        ConvolutionND<double, 2> term;
        double fac = 1.0;
        for (int d = 0; d < 2; ++d) {
            const double w2 = width[d] * width[d];
            const double a = std::exp(std::max(la[d] + std::log(w2), std::log(1e-10)));
            term.setop(d, GaussianConvolution1DCache<double>::get(k, a, 0, lattice[d]));
            fac *= std::sqrt(constants::pi / a) * width[d];
        }
        term.setfac(fac);
        return std::make_shared<operator2dT>(world_, std::vector<ConvolutionND<double, 2>>{term}, k);
    }

    /// refine a box (one level) if it meets rho < rext, its phase B x y/2 varies by more than dphi
    /// across it, and its coefficients are not negligible
    struct PhaseRefineOp {
        double Bh = 0, lo = 0, width = 0, rext = 0, dphi = 0, tol = 0;
        template <typename implT, typename nodeT>
        bool operator()(implT*, const Key<3>& key, const nodeT& t) const {
            const double h = width / double(1L << key.level());
            const double x0 = lo + key.translation()[0] * h, y0 = lo + key.translation()[1] * h;
            auto nearest = [h](double a) { return (a <= 0.0 && a + h >= 0.0) ? 0.0 : std::min(std::abs(a), std::abs(a + h)); };
            const double rx = nearest(x0), ry = nearest(y0);
            if (rx * rx + ry * ry > rext * rext) return false;
            const double mx = std::max(std::abs(x0), std::abs(x0 + h)), my = std::max(std::abs(y0), std::abs(y0 + h));
            if (std::abs(Bh) * (mx + my) * h <= dphi) return false;
            return t.coeff().normf() > tol;
        }
        template <typename Archive> void serialize(Archive& ar) { ar & Bh & lo & width & rext & dphi & tol; }
    };

    /// values at the quadrature points times erfc taper(rho) e^{i Bh x y}
    struct PhaseMulOp {
        double Bh = 0, lo = 0, width = 0, rc = 0, rt = 0;
        int k = 0;
        void operator()(const Key<3>& key, Tensor<double_complex>& v) const {
            const Tensor<double>& qx = FunctionCommonData<double, 3>::get(k).quad_x;
            const double h = width / double(1L << key.level());
            const auto& l = key.translation();
            for (long i = 0; i < v.dim(0); ++i) {
                const double x = lo + (l[0] + qx(i)) * h;
                for (long j = 0; j < v.dim(1); ++j) {
                    const double y = lo + (l[1] + qx(j)) * h;
                    const double rho = std::sqrt(x * x + y * y);
                    const double_complex fac = 0.5 * std::erfc((rho - rc) / rt) * std::exp(double_complex(0.0, Bh * x * y));
                    for (long m = 0; m < v.dim(2); ++m) v(i, j, m) *= fac;
                }
            }
        }
        template <typename Archive> void serialize(Archive& ar) { ar & Bh & lo & width & rc & rt & k; }
    };

    /// f <- taper(rho) e^{sign i B x y/2} f, pointwise on f's own tree after refining it where the
    /// phase is not resolved (no 3D phase function is stored: at large B it would refine everywhere)
    void multiply_phase(complex_function_3d& f, double sign) const {
        PhaseRefineOp rop;
        rop.Bh = 0.5 * B_;
        rop.lo = lo_;
        rop.width = hi_ - lo_;
        rop.rext = rcut_ + 3 * rtaper_;
        rop.dphi = p_.dphi;
        rop.tol = 0.1 * FunctionDefaults<3>::get_thresh();
        for (int pass = 0; pass < 12; ++pass) {
            const std::size_t before = f.tree_size();
            f.refine_general(rop);
            if (f.tree_size() == before) break;
        }
        PhaseMulOp mop;
        mop.Bh = sign * 0.5 * B_;
        mop.lo = lo_;
        mop.width = hi_ - lo_;
        mop.rc = rcut_;
        mop.rt = rtaper_;
        mop.k = f.k();
        f.unaryop(mop);
        f.truncate();
    }

    /// step 2 for the modes j0..j1-1: F_q(x, z) = int dy e^{-iqy} fxzy(x, z, y), as Function::project_out
    /// (fxzy reconstructed, e^{-iqy} redundant), but with one fence for all modes instead of three per mode
    std::vector<complex_function_2d> project_out_y(const complex_function_3d& fxzy, std::size_t j0, std::size_t j1) const {
        std::vector<complex_function_2d> F;
        for (std::size_t j = j0; j < j1; ++j)
            F.push_back(FunctionFactory<double_complex, 2>(world_).k(planewave_conj_[j].k())
                        .thresh(planewave_conj_[j].thresh()).nofence());
        world_.gop.fence();
        for (std::size_t j = j0; j < j1; ++j)
            fxzy.get_impl()->project_out(F[j - j0].get_impl().get(), planewave_conj_[j].get_impl().get(), 1, false);
        world_.gop.fence();
        for (auto& r : F) {
            r.get_impl()->trickle_down(false);
            r.get_impl()->set_tree_state(reconstructed);
        }
        world_.gop.fence();
        return F;
    }

    complex_function_3d apply(const complex_function_3d& f, const std::vector<Node>& nodes, double s0) const {
        const double thresh = FunctionDefaults<3>::get_thresh();
        const double fq_tol = (p_.fq_tol > 0 ? p_.fq_tol : 0.01 * thresh) * f.norm2();
        std::vector<std::shared_ptr<operator2dT>> ops;
        for (const Node& nd : nodes) ops.push_back(op2d(nd));

        // step 1, and (x, y, z) -> (x, z, y)
        complex_function_3d ft = copy(f);
        multiply_phase(ft, -1.0);
        complex_function_3d fxzy = mapdim(ft, std::vector<long>{0, 2, 1});
        fxzy.reconstruct();

        // the modes are processed in batches: every step below is one vector operation over the batch
        // with a single fence (one mode at a time leaves most threads idle: a 2D function is too small
        // to fill the task pool); per mode the arithmetic is the same as one mode at a time
        std::vector<complex_function_2d> left;
        std::vector<complex_function_1d> right;
        long napplied = 0, nskip = 0;
        const double wall0 = wall_time();
        double wall_proj = 0.0;
        const std::size_t nq = q_.size();
        const std::size_t batch = p_.qbatch > 0 ? std::size_t(p_.qbatch) : nq;
        for (std::size_t j0 = 0; j0 < nq; j0 += batch) {
            const std::size_t j1 = std::min(nq, j0 + batch);
            // step 2: F_q(x, z) = int dy e^{-iqy} ft   (project_out does not conjugate)
            const double wp = wall_time();
            std::vector<complex_function_2d> Fall = project_out_y(fxzy, j0, j1);
            const std::vector<double> nall = norm2s(world_, Fall);
            wall_proj += wall_time() - wp;
            std::vector<std::size_t> jq;
            std::vector<complex_function_2d> F;
            std::vector<double> nF;
            for (std::size_t i = 0; i < Fall.size(); ++i) {
                if (nall[i] < fq_tol) continue;
                jq.push_back(j0 + i);
                F.push_back(Fall[i]);
                nF.push_back(nall[i]);
            }
            Fall.clear();
            if (F.empty()) continue;
            // step 3
            std::vector<complex_function_2d> Psi(F.size());
            for (std::size_t i = 0; i < F.size(); ++i) Psi[i] = mul(s0, F[i], false);
            world_.gop.fence();
            compress(world_, Psi);
            for (std::size_t k = 0; k < nodes.size(); ++k) {
                const double m = nodes[k].m, b = nodes[k].b, d = nodes[k].d;
                std::vector<real_function_2d> g;
                std::vector<complex_function_2d> Fk, Pk;
                for (std::size_t i = 0; i < F.size(); ++i) {
                    const double q = q_[jq[i]];
                    auto expo = [m, b, d, q](double x) { return m * x * x + b * q * x + d * q * q; };
                    double emin = std::min(expo(lo_), expo(hi_));
                    if (m > 0) {
                        const double xs = -b * q / (2 * m);
                        if (xs > lo_ && xs < hi_) emin = std::min(emin, expo(xs));
                    }
                    // ||g K g F|| <= max|g|^2 ||K||_1 ||F||, ||K||_1 = sqrt(pi/c) sqrt(2 pi tau)
                    const double log_bound = nodes[k].log_w + 0.5 * (std::log(constants::pi) - nodes[k].log_c)
                                           + 0.5 * std::log(2 * constants::pi * nodes[k].tau) - 2.0 * emin + std::log(nF[i]);
                    if (log_bound < std::log(fq_tol / nodes.size())) { ++nskip; continue; }
                    ++napplied;
                    g.push_back(real_factory_2d(world_).functor(
                        [expo](const coord_2d& r) { return std::exp(-expo(r[0])); }).nofence());
                    Fk.push_back(F[i]);
                    Pk.push_back(Psi[i]);     // shares the implementation: gaxpy below updates Psi[i]
                }
                if (g.empty()) continue;
                world_.gop.fence();
                std::vector<complex_function_2d> t = mul(world_, g, Fk);
                t = madness::apply(world_, *ops[k], t);
                t = mul(world_, g, t);
                compress(world_, t);
                gaxpy(world_, 1.0, Pk, std::exp(nodes[k].log_w), t);
            }
            truncate(world_, Psi);
            for (std::size_t i = 0; i < Psi.size(); ++i) {
                left.push_back(Psi[i]);
                right.push_back(planewave_win_[jq[i]]);
            }
        }
        if (p_.verbose && world_.rank() == 0)
            print("MagneticGreen: kept", left.size(), "of", q_.size(), "modes, applied", napplied,
                  "(node, q) terms, skipped", nskip, " wall: project", wall_proj, "s, modes", wall_time() - wall0, "s");
        // step 4, and (x, z, y) -> (x, y, z)
        complex_function_3d rxzy = hartree_product(left, right);
        complex_function_3d result = mapdim(rxzy, std::vector<long>{0, 2, 1});
        result.scale(1.0 / period_);
        multiply_phase(result, +1.0);
        return result;
    }

    World& world_;
    double B_;
    Params p_;
    double lo_, hi_, period_, rcut_, rtaper_;
    std::vector<double> q_;
    std::vector<complex_function_1d> planewave_conj_, planewave_win_;
};

}  // namespace madness

#endif  // MADNESS_MRA_MAGNETIC_GREEN_H__INCLUDED
