#!/usr/bin/env python3
"""Writes all job files in ../jobs (run: python3 tools/make_jobs.py). Only a text generator: it computes
the nuclear coordinates and writes one line per calculation; all physics is done by magneticatom.

Conventions (all files)
  * B along z, gauge origin = origin = centroid of the nuclei. Atomic units (bohr, hartree).
  * Diatomics at angle theta to B:   +-(R/2) (sin theta, 0, cos theta).
  * H3 linear (along B):             (0,0,-d), (0,0,0), (0,0,+d).
  * H3 equilateral, side a, rho = a/sqrt(3):
      perp      plane perpendicular to B: (rho,0,0), (-rho/2, +-a/2, 0)
      xNN       perp triangle rotated by NN degrees about x: (x, y cos b, y sin b); x90: one edge parallel to B
      y90       perp triangle rotated by 90 degrees about y: one vertex on the B axis, opposite edge perpendicular to B
  * H3 bent, bond d, angle alpha = H-H-H, central atom on the B axis before centring:
      outer atoms d (sin g, 0, +-cos g), g = 90 - alpha/2; alpha = 180 is the linear chain along B,
      alpha = 60 is the x90 triangle with side d (same geometry up to a rotation about B: built-in check).
  * Structure optimization: jobs named ropt_* use Hellmann-Feynman forces along the breathing coordinate
    (all nuclei scaled about the centroid; shape and orientation fixed; R = distance nucleus 0 - nucleus 1).
    The log also prints |F_other|, the force outside this mode (shape change / rotation towards B):
    small |F_other| at the optimum = the symmetric structure is stationary. The R scans check the forces.
  * Every molecular file at B = 425.5 B0 uses the SAME numerical settings (MOL below), and the separated atoms
    are computed with exactly these settings (file 10), so dissociation energies are consistent.
  * Starts ('axial' only selects the start orbitals with; nothing is restricted in the SCF):
      mol   orbitals centred on the B axis with Lz = 0, -1, -2, ... (axial=1, ma/mb or mlist)
      atm   each atom in its own ground configuration, magnetically translated to its nucleus (axial=0, mloc)
    Where both are run, the lower converged E_phys is the ground state at that geometry.
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
JOBS = os.path.join(HERE, "..", "jobs")

B = 425.5
# one set of numerical settings for all molecules and their atoms at B = 425.5 B0 (= 1.0e12 G); only rcut differs:
#   rcut 0.75 when all nuclei are on the B axis, rcut 1.0 when nuclei are up to 0.30 bohr off the axis
#   (rule: rcut >= 0.70 + largest distance of a nucleus from the axis; 0.65 is converged for the H atom,
#   +0.05 for the more extended Lz = -2, -3 orbitals). File 10 computes both atoms at BOTH rcut values:
#   they must agree to < 1 uHa.
MOL0 = ("B=425.5 k=10 thresh=1e-6 eprec=1e-9 L=8 qextra=10 "
        "h=0.25 tmin=1e-8 tmax_factor=45 dphi=3 inner=3 tolscf=1e-5")


def MOL(rc):
    return f"{MOL0} rcut={rc} rtaper={rc / 10:g}"


RC_AXIS, RC_OFF = 0.75, 1.0


def fmt(x):
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def atoms(Z, pts, rc):
    off = max(math.hypot(p[0], p[1]) for p in pts)
    assert off <= rc - 0.70 + 1e-9, f"nucleus {off:.3f} bohr off the B axis: rcut {rc} too small"
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    cz = sum(p[2] for p in pts) / len(pts)
    assert abs(cx) < 1e-12 and abs(cy) < 1e-12 and abs(cz) < 1e-12, "centroid not at the origin"
    return '"atoms=' + ";".join(f"{Z},{fmt(p[0])},{fmt(p[1])},{fmt(p[2])}" for p in pts) + '"'


def centre(pts):
    c = [sum(p[i] for p in pts) / len(pts) for i in range(3)]
    return [tuple(p[i] - c[i] for i in range(3)) for p in pts]


def diatomic(R, th):
    t = math.radians(th)
    h = 0.5 * R
    return [(-h * math.sin(t), 0.0, -h * math.cos(t)), (h * math.sin(t), 0.0, h * math.cos(t))]


def triangle(a, orient):
    r = a / math.sqrt(3)
    p = [(r, 0.0, 0.0), (-r / 2, a / 2, 0.0), (-r / 2, -a / 2, 0.0)]
    if orient == "perp":
        return p
    if orient.startswith("x"):
        b = math.radians(float(orient[1:]))
        return [(x, y * math.cos(b), y * math.sin(b)) for x, y, z in p]
    if orient == "y90":
        return [(z, y, -x) for x, y, z in p]
    raise ValueError(orient)


def bent(d, alpha):
    g = math.radians(90 - alpha / 2)
    return centre([(d * math.sin(g), 0.0, -d * math.cos(g)), (0.0, 0.0, 0.0), (d * math.sin(g), 0.0, d * math.cos(g))])


def write(name, header, common, lines):
    path = os.path.join(JOBS, name + ".jobs")
    with open(path, "w") as f:
        for h in header:
            f.write("# " + h + "\n" if h else "#\n")
        f.write("COMMON " + common + "\n")
        w = max(len(n) for n, _ in lines)
        for n, a in lines:
            f.write(f"{n:<{w}}  {a}\n")
    print(f"{path}: {len(lines)} jobs")


# start orbitals and starting orbital energies (E0 < B/2 = 212.75; only a start, the SCF updates them)
H2_MOL = "nelec=2 state=triplet axial=1 ma=0 mb=-1 E0=203.5 E0b=207.7"
H2_ATM = "nelec=2 state=triplet axial=0 mloc=0 E0=207"
H3_MOL = "nelec=3 state=highspin axial=1 mlist=0,-1,-2 E0list=203,207,209"
H3_ATM = "nelec=3 state=highspin axial=0 mloc=0 E0=207"
HE_ATOM = "nelec=2 state=triplet axial=1 ma=0 mb=-1 E0=200.09 E0b=205.37"
HE2_MOL = "nelec=4 state=highspin axial=1 mlist=0,-1,-2,-3 E0list=201,205.5,207.5,209"
HE2_ATM = "nelec=4 state=highspin axial=0 mloc=0,-1 E0list=202,207,202,207"
OPT = "ropt=1 maxgeo=10 gtol=2e-3 rtol=1e-4 dr0=0.02 maxstep=0.05"

# ---------------------------------------------------------------- hydrogen atom: Fourier cutoff and settings
HATOM = "k=10 eprec=1e-9 tol=1e-7 maxiter=20 h=0.25 tmin=1e-8 tmax_factor=45 dphi=3"
write("01_H_B100_qmax", [
    "H atom, B = 100 B0: Fourier cutoff qmax = B (rcut + 3 rtaper)/2 + qextra = 97 + (qextra - 10) at rcut 1.34.",
    "Exact binding 3.789804236305 Ha (Kravchenko et al., PRA 54, 287 (1996), Table I).",
    "Summary: bash summarize.sh results/01_H_B100_qmax 3.789804236305"],
    f"B=100 E0=48.1 thresh=1e-7 L=8 rcut=1.34 rtaper=0.134 {HATOM} ref=46.210195763695",
    [(f"qextra{q}", f"qextra={q}") for q in (10, 30, 50, 80, 120, 160, 240, 320)])

write("02_H_B100_settings", [
    "H atom, B = 100 B0: every other numerical setting, one change per job, at qextra = 80 (cutoff error -8 uHa)",
    "and rcut 2.13, rtaper 0.4, dphi 1.5 at qextra = 10.",
    "Compare each job with qextra80 / qextra10 of file 01. Summary: bash summarize.sh results/02_H_B100_settings 3.789804236305"],
    f"B=100 E0=48.1 thresh=1e-7 L=8 rcut=1.34 rtaper=0.134 {HATOM} ref=46.210195763695",
    [("q80_k8", "qextra=80 k=8"), ("q80_k12", "qextra=80 k=12"),
     ("q80_th1e-6", "qextra=80 thresh=1e-6"), ("q80_th1e-8", "qextra=80 thresh=1e-8"),
     ("q80_L6", "qextra=80 L=6"), ("q80_L12", "qextra=80 L=12"),
     ("q10_rcut2.13", "qextra=10 rcut=2.13 rtaper=0.213"), ("q10_rtaper0.4", "qextra=10 rtaper=0.4"),
     ("q10_dphi1.5", "qextra=10 dphi=1.5")])

write("03_H_B425_rcut_L_thresh7", [
    "H atom, B = 425.5 B0: r_cut x box L at thresh 1e-7 (Florian's question: is the result independent of both?).",
    "Summary: bash summarize.sh results/03_H_B425_rcut_L_thresh7"],
    f"B=425.5 E0=207 thresh=1e-7 qextra=10 {HATOM}",
    [(f"rcut{rc}_L{L}", f"rcut={rc} rtaper={rc/10:g} L={L}") for rc in (1.0, 1.2) for L in (4, 8, 12)])

write("04_H_B425_qmax", [
    "H atom, B = 425.5 B0: Fourier cutoff at strong field, rcut 0.65: qmax = 190 + (qextra - 10).",
    "compare with files 03 and 10.",
    "Summary: bash summarize.sh results/04_H_B425_qmax"],
    f"B=425.5 E0=207 thresh=1e-7 L=4 rcut=0.65 rtaper=0.065 {HATOM}",
    [(f"qextra{q}", f"qextra={q}") for q in (10, 50, 100, 150, 200)])

write("05_H_B500_qmax", [
    "H atom, B = 500 B0: second exact reference, binding 6.257087674681 Ha (Kravchenko 1996, Table I).",
    "Error expected to fall like qmax^-3 (as at B = 100 B0, file 01).",
    "Summary: bash summarize.sh results/05_H_B500_qmax 6.257087674681"],
    f"B=500 E0=243.7 thresh=1e-7 L=6 rcut=0.6 rtaper=0.06 {HATOM} ref=243.742912325319",
    [(f"qextra{q}", f"qextra={q}") for q in (10, 50, 100, 150, 200)])

write("06_H_B1_qmax", [
    "H atom, B = 1 B0 (weak field: the Fourier cutoff is set by qextra alone, qmax ~ 9 + qextra).",
    "Exact binding 0.831168896733 Ha (Kravchenko 1996).",
    "Settings identical to the Bischoff H3 run (file 07), so this is also the H-atom reference of file 07.",
    "Summary: bash summarize.sh results/06_H_B1_qmax 0.831168896733"],
    "B=1 E0=-0.3 axial=0 k=9 thresh=1e-6 eprec=1e-8 L=20 maxiter=30 tol=1e-7 ref=-0.331168896733",
    [(f"qextra{q}", f"qextra={q}") for q in (10, 20, 40, 60)])

a_b = 1.445 / 0.529177210903
write("07_H3_B1_Bischoff", [
    "Cyclic H3 quartet (3 spins anti-parallel to B), B = 1 B0, equilateral triangle, side 144.5 pm = %.4f bohr," % a_b,
    "plane perpendicular to B. Same level of theory (UHF) as Bischoff, PRA 101, 053413 (2020), Table II:",
    "E = -2.50237 Eh (3 H atoms: -2.49349 Eh). Validation of the 3-electron code on a cyclic structure.",
    "Expensive: ~600 modes, L = 20, 3 electrons.",
    "Summary: bash summarize.sh results/07_H3_B1_Bischoff  ;  difference to Bischoff: grep 'E - ref' results/07_H3_B1_Bischoff/*.log"],
    "B=1 axial=0 k=9 thresh=1e-6 eprec=1e-8 L=20 nelec=3 state=highspin mloc=0 E0=-0.3 inner=3 maxscf=40 tolscf=1e-6 ref=-2.50237",
    [("H3_q20", "qextra=20 " + '"atoms=1,1.577089,0,0;1,-0.788545,1.365800,0;1,-0.788545,-1.365800,0"')])

# ---------------------------------------------------------------- B = 425.5 B0: atoms with the molecular settings
write("10_atoms_B425", [
    "Separated atoms at B = 425.5 B0 with EXACTLY the numerical settings of the molecule files 11-20, at both",
    "rcut values used there (0.75: nuclei on the B axis; 1.0: nuclei off the axis). The two must agree to < 1 uHa.",
    "E_phys (spin Zeeman included) of these runs is the zero of every dissociation energy:",
    "  D(H2) = 2 E(H) - E(H2),  D(H3) = 3 E(H) - E(H3),  D(He2) = 2 E(He) - E(He2)   (positive = bound).",
    "H: E_phys = -(binding).",
    "He: spin-polarized 1s0 2p-1 (Lz 0, -1), the He ground state above B ~ 0.75 B0.",
    "Summary: bash summarize.sh results/10_atoms_B425"],
    MOL0,
    [(f"{n}_rcut{rc}", f"{a} rcut={rc} rtaper={rc / 10:g}") for rc in (RC_AXIS, RC_OFF)
     for n, a in (("H_atom", "nelec=1 E0=207 tol=1e-7 maxiter=20"), ("He_atom", HE_ATOM + " maxscf=30"))])

# ---------------------------------------------------------------- H2
RS0 = (0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.32, 0.40, 0.50)
write("11_H2_theta0_B425", [
    "H2 spin triplet (Lz 0, -1 start), bond ALONG B: potential curve E(R) and a bond-length search (ropt).",
    "Lai, RMP 73, 629 (2001), Table I, B = 425.5 B0 (his B12 = 1, i.e. 1e12 G): Q2(inf) = 2 E(H) - E(H2) = 46 eV (accuracy 10 %).",
    "Summary: bash summarize.sh results/11_H2_theta0_B425 <2 E(H) from file 10>"],
    MOL(RC_AXIS) + " maxscf=30",
    [(f"R{fmt(R)}", f"{H2_MOL} {atoms(1, diatomic(R, 0), RC_AXIS)}") for R in RS0]
    + [("ropt_R0.24", f"{H2_MOL} {OPT} {atoms(1, diatomic(0.24, 0), RC_AXIS)}")])

RS = (0.20, 0.24, 0.28, 0.33, 0.40, 0.50, 0.60)
for th, fn in ((90, "12_H2_theta90_B425"), (60, "13_H2_theta60_B425"), (30, "14_H2_theta30_B425")):
    write(fn, [
        f"H2 spin triplet, bond at theta = {th} deg to B (bond in the x-z plane), potential curve E(R).",
        "Two starts per R: mol (Lz 0, -1 about the B axis) and atm (two magnetically translated H 1s).",
        "The lower converged E_phys is the ground state at that R. theta = 90: Lange et al., Science 337, 327",
        "(2012) found the H2 triplet bound perpendicular to B (paramagnetic bonding, B = 0.75-2.25 B0 there); no reference exists at this field.",
        "R up to 0.60 to see whether this orientation is bound at all.",
        f"Summary: bash summarize.sh results/{fn} <2 E(H) from file 10>"],
        MOL(RC_OFF) + " maxscf=30",
        [(f"R{fmt(R)}_{s}", f"{st} {atoms(1, diatomic(R, th), RC_OFF)}") for R in RS for s, st in (("mol", H2_MOL), ("atm", H2_ATM))]
        + [(f"ropt_R0.26_{s}", f"{st} {OPT} {atoms(1, diatomic(0.26, th), RC_OFF)}") for s, st in (("mol", H2_MOL), ("atm", H2_ATM))])

# ---------------------------------------------------------------- H3
write("15_H3_linear_B425", [
    "H3 quartet (3 spins anti-parallel to B), linear symmetric chain ALONG B, spacing d, start Lz 0, -1, -2.",
    "Nader, Lopez Vieyra, Turbiner, PRA 100, 012508 (2019), Table I, B = 425.5 B0 (their 1000 x 1e9 G = 1e12 G):",
    "E_T = -20.90 Ha, R_eq = 0.219 bohr. Both theirs (correlated trial function) and ours (UHF) are upper bounds to the",
    "exact energy, so no strict order: expect agreement to roughly 0.1 Ha and in R_eq to a few %.",
    "Summary: bash summarize.sh results/15_H3_linear_B425 <3 E(H) from file 10>"],
    MOL(RC_AXIS) + " maxscf=40",
    [(f"d{fmt(d)}", f"{H3_MOL} {atoms(1, [(0, 0, -d), (0, 0, 0), (0, 0, d)], RC_AXIS)}")
     for d in (0.17, 0.19, 0.21, 0.22, 0.23, 0.25, 0.28, 0.32)]
    + [("ropt_d0.22", f"{H3_MOL} {OPT} {atoms(1, [(0, 0, -0.22), (0, 0, 0), (0, 0, 0.22)], RC_AXIS)}")])

AS = (0.22, 0.26, 0.30, 0.36, 0.45)
write("16_H3_triangle_perp_B425", [
    "Cyclic H3 quartet, equilateral triangle of side a, plane PERPENDICULAR to B (C3 axis along B).",
    "Two starts: mol (Lz 0, -1, -2 about the B axis) and atm (three magnetically translated H 1s).",
    "Compare with the linear chain (file 15) and with H2 + H (files 10, 11).",
    "Summary: bash summarize.sh results/16_H3_triangle_perp_B425 <3 E(H) from file 10>"],
    MOL(RC_OFF) + " maxscf=40",
    [(f"a{fmt(a)}_{s}", f"{st} {atoms(1, triangle(a, 'perp'), RC_OFF)}") for a in AS for s, st in (("mol", H3_MOL), ("atm", H3_ATM))]
    + [(f"ropt_a0.30_{s}", f"{st} {OPT} {atoms(1, triangle(0.30, 'perp'), RC_OFF)}") for s, st in (("mol", H3_MOL), ("atm", H3_ATM))])

write("17_H3_triangle_tilt_B425", [
    "Cyclic H3 quartet, equilateral triangle of side a, plane TILTED against B (atomic start, mloc=0):",
    "  x30, x60  perp triangle rotated by 30, 60 deg about x;  x90  plane contains B, one edge parallel to B;",
    "  y90  plane contains B, one vertex on the B axis (median parallel to B). Angle 0 = file 16.",
    "Check: x90 with a = 0.22 is the same geometry as bent alpha = 60 with d = 0.22 (file 18).",
    "Summary: bash summarize.sh results/17_H3_triangle_tilt_B425 <3 E(H) from file 10>"],
    MOL(RC_OFF) + " maxscf=40",
    [(f"{o}_a{fmt(a)}", f"{H3_ATM} {atoms(1, triangle(a, o), RC_OFF)}") for o in ("x30", "x60", "x90", "y90") for a in AS]
    + [(f"ropt_{o}_a0.30", f"{H3_ATM} {OPT} {atoms(1, triangle(0.30, o), RC_OFF)}") for o in ("x30", "x60", "x90", "y90")])

write("18_H3_bent_B425", [
    "H3 quartet, bending the linear chain: bond d = 0.22 (Nader R_eq), angle alpha = H-H-H in the x-z plane,",
    "central atom on the B axis before centring; alpha = 180 linear along B (file 15), alpha = 60 the x90 triangle.",
    "Two starts per angle (mol: Lz 0, -1, -2; atm: three translated 1s); the lower E_phys counts.",
    "Summary: bash summarize.sh results/18_H3_bent_B425 <3 E(H) from file 10>"],
    MOL(RC_OFF) + " maxscf=40",
    [(f"alpha{al}_{s}", f"{st} {atoms(1, bent(0.22, al), RC_OFF)}") for al in (180, 150, 120, 90, 60)
     for s, st in (("mol", H3_MOL), ("atm", H3_ATM))])

# ---------------------------------------------------------------- He2
write("19_He2_theta0_B425", [
    "He2, 4 spins anti-parallel to B, bond ALONG B: potential curve E(R) + bond-length search.",
    "mol start: Lz 0, -1, -2, -3 about the axis: all 4 electrons tightly bound (nu = 0), as found by Medin & Lai 2006, Table III.",
    "atm start at 3 R (two He atoms in 1s0 2p-1): checks that both starts reach the same state.",
    "References at B = 425.5 B0 (bond along B), energies per atom:",
    "  Hartree-Fock (Lai, RMP 73, 629 (2001), as quoted by Medin & Lai, PRA 74, 062507 (2006)): He -575.5 eV, He2 -601.2 eV",
    "    -> D(He2) = 2 x 25.7 = 51.4 eV. (Lai's HF is the adiabatic approximation, so exact HF is somewhat lower.)",
    "  DFT (Medin & Lai 2006, Table III): He -603.5 eV, He2 -641.2 eV, R_eq = 0.25 bohr (DFT overbinds by ~10 %).",
    "Summary: bash summarize.sh results/19_He2_theta0_B425 <2 E(He) from file 10>"],
    MOL(RC_AXIS) + " maxscf=40",
    [(f"R{fmt(R)}", f"{HE2_MOL} {atoms(2, diatomic(R, 0), RC_AXIS)}") for R in (0.22, 0.26, 0.30, 0.34, 0.38, 0.45, 0.55, 0.70, 0.90)]
    + [(f"R{fmt(R)}_atm", f"{HE2_ATM} {atoms(2, diatomic(R, 0), RC_AXIS)}") for R in (0.30, 0.45, 0.70)]
    + [("ropt_R0.26", f"{HE2_MOL} {OPT} {atoms(2, diatomic(0.26, 0), RC_AXIS)}")])

write("20_He2_angles_B425", [
    "He2 (4 spins anti-parallel to B), bond at theta = 45 and 90 deg to B, two starts (mol, atm) per point.",
    "Most expensive file (4 electrons, rcut 1.0): run after file 19; add R values near its R_eq if needed.",
    "No reference at this field. Austad, Borgoo, Tellgren, Helgaker, PCCP (2020) (doi 10.1039/D0CP03259J), FCI, Table V,",
    "B = 1 B0 only: quintet ground state UNBOUND along B (R_eq = inf) but bound perpendicular (R_eq 3.62, D 0.0198 Ha).",
    "At 425.5 B0 He2 is strongly bound along B (file 19 refs): whether the perpendicular preference survives is open.",
    "Summary: bash summarize.sh results/20_He2_angles_B425 <2 E(He) from file 10>"],
    MOL(RC_OFF) + " maxscf=40",
    [(f"th{th}_R{fmt(R)}_{s}", f"{st} {atoms(2, diatomic(R, th), RC_OFF)}")
     for th, Rs in ((90, (0.30, 0.40, 0.55)), (45, (0.30, 0.40, 0.55)))
     for R in Rs for s, st in (("mol", HE2_MOL), ("atm", HE2_ATM))]
    + [(f"ropt_th{th}_R0.40_{s}", f"{st} {OPT} {atoms(2, diatomic(0.40, th), RC_OFF)}")
       for th in (90, 45) for s, st in (("mol", HE2_MOL), ("atm", HE2_ATM))])
