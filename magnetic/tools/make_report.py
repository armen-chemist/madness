#!/usr/bin/env python3
"""Progress report (PDF) of the magnetic-field campaign: energies are read from results/*/*.log.
Usage:  ~/miniforge3/envs/report/bin/python tools/make_report.py [out.pdf]
"""
import datetime, glob, os, re, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, "results")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "docs", "report_%s.pdf" % datetime.date.today())
WORK = os.path.join(HERE, "docs", "report_figs")
HA_EV = 27.211386

# ------------------------------------------------------------------ data from the logs
def e_phys(log):
    """E_phys of a finished log (one-electron runs: E_phys = -binding), else None."""
    txt = open(log).read()
    m = re.findall(r"^final UHF.*incl\.\) = *(-?[0-9.]+)", txt, re.M)
    if m:
        return float(m[-1])
    m = re.findall(r"^final E .*B/2 - E = *(-?[0-9.]+)", txt, re.M)
    return -float(m[-1]) if m else None

def ropt(log):
    """(R, E_phys, dE/dR) of the last geometry step of an optimization log."""
    rows = re.findall(r"^geo +\d+ +R = *([0-9.]+) +E_phys = *(-?[0-9.]+) +dE/dR \(HF\) = *([-+0-9.e]+)",
                      open(log).read(), re.M)
    return tuple(float(x) for x in rows[-1]) if rows else None

def scan(folder, prefix):
    pts = []
    for log in glob.glob(os.path.join(RES, folder, prefix + "*.log")):
        name = os.path.basename(log)[:-4]
        if name.startswith("ropt"):
            continue
        e = e_phys(log)
        if e is not None:
            pts.append((float(name[len(prefix):]), e))
    return sorted(pts)

EH = e_phys(os.path.join(RES, "10_atoms_B425", "H_atom_rcut0.75.log"))
EH10 = e_phys(os.path.join(RES, "10_atoms_B425", "H_atom_rcut1.0.log"))
EHe = e_phys(os.path.join(RES, "10_atoms_B425", "He_atom_rcut0.75.log"))
EHe10 = e_phys(os.path.join(RES, "10_atoms_B425", "He_atom_rcut1.0.log"))
h2 = scan("11_H2_theta0_B425", "R")
h2opt = ropt(os.path.join(RES, "11opt_H2_theta0_B425", "ropt_R0.24.log"))
h3 = scan("15_H3_linear_B425", "d")
h3opt = ropt(os.path.join(RES, "15opt_H3_linear_B425", "ropt_d0.22.log"))
bis = e_phys(os.path.join(RES, "07_H3_B1_Bischoff", "H3_q20.log"))
bis_geo = ropt(os.path.join(RES, "07_H3_B1_Bischoff", "H3_q20.log"))

# exact H binding (Kravchenko, Liberman & Johansson 1996)
EXACT = {1: 0.831168896733, 10: 1.747797165, 100: 3.789804236305, 500: 6.257087674681, 1000: 7.662423245}
hlit = []
for folder, b in (("21_H_B1_lit", 1), ("22_H_B10_lit", 10), ("23_H_B100_lit", 100),
                  ("24_H_B500_lit", 500), ("25_H_B1000_lit", 1000)):
    logs = glob.glob(os.path.join(RES, folder, "*.log"))
    if logs and e_phys(logs[0]) is not None:
        hlit.append((b, -e_phys(logs[0]), os.path.basename(logs[0])[:-4]))
EH_B1 = -hlit[0][1]          # our H atom at B = 1 (run 21, same 631 modes as run 07)

# ------------------------------------------------------------------ figures
os.makedirs(WORK, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 2.9))
R = [p[0] for p in h2]; D = [(2 * EH - p[1]) * HA_EV for p in h2]
a1.plot(R, D, "o-", color="#2a6fb0", ms=4, label="this work, scan (UHF)")
a1.plot([h2opt[0]], [(2 * EH - h2opt[1]) * HA_EV], "*", color="#c0392b", ms=11, label="this work, optimized")
a1.axhspan(46 * 0.9, 46 * 1.1, color="0.85", zorder=0, label="Lai 2001: 46 eV ± 10 %")
a1.plot([0.243], [(2 * EH + 13.576) * HA_EV], "s", mfc="none", color="k", ms=6, label="Nader et al. 2020")
a1.set_xlabel("R (bohr)"); a1.set_ylabel("D = 2E(H) − E(H$_2$)  (eV)")
a1.set_title("H$_2$ triplet, bond ∥ B, B = 425 a.u.", fontsize=9); a1.legend(fontsize=7, frameon=False)
d = [p[0] for p in h3]; D3 = [(3 * EH - p[1]) * HA_EV for p in h3]
a2.plot(d, D3, "o-", color="#2a6fb0", ms=4, label="this work, scan (UHF)")
a2.plot([h3opt[0]], [(3 * EH - h3opt[1]) * HA_EV], "*", color="#c0392b", ms=11, label="this work, optimized")
a2.plot([0.219], [(3 * EH + 20.90) * HA_EV], "s", mfc="none", color="k", ms=6, label="Nader et al. 2019")
a2.plot([0.22], [(3 * EH + 20.81) * HA_EV], "^", mfc="none", color="0.4", ms=6, label="Lai et al. 1992 (HF)")
a2.set_xlabel("d, neighbour distance (bohr)"); a2.set_ylabel("D = 3E(H) − E(H$_3$)  (eV)")
a2.set_title("linear H$_3$ quartet, chain ∥ B, B = 425 a.u.", fontsize=9); a2.legend(fontsize=7, frameon=False)
fig.tight_layout()
FIG1 = os.path.join(WORK, "curves.png"); fig.savefig(FIG1, dpi=200); plt.close(fig)

# ------------------------------------------------------------------ document
FT = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf")
pdfmetrics.registerFont(TTFont("DV", os.path.join(FT, "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("DVB", os.path.join(FT, "DejaVuSans-Bold.ttf")))
pdfmetrics.registerFont(TTFont("DVI", os.path.join(FT, "DejaVuSans-Oblique.ttf")))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DVB", italic="DVI", boldItalic="DVB")
ss = getSampleStyleSheet()
body = ParagraphStyle("b", parent=ss["Normal"], fontName="DV", fontSize=9.2, leading=12.5, spaceAfter=4)
small = ParagraphStyle("s", parent=body, fontSize=7.8, leading=10)
cell = ParagraphStyle("c", parent=body, fontSize=8, leading=10, spaceAfter=0)
h1 = ParagraphStyle("h1", parent=body, fontName="DVB", fontSize=12.5, leading=16, spaceBefore=10, spaceAfter=4,
                    textColor=colors.HexColor("#1d3d5c"))
title = ParagraphStyle("t", parent=body, fontName="DVB", fontSize=16, leading=20, spaceAfter=2)
P = lambda t, s=body: Paragraph(t, s)

def table(rows, widths, hl=None):
    data = [[P(re.sub(r"(^|>)-(?=\d)", "\\1\u2212", str(c)), cell) for c in r] for r in rows]
    t = Table(data, colWidths=[w * cm for w in widths], repeatRows=1)
    st = [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eef4")),
          ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#1d3d5c")),
          ("LINEBELOW", (0, -1), (-1, -1), 0.4, colors.grey),
          ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("TOPPADDING", (0, 0), (-1, -1), 2),
          ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    for r in (hl or []):
        st.append(("BACKGROUND", (0, r), (-1, r), colors.HexColor("#fdf3e1")))
    t.setStyle(TableStyle(st))
    return t

f6 = lambda x: "%.6f" % x
s = []
s += [P("Molecules in a strong magnetic field: progress report", title),
      P("MADNESS <i>magneticatom</i> (Fourier route), branch magnetic-fourier · hawaii cluster · "
        "status of %s" % datetime.datetime.now().strftime("%d %B %Y, %H:%M"), small), Spacer(1, 6)]

s += [P("1. Summary", h1),
      P("The code is validated on the hydrogen atom against exact values (errors 8–26 µHa from B = 1 to 1000 a.u.) "
        "and on the H<sub>3</sub> molecule against Bischoff’s Hartree–Fock benchmark at B = 1 a.u. (0.09 mHa). "
        "At B = 425 a.u. (10<super>12</super> G) the two parallel structures are finished: "
        "<b>H<sub>2</sub></b> (bond ∥ B) with R<sub>e</sub> = %.3f bohr and D = %.2f eV, matching Lai’s 46 eV, and "
        "<b>linear H<sub>3</sub></b> (chain ∥ B) with d = %.3f bohr and D = %.2f eV relative to three H atoms, "
        "within 0.07 %% of Nader et al. in total energy. The tilted and triangular geometries and He<sub>2</sub> are "
        "running or queued." % (h2opt[0], (2 * EH - h2opt[1]) * HA_EV, h3opt[0], (3 * EH - h3opt[1]) * HA_EV))]

s += [P("2. Method and conventions", h1),
      P("Unrestricted Hartree–Fock in a uniform field <b>B along z</b>; the centroid of the nuclei is the gauge origin. "
        "All electron spins are anti-parallel to B (high-spin states), and every energy is "
        "E<sub>phys</sub> = E<sub>tot</sub> − N·B/2, i.e. it includes the spin-Zeeman term. Dissociation energies are "
        "D = Σ E(atoms) − E(molecule) (positive = bound), with the atoms computed in run 10 at the same numerical "
        "settings. Nuclei are clamped (no zero-point energy). Atomic units: B<sub>0</sub> = 2.3505×10<super>9</super> G, "
        "so 10<super>12</super> G = 425.4 B<sub>0</sub>; 1 Ha = 27.2114 eV."),
      P("Settings at B = 425: k = 10, thresh 10<super>−6</super>, L = 8, qextra = 10; rcut = 0.75 bohr when all nuclei "
        "are on the field axis (271 Fourier modes), 1.0 bohr otherwise (477 modes). Jobs run on single 16-core nodes "
        "through <i>submit_next.sh</i>, which yields to other users and stops each job after 12 h.")]

s += [P("3. Validation", h1),
      P("<b>3a. Hydrogen atom vs. exact values</b> (Kravchenko, Liberman &amp; Johansson, PRA 54, 287 (1996)). "
        "Binding energy B/2 − E in Ha; our values are slightly underbound at every field.")]
rows = [["B (a.u.)", "run", "this work", "exact", "error (µHa)"]]
for b, bind, job in hlit:
    rows.append([b, job, "%.7f" % bind, "%.7f" % EXACT[b], "%+.1f" % ((bind - EXACT[b]) * 1e6)])
s += [table(rows, [2, 2.6, 3, 3, 2.4]), Spacer(1, 5)]
s += [P("<b>3b. Separated atoms at B = 425 a.u.</b> (zero of all dissociation energies)")]
s += [table([["atom", "rcut 0.75", "rcut 1.0", "difference"],
             ["H (E<sub>phys</sub> = −binding)", f6(EH), f6(EH10), "%.1f µHa" % ((EH10 - EH) * 1e6)],
             ["He (UHF high spin)", f6(EHe), f6(EHe10), "%.1f µHa" % ((EHe10 - EHe) * 1e6)]],
            [4.4, 3, 3, 2.6]),
      P("H binding: %.2f eV (Lai 2001, Table I: Q<sub>1</sub> = 161 eV at 10<super>12</super> G)."
        % (-EH * HA_EV), small)]

s += [P("<b>3c. H<sub>3</sub> at B = 1 a.u. vs. Bischoff, PRA 101, 053413 (2020)</b>"),
      P("Quartet H<sub>3</sub> (UHF), equilateral triangle with side 144.5 pm = %.4f bohr, lying in the xy plane: "
        "the field is perpendicular to the molecular plane, i.e. the normal of the triangle is parallel to B "
        "(“angle (normal, B) = 0°” in Bischoff’s Table II). 631 Fourier modes, 12 SCF cycles, 9.9 h on one node."
        % bis_geo[0])]
D07 = (3 * EH_B1 - bis) * 1e3
s += [table([["", "this work", "Bischoff (Table II)", "difference"],
             ["E(H<sub>3</sub>) (Ha)", "%.7f" % bis, "−2.50237", "%+.1f µHa" % ((bis + 2.50237) * 1e6)],
             ["3 E(H) (Ha)", "%.7f" % (3 * EH_B1), "−2.49349", "%+.1f µHa" % ((3 * EH_B1 + 2.49349) * 1e6)],
             ["D vs. 3 H (mHa)", "%.3f" % D07, "8.88", "%+.0f µHa" % ((D07 - 8.88) * 1e3)]],
            [3.6, 3, 3.4, 2.6]),
      P("Most of the 88 µHa offset in the total energy is our H atom at B = 1 (25.7 µHa per atom, run 21). "
        "Remaining force along the triangle size: dE/dR = %+.1e Ha/bohr (Bischoff optimized to |g| &lt; 10<super>−5</super>); "
        "relaxing it would lower E by only a few µHa. Out-of-plane forces are ~10<super>−15</super>, so the "
        "perpendicular orientation is a stationary point. Bischoff finds the triangle tilting with field "
        "(0° at B = 1, 25° at B = 3, 42° at B = 5); he cites linear H<sub>3</sub> only for B &gt; 50 a.u."
        % bis_geo[2], small)]

s += [KeepTogether([P("4. Results at B = 425 a.u. (10<super>12</super> G)", h1), Image(FIG1, width=17.4 * cm, height=17.4 * cm * 2.9 / 7.2),
                    P("Figure 1. Dissociation energies along the potential curves. Literature points use their own "
                      "E(H<sub>2</sub>), E(H<sub>3</sub>) and our E(H); 10<super>12</super> G is 425.4 a.u., "
                      "0.1 % above our field.", small)])]
s += [P("<b>4a. H<sub>2</sub> spin triplet, bond parallel to B</b> (runs 11, 11opt; state ³Π<sub>u</sub>, "
        "m = 0 and −1 orbitals)")]
rows = [["source", "E(H<sub>2</sub>) (Ha)", "R<sub>e</sub> (bohr)", "D (Ha)", "D (eV)"],
        ["<b>this work, UHF</b>", "<b>%.4f</b>" % h2opt[1], "<b>%.3f</b>" % h2opt[0],
         "<b>%.4f</b>" % (2 * EH - h2opt[1]), "<b>%.2f</b>" % ((2 * EH - h2opt[1]) * HA_EV)],
        ["Lai, RMP 73, 629 (2001), Table I", "—", "—", "—", "46 (±10 %)"],
        ["Lai, Salpeter &amp; Shapiro, PRA 45, 4832 (1992), HF", "−13.55", "0.24", "", ""],
        ["Lai &amp; Salpeter, PRA 53, 152 (1996), HF", "−13.6", "0.25", "", ""],
        ["Nader, López Vieyra &amp; Turbiner (2020), variational", "−13.576", "0.243",
         "%.4f*" % (2 * EH + 13.576), "%.2f*" % ((2 * EH + 13.576) * HA_EV)]]
s += [table(rows, [7.2, 2.6, 2.2, 2.2, 2.2], hl=[1]),
      P("* with our E(H). Lai’s second value, 32 eV, includes the proton zero-point energy and is not comparable "
        "with clamped nuclei. Our energy is 25 mHa below Nader’s 3-term trial function, which those authors place "
        "~0.7 %% above the large-basis results of Detmer et al. (1998). Optimizer stopped at dE/dR = %.1e Ha/bohr; "
        "the 11-point scan (Fig. 1) confirms the minimum near R = 0.24." % h2opt[2], small)]
s += [P("<b>4b. Linear H<sub>3</sub> quartet, chain parallel to B</b> (runs 15, 15opt; state ⁴(−3)<super>+</super>)")]
D_h2h = EH + h2opt[1] - h3opt[1]
rows = [["source", "E(H<sub>3</sub>) (Ha)", "d (bohr)", "D vs. 3 H (eV)", "D(H<sub>3</sub>→H<sub>2</sub>+H) (Ha)"],
        ["<b>this work, UHF</b>", "<b>%.4f</b>" % h3opt[1], "<b>%.3f</b>" % h3opt[0],
         "<b>%.2f</b>" % ((3 * EH - h3opt[1]) * HA_EV), "<b>%.3f</b>" % D_h2h],
        ["Nader, López Vieyra &amp; Turbiner, PRA 100, 012508 (2019)", "−20.90", "0.219",
         "%.2f*" % ((3 * EH + 20.90) * HA_EV), "1.382"],
        ["Lai, Salpeter &amp; Shapiro (1992), HF", "−20.81", "0.22", "%.2f*" % ((3 * EH + 20.81) * HA_EV), ""]]
s += [table(rows, [6.6, 2.4, 1.8, 2.6, 3.0], hl=[1]),
      P("* with our E(H). Our UHF energy lies 15 mHa above the correlated variational value and 75 mHa below the "
        "adiabatic-approximation HF of Lai et al. The H<sub>3</sub>→H<sub>2</sub>+H energy is 4 %% below Nader’s, "
        "because our H<sub>2</sub> is lower and our H<sub>3</sub> higher than their respective trial functions. "
        "Scan points finished: %d (the last point, d = 0.32, was still running when this report was made)." % len(h3),
        small)]

s += [P("5. Status of the remaining calculations", h1)]
rows = [["file", "system", "geometry", "state"],
        ["12opt / 12", "H<sub>2</sub>", "bond at 90° to B", "optimization running · scan queued (16 points)"],
        ["13opt / 13", "H<sub>2</sub>", "bond at 60° to B", "optimization running · scan 1 of 14"],
        ["14opt / 14", "H<sub>2</sub>", "bond at 30° to B", "optimization running · scan started"],
        ["16opt / 16", "H<sub>3</sub>", "triangle, plane ⊥ B (Bischoff orientation)", "optimization step 1 · scan queued"],
        ["17opt / 17", "H<sub>3</sub>", "triangle tilted 30°, 60°, 90° (two ways)", "optimization 1 of 4 · scan queued"],
        ["18", "H<sub>3</sub>", "bent chain, α = 180°…60°", "queued"],
        ["19opt / 19", "He<sub>2</sub>", "bond ∥ B", "queued"],
        ["20opt / 20", "He<sub>2</sub>", "bond at 45°, 90°", "queued"],
        ["01–06", "H", "convergence refinements (qextra, settings)", "partly done, paused for speed"]]
s += [table(rows, [2.3, 1.6, 6.0, 6.6]),
      P("Optimizations relax only the size (R, d or triangle side) at a fixed orientation; the preferred angle "
        "follows from comparing the fixed-angle energies. The off-axis runs need rcut = 1.0 and take ~9–10 h per "
        "calculation (vs. ~1.7 h on axis). Expected: H<sub>2</sub>/H<sub>3</sub> in 4–5 days and He<sub>2</sub> "
        "1–2 weeks more, if no other user needs the nodes.")]

s += [P("6. Open points", h1),
      P("• Run 06 at qextra = 20 gives the exact H reference for the Bischoff comparison (settings of run 07).<br/>"
        "• Test rcut 0.75 vs. 1.0 on one tilted geometry; if they agree, the off-axis work becomes ~5× cheaper.<br/>"
        "• At B = 425 the preferred tilt of the H<sub>3</sub> triangle is unknown; 0/30/60/90° may need refinement.<br/>"
        "• He<sub>2</sub> reference: Lai (HF) D = 51.4 eV; Medin &amp; Lai (DFT) R<sub>eq</sub> = 0.25 bohr.")]

doc = SimpleDocTemplate(OUT, pagesize=A4, leftMargin=1.8 * cm, rightMargin=1.8 * cm, topMargin=1.6 * cm,
                        bottomMargin=1.6 * cm, title="Magnetic-field campaign: progress report")
def footer(c, d):
    c.setFont("DV", 7.5); c.setFillColor(colors.grey)
    c.drawRightString(A4[0] - 1.8 * cm, 1 * cm, "page %d" % d.page)
doc.build(s, onFirstPage=footer, onLaterPages=footer)
print(OUT)
