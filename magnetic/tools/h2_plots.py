"""Quick-look PNGs of the H2 runs at 10^12 G, made from scratch from the logs (not paper figures).

  h2_theta<angle>.png   binding curve per angle: E(H2) - 2E(H) vs R
  h2_vs_angle.png       the same numbers vs angle, one line per R
  h2_iterations.png     SCF iterations and last |dE| per scan point and start
  h2_ropt_theta0.png    0 deg bond-length optimization

Run with the report env:  ~/miniforge3/envs/report/bin/python -I tools/h2_plots.py [outdir]   (default: ~)
"""
import glob, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "paper"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from logs import R, JOBS, HA_EV, read, finished, e_phys, e_atom

OUT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~")
RUNS = [("0", "11_H2_theta0_B425", "0.75"), ("30", "14_H2_theta30_B425", "1.0"), ("40", "30_H2_theta40_B425", "1.0"),
        ("50", "31_H2_theta50_B425", "1.0"), ("60", "13_H2_theta60_B425", "1.0"), ("90", "12_H2_theta90_B425", "1.0"),
        ("10", "32_H2_theta10_B425", "1.0"), ("20", "33_H2_theta20_B425", "1.0"), ("70", "34_H2_theta70_B425", "1.0"),
        ("80", "35_H2_theta80_B425", "1.0")]
NAME = re.compile(r"R(?P<x>[0-9.]+?)(?:_(?P<s>mol|atm))?")
COL = {"mol": "#2a78d6", "atm": "#eb6834"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e6e6", "savefig.dpi": 200})


def points(folder):
    """{R: {start: log}} of every planned scan point: a finished log, or a job file active in order.list
    (not #HOLD/#CUT/#OOM/#REPLACED)."""
    out = {}
    order = open(os.path.join(JOBS, "..", "batches", "order.list")).read()
    names = re.findall(r"^jobs/split/%s__(R\S+)\.jobs" % re.escape(folder), order, re.M)
    names += [os.path.basename(p)[:-4] for p in glob.glob(R(folder, "R*.log")) if finished(p)]
    for n in set(names):
        m = NAME.fullmatch(n)
        if m:
            out.setdefault(float(m.group("x")), {})[m.group("s") or "mol"] = R(folder, n + ".log")
    return out


def scan(folder, rc):
    """[(R, dE eV, all starts finished)] and whether every planned point is finished."""
    e2h, rows, complete = 2 * e_atom("H", rc), [], True
    for x, starts in sorted(points(folder).items()):
        es = [e for e in map(e_phys, starts.values()) if e is not None]
        complete &= len(es) == len(starts)
        if es:
            rows.append((x, (min(es) - e2h) * HA_EV, len(es) == len(starts)))
    return rows, complete


def scf(log):
    """(number of SCF iterations, last |dE| in uHa) of the last SCF of a finished log."""
    t = read(log).split("geo: new R")[-1]
    de = re.findall(r"^ *scf +\d+ .*dE = *([-+0-9.e]+)", t, re.M)
    return len(de), abs(float(de[-1])) * 1e6


def zero_line(ax):
    ax.axhline(0, color="#666666", ls="--", lw=1)
    ax.annotate("2 H atoms", (0.02, 0), xycoords=("axes fraction", "data"), xytext=(0, 4),
                textcoords="offset points", color="#666666", fontsize=9)


data = {th: (folder, scan(folder, rc)) for th, folder, rc in RUNS}
data = {th: v for th, v in data.items() if v[1][0]}
angles = sorted(data, key=float)

# binding curve per angle
for th in angles:
    rows, complete = data[th][1]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot([r[0] for r in rows], [r[1] for r in rows], "-", color=COL["mol"], lw=2, zorder=1)
    for full, kw, lab in ((True, dict(color=COL["mol"]), "all starts finished"),
                          (False, dict(color="white", markeredgecolor=COL["mol"], markeredgewidth=2), "one start only")):
        pts = [r for r in rows if r[2] == full]
        if pts:
            ax.plot([r[0] for r in pts], [r[1] for r in pts], "o", ms=8, label=lab, zorder=2, **kw)
    zero_line(ax)
    rmin = min(rows, key=lambda r: r[1])
    lines = ["minimum %.2f eV at R = %.2f" % (rmin[1], rmin[0]) if rmin[1] < 0 else "unbound at every R"]
    if not complete:
        lines.append("scan not finished: %d points so far" % len(rows))
    ax.set_xlabel("R (bohr)")
    ax.set_ylabel("E(H$_2$) − 2E(H)  (eV)")
    ax.set_title("H$_2$, θ = %s°, B = 425.5 (10$^{12}$ G); negative = bound\n%s" % (th, "; ".join(lines)))
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "h2_theta%s.png" % th))
    plt.close(fig)

# vs angle at fixed R
fig, ax = plt.subplots(figsize=(7, 4.5))
cmap = plt.get_cmap("Blues")
RS = [0.16, 0.2, 0.24, 0.3, 0.36]
for i, x in enumerate(RS):
    pts = [(float(th), d, full) for th in angles for (r, d, full) in data[th][1][0] if abs(r - x) < 1e-9]
    c = cmap(0.4 + 0.6 * i / (len(RS) - 1))
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", color=c, lw=2, label="R = %.2f" % x)
    for p in pts:
        ax.plot(p[0], p[1], "o", ms=7, color=c if p[2] else "white", markeredgecolor=c, markeredgewidth=2)
zero_line(ax)
ax.set_xticks([float(t) for t in angles])
ax.set_xlabel("θ (deg), angle between bond and B")
ax.set_ylabel("E(H$_2$) − 2E(H)  (eV)")
ax.set_title("H$_2$ vs orientation, 10$^{12}$ G (R in bohr; open = one start only)")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "h2_vs_angle.png"))
plt.close(fig)

# SCF iterations and last |dE| per scan point
fig, axs = plt.subplots(2, len(angles), figsize=(3.4 * len(angles), 6.5), sharex=True, sharey="row", squeeze=False)
for j, th in enumerate(angles):
    pts = points(data[th][0])
    for s in ("mol", "atm"):
        rows = sorted((x, *scf(st[s])) for x, st in pts.items() if s in st and finished(st[s]))
        if rows:
            lab = "molecular start" if s == "mol" else "atomic start"
            axs[0, j].plot([r[0] for r in rows], [r[1] for r in rows], "o-", color=COL[s], label=lab)
            axs[1, j].plot([r[0] for r in rows], [r[2] for r in rows], "o-", color=COL[s])
    axs[0, j].set_title("θ = %s°" % th)
    axs[1, j].axhline(10, color="#666666", ls="--", lw=1)
    axs[1, j].set_yscale("log")
    axs[1, j].set_xlabel("R (bohr)")
axs[1, 0].annotate("tolscf 10 µHa", (0.02, 10), xycoords=("axes fraction", "data"), xytext=(0, 3),
                   textcoords="offset points", fontsize=8, color="#666666")
axs[0, 0].set_ylabel("SCF iterations")
axs[1, 0].set_ylabel("last |dE| (µHa)")
axs[0, -1].legend(frameon=False, fontsize=9)
fig.suptitle("H$_2$, B = 425.5: SCF iterations and final energy change per scan point (finished runs)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "h2_iterations.png"))
plt.close(fig)

# 0 deg bond-length optimization
log = R("11opt_H2_theta0_B425", "ropt_R0.24.log")
steps = read(log).split("geo: new R")
geo = re.findall(r"^geo +(\d+) +R = *([0-9.]+) .*dE/dR \(HF\) = *([-+0-9.e]+)", read(log), re.M)
fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
for k, chunk in enumerate(steps[:len(geo)]):
    de = [abs(float(d)) * 1e6 for d in re.findall(r"^ *scf +\d+ .*dE = *([-+0-9.e]+)", chunk, re.M)[1:]]
    axs[0].plot(range(1, len(de) + 1), de, "o-", label="step %s: R = %.4f" % (geo[k][0], float(geo[k][1])))
axs[0].axhline(10, color="#666666", ls="--", lw=1)
axs[0].set(yscale="log", xlabel="SCF iteration", ylabel="|dE| (µHa)", title="SCF convergence per geometry step")
axs[0].xaxis.get_major_locator().set_params(integer=True)
axs[0].legend(frameon=False, fontsize=8)
axs[1].plot([int(g[0]) for g in geo], [float(g[1]) for g in geo], "o-", color=COL["mol"])
axs[1].set(xlabel="geometry step", ylabel="R (bohr)", title="bond length")
axs[2].plot([int(g[0]) for g in geo], [abs(float(g[2])) for g in geo], "o-", color=COL["mol"])
axs[2].set(yscale="log", xlabel="geometry step", ylabel="|dE/dR| (Ha/bohr)", title="gradient")
for ax in axs[1:]:
    ax.set_xticks([int(g[0]) for g in geo])
fig.suptitle("H$_2$, θ = 0°, B = 425.5: bond-length optimization%s" % ("" if finished(log) else " (running)"))
fig.tight_layout()
fig.savefig(os.path.join(OUT, "h2_ropt_theta0.png"))
plt.close(fig)

print("written to", OUT, ":", ", ".join("h2_theta%s" % t for t in angles), "+ h2_vs_angle, h2_iterations, h2_ropt_theta0")
