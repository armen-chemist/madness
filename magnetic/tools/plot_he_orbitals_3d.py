#!/usr/bin/env python
"""He atom orbitals in strong B (run 42) as 3D isosurfaces, laid out like Fig. 1 of
Farhaz et al., JCP 163, 034308 (2025): column groups = states (3P0, 3P1), columns = B,
rows = Re orbital-2, Re orbital-1 (real part) and |orbital-2|, |orbital-1| (absolute values).

Real part: two isosurfaces at +/- ISO * max|Re| (blue = +, orange = -). Absolute value: one
isosurface at ISO * max|.|. Each panel is scaled to its own orbital (orbitals shrink with B),
B points up. Reads the full 3D grids in results/42_He_orbitals/<state>_B<B>.orb.

Usage (conda env report, needs scikit-image):  python tools/plot_he_orbitals_3d.py [out_dir]
    -> he_orbitals_3d.png
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_he_orbitals import FIELDS  # noqa: E402

_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.environ.get("HE_ORB_DIR", os.path.join(_HERE, "results", "43_He_orbitals_hires"))

ISO = float(os.environ.get("HE_ISO", "0.25"))
STATES = [("S1", r"$^1S_0$"), ("S0T", r"$^3S_0$"), ("P0", r"$^3P_0$"), ("P1", r"$^3P_1$")]   # S1 = singlet 1s^2: one orbital
POS, NEG = to_rgb("#f2b36b"), to_rgb("#8fc3ea")
LIGHT = LightSource(azdeg=315, altdeg=40)
ALPHA = float(os.environ.get("HE_ALPHA", "0.45"))   # see-through surfaces, nucleus visible inside


def read_grid(path):
    with open(path) as f:
        head = f.readline().split()
    kv = dict(x.split("=") for x in head if "=" in x)
    n, nxy, nz = int(kv["N"]), int(kv["nxy"]), int(kv["nz"])
    cache = path + ".npy"                              # parsing the large text files once is the slow part
    if os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(path):
        d = np.load(cache)
    else:
        with open(path) as f:
            txt = "".join(l for l in f if not l.startswith("#"))
        d = np.fromstring(txt, sep=" ").reshape(-1, 5)
        np.save(cache, d)
    m = nxy * nxy * nz
    x = d[:m, 0].reshape(nxy, nxy, nz)[:, 0, 0]
    y = d[:m, 1].reshape(nxy, nxy, nz)[0, :, 0]
    z = d[:m, 2].reshape(nxy, nxy, nz)[0, 0, :]
    orbs = []
    for a in range(n):
        blk = d[a * m:(a + 1) * m]
        v = (blk[:, 3] + 1j * blk[:, 4]).reshape(nxy, nxy, nz)
        # global phase: largest value on the half-plane y = 0, x >= 0 made real and positive. For 2p-1 ~ (x - i y)
        # this puts the Re lobes along x (+x blue, -x orange) at every B; before, the phase came from the largest
        # value anywhere on the ring, so the Re lobe pair turned about B arbitrarily from panel to panel.
        h = v[nxy // 2:, nxy // 2, :]
        i, l = np.unravel_index(np.argmax(np.abs(h)), h.shape)
        k = (nxy // 2 + i, nxy // 2, l)
        orbs.append(v * np.conj(v[k]) / np.abs(v[k]))
    return x, y, z, orbs


def surface(ax, f, level, x, y, z, color):
    if not (f.min() < level < f.max()):
        return None
    sp = (x[1] - x[0], y[1] - y[0], z[1] - z[0])
    verts, faces, normals, _ = marching_cubes(f, level, spacing=sp)
    verts += (x[0], y[0], z[0])
    tri = verts[faces]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-30
    if level < 0:
        n = -n
    shade = LIGHT.shade_normals(n, fraction=1.0)
    cols = np.clip(np.array(color)[None, :] * (0.45 + 0.65 * shade[:, None]), 0, 1)
    ax.add_collection3d(Poly3DCollection(tri, facecolors=np.c_[cols, np.full(len(cols), ALPHA)], edgecolor="none"))
    return verts


def panel(ax, field, kind, x, y, z):
    pts = []
    if kind == "Re":
        f = field.real
        m = np.abs(f).max()
        for lev, col in ((ISO * m, NEG), (-ISO * m, POS)):   # + blue, - orange (1s blue as in Farhaz Fig. 1)
            v = surface(ax, f, lev, x, y, z, col)
            if v is not None:
                pts.append(v)
    else:
        f = np.abs(field)
        v = surface(ax, f, ISO * f.max(), x, y, z, NEG)
        if v is not None:
            pts.append(v)
    p = np.vstack(pts)
    r = np.abs(p).max(axis=0)
    h = max(r) * 1.05                                 # one cube per panel: true proportions
    ax.set_xlim(-h, h); ax.set_ylim(-h, h); ax.set_zlim(-h, h)
    ax.set_box_aspect((1, 1, 1), zoom=1.45)   # 3D axes leave a wide margin; zoom in so panels sit close
    u, w = np.mgrid[0:2 * np.pi:24j, 0:np.pi:12j]               # nucleus marker at the origin
    rn = 0.06 * h
    ax.plot_surface(rn * np.cos(u) * np.sin(w), rn * np.sin(u) * np.sin(w), rn * np.cos(w),
                    color="#3d5f86", shade=True, linewidth=0, zorder=10)
    ax.view_init(elev=22, azim=-90)   # look along +y: x to the right, B up
    ax.set_axis_off()
    # truncated by the plotted box?
    edge = (abs(abs(p[:, 2]).max() - abs(z).max()) < 1.5 * (z[1] - z[0])) or \
           (abs(abs(p[:, 0]).max() - abs(x).max()) < 1.5 * (x[1] - x[0]))
    return edge


def main(out):
    cols = [(s, lab, b) for s, lab in STATES for b in FIELDS
            if os.path.exists(os.path.join(RES, f"{s}_B{b}.orb"))]
    rows = [("Re", 1), ("Re", 0), ("abs", 1), ("abs", 0)]
    fig = plt.figure(figsize=(0.8 * len(cols) + 1.0, 1.0 * 4 + 1.1))
    gs = fig.add_gridspec(4, len(cols), left=0.07, right=0.95, top=0.88, bottom=0.08,
                          wspace=-0.05, hspace=0.0)
    cut = []
    for c, (s, lab, b) in enumerate(cols):
        x, y, z, orbs = read_grid(os.path.join(RES, f"{s}_B{b}.orb"))
        for r, (kind, a) in enumerate(rows):
            ax = fig.add_subplot(gs[r, c], projection="3d")
            ax.set_facecolor("none")
            if a >= len(orbs):                         # singlet: only one (doubly occupied) orbital
                ax.set_axis_off()
                continue
            if panel(ax, orbs[a], kind, x, y, z):
                cut.append(f"{s} B={b} {kind} orb-{a + 1}")
        fig.text((gs[3, c].get_position(fig).x0 + gs[3, c].get_position(fig).x1) / 2, 0.055, b,
                 ha="center", fontsize=9)
    for s, lab in STATES:
        idx = [i for i, cc in enumerate(cols) if cc[0] == s]
        if idx:
            x0 = gs[0, idx[0]].get_position(fig).x0
            x1 = gs[0, idx[-1]].get_position(fig).x1
            fig.text((x0 + x1) / 2, 0.92, lab, ha="center", fontsize=15)
    for r, lab in enumerate(["orbital-2", "orbital-1", "orbital-2", "orbital-1"]):
        pos = gs[r, 0].get_position(fig)
        fig.text(0.03, (pos.y0 + pos.y1) / 2, lab, rotation=90, ha="center", va="center",
                 fontsize=10, fontweight="bold", family="serif")
    for r0, lab in ((0, "Real part"), (2, "Absolute values")):
        y0 = gs[r0 + 1, -1].get_position(fig).y0
        y1 = gs[r0, -1].get_position(fig).y1
        fig.text(0.975, (y0 + y1) / 2, lab, rotation=-90, ha="center", va="center",
                 fontsize=11, fontweight="bold", family="serif")
    ymid = gs[2, 0].get_position(fig).y1
    fig.add_artist(plt.Line2D([0.05, 0.96], [ymid, ymid], color="k", lw=0.6))
    for i in range(1, len(cols)):                      # vertical line between state groups
        if cols[i][0] != cols[i - 1][0]:
            xs = (gs[0, i - 1].get_position(fig).x1 + gs[0, i].get_position(fig).x0) / 2
            fig.add_artist(plt.Line2D([xs, xs], [0.07, 0.95], color="k", lw=0.6))
    fig.text(0.5, 0.015, r"$B/B_0$", ha="center", fontsize=11, fontweight="bold")
    for ext in ("png",):                              # user: png only
        fig.savefig(os.path.join(out, f"he_orbitals_3d.{ext}"), dpi=300)
    print(f"{len(cols)} columns -> {out}/he_orbitals_3d.png  (iso = {ISO} x max)")
    if cut:
        print("isosurface reaches the plotted box (cut):", "; ".join(cut))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "docs")
    main(out)
