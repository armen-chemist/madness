#!/usr/bin/env python
"""He atom orbitals in strong B (run 42), laid out like Fig. 1 of Farhaz et al., JCP 163, 034308 (2025).

Reads results/42_He_orbitals/<state>_B<B>.orb (written by magneticatom_plot, plot=...), cuts every orbital in the
x-z plane (y = 0, B along z) and draws, per state, columns = B and rows = Re orbital-2, Re orbital-1,
|orbital-2|, |orbital-1|. The global phase of each orbital is fixed so that its largest value is real and positive.
Each panel has its own length scale (orbitals shrink with B); the axis ticks are in bohr.

Usage (conda env report):  python tools/plot_he_orbitals.py [out_dir]     -> he_orbitals_P1.pdf/png, he_orbitals_P0.pdf/png
"""
import os
import sys

import numpy as np
import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.environ.get("HE_ORB_DIR", os.path.join(HERE, "results", "42_He_orbitals"))
FIELDS = ["10", "100", "425.5", "1000", "4255"]
STATES = {"P1": r"$^3P_1$: 1s 2p$_{-1}$", "P0": r"$^3P_0$: 1s 2p$_0$"}


def read_orb(path):
    with open(path) as f:
        head = f.readline().split()
        eps = [float(x) for x in f.readline().split()[2:]]
    kv = dict(x.split("=") for x in head if "=" in x)
    n, nxy, nz = int(kv["N"]), int(kv["nxy"]), int(kv["nz"])
    d = np.loadtxt(path, comments="#")
    if d.shape[0] != n * nxy * nxy * nz:
        raise ValueError(f"{path}: {d.shape[0]} points, expected {n * nxy * nxy * nz}")
    x = d[: nxy * nxy * nz, 0].reshape(nxy, nxy, nz)[:, 0, 0]
    z = d[: nxy * nxy * nz, 2].reshape(nxy, nxy, nz)[0, 0, :]
    orbs = []
    for a in range(n):
        blk = d[a * nxy * nxy * nz:(a + 1) * nxy * nxy * nz]
        v = (blk[:, 3] + 1j * blk[:, 4]).reshape(nxy, nxy, nz)
        # global phase: largest value on the half-plane y = 0, x >= 0 made real and positive. For 2p-1 ~ (x - i y)
        # this puts the Re lobes along x (+x blue, -x orange) at every B; before, the phase came from the largest
        # value anywhere on the ring, so the Re lobe pair turned about B arbitrarily from panel to panel.
        h = v[nxy // 2:, nxy // 2, :]
        i, l = np.unravel_index(np.argmax(np.abs(h)), h.shape)
        k = (nxy // 2 + i, nxy // 2, l)
        v = v * np.conj(v[k]) / np.abs(v[k])       # global phase: largest value real and positive
        orbs.append(v[:, nxy // 2, :])               # x-z cut at y = 0
    return float(kv["B"]), x, z, orbs, eps


def nice(v):
    """largest 1/2/5 x 10^n not above v"""
    e = 10 ** np.floor(np.log10(v))
    for m in (5, 2, 1):
        if m * e <= v:
            return m * e
    return e


def figure(state, out):
    cols = [(b, os.path.join(RES, f"{state}_B{b}.orb")) for b in FIELDS]
    cols = [(b, p) for b, p in cols if os.path.exists(p)]
    if not cols:
        print(f"{state}: no .orb files yet")
        return
    data = [read_orb(p) for _, p in cols]
    rows = [("Re", 1), ("Re", 0), ("abs", 1), ("abs", 0)]
    rowlab = [r"Re $\varphi_2$ (2p)", r"Re $\varphi_1$ (1s)", r"$|\varphi_2|$ (2p)", r"$|\varphi_1|$ (1s)"]
    ratios = [(x[-1] - x[0]) / (z[-1] - z[0]) for _, x, z, _, _ in data]
    ph = 2.3                                            # panel height, inch
    fig = plt.figure(figsize=(sum(ratios) * ph + 0.75 * len(cols) + 0.8, 4 * ph + 1.2), layout="constrained")
    gs = fig.add_gridspec(4, len(cols), width_ratios=ratios)
    for c, (b, _) in enumerate(cols):
        B, x, z, orbs, eps = data[c]
        tx, tz = nice(0.85 * x[-1]), nice(0.6 * z[-1])
        for r, (kind, a) in enumerate(rows):
            ax = fig.add_subplot(gs[r, c])
            v = orbs[a]
            if kind == "Re":
                f = v.real
                m = np.abs(f).max()
                ax.pcolormesh(x, z, f.T, shading="gouraud", cmap="RdBu_r", vmin=-m, vmax=m, rasterized=True)
            else:
                f = np.abs(v)
                ax.pcolormesh(x, z, f.T, shading="gouraud", cmap="Blues", vmin=0, vmax=f.max(), rasterized=True)
            ax.contour(x, z, f.T, levels=np.linspace(-1, 1, 11)[1:-1] * np.abs(f).max(),
                       colors="k", linewidths=0.3, alpha=0.45)
            ax.set_aspect("equal")
            ax.set_xticks([-tx, 0, tx])
            ax.set_xticklabels(["", "", f"{tx:g}"])
            ax.set_yticks([-tz, 0, tz])
            ax.tick_params(labelsize=8, length=2.5)
            if r == 3:
                ax.set_xlabel(r"$x$ (bohr)", fontsize=9)
            if c == 0:
                ax.set_ylabel(rowlab[r] + "\n" + r"$z$ (bohr)", fontsize=9)
            if r == 0:
                ax.set_title(f"$B$ = {b} $B_0$", fontsize=10)
    fig.suptitle(f"He {STATES[state]}, UHF, MRA  (x–z plane, y = 0; $B \\parallel z$)", fontsize=11)
    fig.savefig(os.path.join(out, f"he_orbitals_{state}.png"), dpi=300)
    print(f"{state}: {len(cols)} fields -> {out}/he_orbitals_{state}.png")


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "docs")
    os.makedirs(out, exist_ok=True)
    for s in STATES:
        figure(s, out)
