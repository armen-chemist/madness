# Magnetic-field MADNESS: Fourier approach

Branch `magnetic-fourier` of MADNESS (upstream commit 37fa368cd) with a solver for atoms and small
molecules in a strong uniform magnetic field B (along z). The Green's function of the kinetic
operator (p + A)^2/2 is applied exactly, by a Fourier separation in y, a mask-convolution-mask in
(x, z) and a proper-time quadrature; the field never appears in the potential.

| Where | What |
|---|---|
| `src/madness/mra/magnetic_green.h` | magnetic Green's function (Fourier approach) |
| `src/madness/chem/magneticatom.cc` | the program `magneticatom` (H, He, H2, H3, He2, ...; RHF/UHF, forces) |
| `src/madness/mra/funcimpl.h` | changed: `hartree_product` / `project_out` with mixed dimensions (2D x 1D <-> 3D) |
| `magnetic/` | this folder: build script, job files, run and summary scripts |
| `magnetic/docs/parameters_physics.pdf` | meaning of every parameter |

The changes against upstream MADNESS (once: `git remote add upstream https://github.com/m-a-d-n-e-s-s/madness.git`):

```
git fetch upstream
git diff upstream/master...magnetic-fourier -- src
```

To bring in the latest MADNESS: `git fetch upstream && git merge upstream/master` on `magnetic-fourier`,
then rebuild (section 2) and push (`git push origin magnetic-fourier`).


## 1. Get the code

```
git clone -b magnetic-fourier https://github.com/armen-chemist/madness.git ~/madness
cd ~/madness/magnetic
```


## 2. Build the program (once, 30-60 minutes)

### 2a. Normal case (the machine can reach the internet)

```
cd ~/madness/magnetic
bash install.sh
```

What it does, in 3 steps (each prints `=== N.`):

1. installs Miniforge (a free conda) in `~/miniforge3` (only if there is no conda yet) and creates the environment `madness`
   with gcc 14, OpenMPI, OpenBLAS and cmake 3 (about 1 GB download);
2. compiles in `~/madness/build` (logs: `cmake.log` and `build.log` there);
3. runs a short test (4-10 minutes) and prints **`PASSED`** or **`FAILED`**.

You must see `PASSED` at the end. The program is then

```
~/madness/build/src/madness/chem/magneticatom
```

Options: `bash install.sh /other/build 24` builds in `/other/build` with 24 cores.
Then put `export MAGNETICATOM=/other/build/src/madness/chem/magneticatom` into your `~/.bashrc`.

### 2b. Without internet (conda cannot download)

Use the compilers of the machine: a C++ compiler with MPI (`mpicxx`), cmake >= 3.24, and LAPACK/BLAS. For example:

```
module load gcc openmpi openblas cmake                 # names differ between servers
USE_CONDA=0 LAPACK_LIBS="/path/to/liblapack.so;/path/to/libopenblas.so" bash install.sh
```

### 2c. If something fails

- `ERROR: cmake failed` -> look at the last 30 lines of `~/madness/build/cmake.log`.
- `ERROR: build failed` -> look for the first `error:` in `~/madness/build/build.log`.
- To start fresh: `rm -rf ~/madness/build` and run `bash install.sh` again.


## 3. Run calculations

### 3a. Start a job file

```
cd ~/madness/magnetic
nohup bash run_queue.sh jobs/01_H_B100_qmax.jobs > queue_01.out 2>&1 &
```

- `nohup ... &` = runs in the background and continues after you log out (you can close the terminal).
- The jobs of the file run **one after another**.
- Each job writes `results/01_H_B100_qmax/<job name>.log`.
- When the whole file is done, the file `results/01_H_B100_qmax/DONE` appears.
- You do **not** need to activate conda before running (tested: the program finds its libraries by itself).

### 3b. Check progress

```
bash status.sh                              # what runs, memory, newest lines of every queue
cat queue_01.out                            # start/end time and result of each job
grep iter results/01_H_B100_qmax/qextra120.log     # all iterations of one job
```

### 3c. Results table

```
bash summarize.sh results/01_H_B100_qmax 3.789804236305
```

The second number is the exact binding energy (optional); the column `err(uHa)` is then
(our binding − exact) in µHa.

### 3d. Stop and restart

```
bash stop.sh                                # stops ALL queues and calculations
```

Starting the same command of 3a again continues the file: jobs that ended normally
(their log ends with `# JOB FINISHED`) are skipped, all others start again from the beginning.

### 3e. Memory protection

While a job runs, the queue checks the free memory every 10 s. If less than 3 GB is free, it kills the
job, writes `KILLED: low memory` into its log, and continues with the next job.
Change the limit with `MINMEM_GB=8 nohup bash run_queue.sh ...`.

### 3f. Two queues at the same time

Only if both fit into memory (see section 6). Give each half of the cores:

```
MAD_NUM_THREADS=24 nohup bash run_queue.sh jobs/01_H_B100_qmax.jobs > queue_01.out 2>&1 &
MAD_NUM_THREADS=24 nohup bash run_queue.sh jobs/05_H_B500_qmax.jobs  > queue_05.out 2>&1 &
```


## 4. The job files (`jobs/`)

All job files are written by `tools/make_jobs.py` (geometries computed and checked there; edit it and run
`python3 tools/make_jobs.py` to change a scan). The comment lines at the top of each file say what it is for
and what to compare with.

### 4a. Hydrogen atom: Fourier cutoff and all other settings (exact references)

| File | What | Compare with |
|---|---|---|
| `01_H_B100_qmax` | B = 100 B0, qextra = 10 … 320 (qmax 97 … 407) | exact binding 3.789804236305 Ha |
| `02_H_B100_settings` | B = 100 B0, one setting changed per job: k 8/12, thresh 1e-6/1e-8, L 6/12 (at qextra 80); rcut 2.13, rtaper 0.4, dphi 1.5 (at qextra 10) | file 01 (same qextra) |
| `03_H_B425_rcut_L_thresh7` | B = 425.5 B0: rcut 1.0, 1.2 × L 4, 8, 12 at thresh 1e-7 (Florian's question) | each other (independent of rcut and L?) |
| `04_H_B425_qmax` | B = 425.5 B0, qextra 10 … 200 | files 03, 10 |
| `05_H_B500_qmax` | B = 500 B0, qextra 10 … 200 | exact 6.257087674681 Ha |
| `06_H_B1_qmax` | B = 1 B0, qextra 10 … 60 (weak field: qmax set by qextra) | exact 0.831168896733 Ha |
| `07_H3_B1_Bischoff` | cyclic H3 quartet, B = 1 B0, Bischoff's geometry (validation of 3 electrons) | Bischoff, PRA 101, 053413 (2020): −2.50237 Eh (same theory: UHF) |

Exact H values: Kravchenko, Liberman, Johansson, PRA 54, 287 (1996), Table I.

### 4b. Molecules at B = 425.5 B0 (= 1.0×10¹² G)

**One set of numerical settings for all of them** (k 10, thresh 1e-6, L 8, eprec 1e-9, qextra 10, inner 3,
tolscf 1e-5). Only `rcut` differs: 0.75 when all nuclei are on the B axis, 1.0 when nuclei are off the axis
(rule: rcut ≥ 0.70 + largest distance of a nucleus from the axis; checked by `make_jobs.py`).
File 10 computes the separated atoms at **both** rcut values with the same settings: they must agree to < 1 µHa,
and they are the zero of every dissociation energy D = E(atoms) − E(molecule) (D > 0: bound).

Molecule files 11–17, 19, 20 end with `ropt_*` jobs: structure optimization with Hellmann–Feynman forces (section 5a),
one per orientation and start. The R scans in the same file are the independent check (the minimum of the
scan and the optimized R must agree within the scan spacing). Optimization jobs take 4–8 SCFs, i.e. 3–6× a scan point.

Orientation convention: B along z; the centroid of the nuclei at the origin (= gauge origin); bonds and
planes tilted in the x–z plane. Spins: all anti-parallel to B (spin-polarized UHF). This is the ground state at this field
for H2 (Lai 2001), He (Al-Hujaj & Schmelcher, PRA 67, 023403 (2003)), H3 (Nader 2019) and He2 (Medin & Lai,
PRA 74, 062507 (2006): all electrons in tightly bound orbitals).

Two kinds of start orbitals (they only start the SCF, nothing is restricted afterwards):
`mol` = orbitals around the B axis with Lz = 0, −1, −2, … ; `atm` = every atom in its own ground configuration,
moved to its nucleus with the magnetic phase. Where both are run, **the lower converged E_phys is the result**
for that geometry; if both converge to the same energy the state is unambiguous.

| File | System | Geometries | Start | Reference |
|---|---|---|---|---|
| `10_atoms_B425` | H atom, He atom (1s0 2p−1) | rcut 0.75 and 1.0 | | each other (rcut 0.75 vs 1.0) |
| `11_H2_theta0_B425` | H2 triplet, bond ∥ B | R = 0.16 … 0.50 (10) + bond-length search | mol | Lai RMP 73, 629 (2001): Q2 = 46 eV (10 %) |
| `12_H2_theta90_B425` | H2, bond ⟂ B | R = 0.20 … 0.60 (7) | mol + atm | none (is it bound?) |
| `13_H2_theta60_B425` | H2, 60° | same R | mol + atm | |
| `14_H2_theta30_B425` | H2, 30° | same R | mol + atm | |
| `15_H3_linear_B425` | H3 quartet, linear chain ∥ B | spacing d = 0.17 … 0.32 (8) | mol | Nader et al., PRA 100, 012508 (2019): −20.90 Ha, R_eq 0.219 (correlated variational; both upper bounds, expect ~0.1 Ha agreement) |
| `16_H3_triangle_perp_B425` | cyclic H3, plane ⟂ B | side a = 0.22 … 0.45 (5) | mol + atm | |
| `17_H3_triangle_tilt_B425` | cyclic H3, plane tilted 30°, 60°, 90° (edge ∥ B), 90° (vertex on axis) | a = 0.22 … 0.45 | atm | |
| `18_H3_bent_B425` | H3 bent, bond 0.22, angle H–H–H 180° … 60° | 5 angles | mol + atm | 180° = file 15 at d 0.22; 60° = file 17 x90 at a 0.22 (built-in checks) |
| `19_He2_theta0_B425` | He2, 4 spins, bond ∥ B | R = 0.22 … 0.90 (9) + 3 atm checks + bond-length search | mol (Lz 0…−3) | HF (Lai 2001): D = 51.4 eV; DFT (Medin & Lai 2006): R_eq 0.25 |
| `20_He2_angles_B425` | He2, 45° and 90° | R = 0.30, 0.40, 0.55 | mol + atm | none at this field; Austad et al. PCCP 2020 (FCI, B = 1 B0): quintet bound only ⟂ B |

### 4c. Recommended order and cost (rough estimates for 48 cores, from test runs on a 19-thread machine)

| Order | File | Jobs | Time per job | Memory per job |
|---|---|---|---|---|
| 1 | `10` atoms | 4 | H 0.5 h, He 1–2 h | 10–40 GB |
| 2 | `11` H2 ∥ B | 11 | 1–2 h (search 5–8 h) | ~30 GB |
| 3 | `01`, `05` H exact | 13 | 0.5–2 h | 10–60 GB (rises with qextra) |
| 4 | `15` H3 linear | 9 | 2–4 h | ~30 GB |
| 5 | `12`, `13`, `14` H2 tilted | 48 | 1–3 h | 40–60 GB |
| 6 | `16`, `17`, `18` H3 cyclic/bent | 46 | 3–6 h | 50–80 GB |
| 7 | `19` He2 ∥ B | 13 | 6–12 h | ~40 GB |
| 8 | `20` He2 tilted | 16 | 12–24 h | 80–120 GB |
| any time | `02`, `03`, `04`, `06`, `07` | | 0.5–20 h | up to ~100 GB (`07`) |

He2 test: one SCF cycle with 1 inner step took 19 min on 19 threads even at low precision (k 6, thresh 1e-4).
These are estimates; the first job of each file tells the real numbers (`wall` in the summary, memory in `status.sh`).
Measured: H2 with rcut 1.0 (477 modes) already needs 20 GB at LOW precision (k 6, thresh 1e-4). Files with
rcut 1.0 (12–14, 16–18, 20): run ONE queue at a time. Two queues only for rcut 0.75 files and H-atom files.

### 4d. Results tables

```
bash summarize_all.sh > all_results.txt     # every folder; D in Ha and eV filled in automatically from file 10
bash summarize.sh results/11_H2_theta0_B425 <REF>     # one folder; REF = 2 x E_phys(H) from file 10, by hand
```

`status` column: `ok` = converged; `NOT-CONVERGED` = stopped at maxiter/maxscf (do not use the number, run
it again with a larger `maxscf`); `running/stopped`; `KILLED: low memory`; `ERROR`.
A value in parentheses is the last SCF energy of an unfinished job.

### Make your own job file

Copy one and edit it (or edit `tools/make_jobs.py`). Format:

```
# comment
COMMON  <arguments shared by all jobs of this file>
<name>  <arguments of this job>          # a later argument overrides the same one in COMMON
```

Arguments containing `;` (the `atoms=` list) must be in double quotes: `"atoms=1,0,0,-0.12;1,0,0,0.12"`.


## 5. The program `magneticatom`

All arguments are `name=value`. Atomic units: lengths in bohr, energies in Hartree,
B in units of B0 = 2.35 × 10⁵ T = 2.35 × 10⁹ G (the atomic unit). Everywhere in this package B is given in B0;
B = 425.5 B0 = 1.0 × 10⁸ T = 1.0 × 10¹² G.

| Argument | Meaning | Values we use |
|---|---|---|
| `B` | magnetic field, along z | 1, 100, 425.5, 500 |
| `atoms` | nuclei: `Z,x,y,z;Z,x,y,z;...` | default `1,0,0,0` (one H at the origin) |
| `k` | polynomial order of MADNESS | 10 |
| `thresh` | MADNESS accuracy | 1e-7 (H convergence files), 1e-6 (molecules and their atoms) |
| `L` | box from −L to +L in every direction (bohr) | 4–8 |
| `eprec` | smoothing of the nuclear −Z/r | 1e-9 |
| `rcut`, `rtaper` | radius around the B axis where the magnetic phase is kept, width of its edge | see job files |
| `qextra` | extra Fourier modes: qmax = B·(rcut + 3·rtaper)/2 + qextra (if the box is larger than rcut + 3·rtaper) | 10 … 320 |
| `h`, `tmin`, `tmax_factor` | time grid of the Green's function (no effect found) | 0.25, 1e-8, 45 |
| `E0` | starting energy; **must be below B/2** | B/2 − (rough binding) |
| `tol`, `maxiter` | stop when the energy changes by less than tol | 1e-7, 20 |
| `ref` | reference energy: 1 electron = exact total E (log prints `E - ref`, `<H> - ref`); UHF = reference E_phys (log prints `E - ref`) | |
| `nelec`, `state` | number of electrons; `state=triplet` (2) or `state=highspin` (3 or more): all spins parallel | |
| `axial` | starting orbitals only: 1 = around the B axis (Lz = `ma`, `mb` or `mlist`), 0 = per atom (`mloc`) or σg/σu | |
| `ma`, `mb`, `E0b` | 2 electrons, `axial=1`: Lz and starting energy of the 2nd orbital | 0, −1 |
| `mlist`, `E0list` | any number of electrons, `axial=1`: Lz of each start orbital, and their starting energies | He2: `0,-1,-2,-3` |
| `mloc` | `axial=0`: Lz of the start orbitals on EVERY nucleus (atomic start); electrons = nuclei × entries | H: `0`; He: `0,-1` |
| `inner`, `maxscf`, `tolscf` | molecule SCF: inner steps, max. cycles, convergence of E_phys | 3, 30–40, 1e-5 |
| `ropt`, `maxgeo`, `gtol`, `dr0`, `maxstep`, `rtol` | structure optimization with Hellmann–Feynman forces along the breathing coordinate (see 5a) | 1, 10, 2e-3, 0.02, 0.05, 1e-4 |
| `dphi` | grid refinement for the magnetic phase | 3 (default) |

The physical meaning of every parameter: `docs/parameters_physics.pdf`.

### 5a. Structure optimization (`ropt=1`, jobs named `ropt_*`)

Forces: Hellmann–Feynman, g_A = ⟨ρ| ∂V_nuc/∂R_A⟩ + ∂E_nn/∂R_A. This is exact for the converged UHF here: with a
fixed gauge origin only V_nuc and E_nn depend on the nuclei, and MADNESS has no atom-centred basis (no Pulay
terms). Earlier test (H2, B = 425.5 B0): force 0.164 vs finite difference of energies 0.165.

Coordinate: all nuclei are scaled about their centroid (shape and orientation kept); R = distance between the
first two nuclei (bond length for diatomics, spacing d for linear H3, side a for triangles).
dE/dR = Σ_A g_A·(R_A − c)/R; steps by the secant method, the SCF restarts from the current orbitals.
Stops when |dE/dR| < gtol, or when a nucleus would come closer than 0.70 to the rcut edge (`geo: STOP`):
then there is no minimum inside the allowed range (e.g. an unbound orientation) — look at the R scan.

`|F_other|` in the `geo` line = the part of the forces outside the breathing coordinate. At the optimum, small
`|F_other|` means the symmetric structure is stationary; a large one means it wants to change shape or turn
(e.g. a tilted bond rotating towards B). This is reported, not optimized.

### Reading a log

One electron:

```
  iter  4  E = ...  dE = ...  |psi - psi~| = ...       E = total energy (includes B/2)
final E = ...   binding B/2 - E = ...   iterations = ...   wall = ... s   converged = yes
<H> abgv = ...  <H> - E = ...                         independent check; should be small (~1e-7)
```

The result is the **binding** (= B/2 − E).

Molecules:

```
geometry (bohr):  (0, 0, -0.12)  (0, 0, 0.12)            nuclei as used
  scf  3  eps = ...  E_phys = ...  dE = ...       E_phys = physical energy (spin included)
final UHF ... E_phys = ... converged = yes              converged = NO: maxscf reached, number not final
geo  3  R = ...  E_phys = ...  dE/dR (HF) = ...  |F_other| = ...  forces: (...) (...)
```

The search ends when |dE/dR| < gtol, when the step is smaller than rtol, at `geo: STOP` (rcut range), or after maxgeo
steps; the last `geo` line has the final R, energy and force — check that |dE/dR| < gtol there.


## 6. Memory and time of test runs (19 threads, 30 GB machine) — for planning only

| Calculation | Fourier modes | Memory | Time |
|---|---|---|---|
| install test (H, B = 100 B0, thresh 1e-5, 2 iterations) | 217 | 4 GB | 4–9 min |
| H, B = 100 B0, thresh 1e-7, qextra 10 … 80 | 217–373 | 9–20 GB | 25–40 min |
| H, B = 425.5 B0, thresh 1e-7, rcut 0.65 | 207 | ~10 GB | 17–30 min |
| H, thresh 1e-7, 450+ modes | 445–681 | more than 30 GB | — |
| He atom, B = 425.5 B0, thresh 1e-6, rcut 0.65 | 207 | ~15 GB | 80 min |
| H2 triplet, B = 425.5 B0, thresh 1e-6, rcut 0.65, 4 bond lengths | 207 | ~20 GB | ~4 h |

Cost of one SCF cycle ∝ (number of electrons) × (Fourier modes). Modes = 2J + 1 with J = ⌈qmax·P/(2π)⌉,
qmax = B·rext/2 + qextra, period P = 4·rext, rext = rcut + 3·rtaper: at strong field modes grow like B·rext²
(B = 425.5 B0: rcut 0.65 → 207, 0.75 → 271, 1.0 → 477 modes, all measured). The watchdog protects the machine (section 3e).


## 7. Known problems

- One B = 100 B0 test run hung after iteration 2. If a log does not grow for much longer than one
  iteration takes, run `bash stop.sh` and start the queue again (the hung job starts from the beginning).
- `E0` (or any `E0list` entry) above B/2 → the program stops immediately.
- A nucleus farther than rcut − 0.70 from the B axis → results not trustworthy. The job files respect this.
- Several electrons converge slowly (plain Green-function iteration, no accelerators): 10–30 SCF cycles.
  `NOT-CONVERGED` → rerun that job with larger `maxscf` (it starts from scratch).
