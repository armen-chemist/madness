#!/bin/bash
# Check that the hawaii build (shared home, conda env "madness") runs on the bra08 Slurm nodes: same test as install.sh.
# Usage (on bra08):   cd ~/madness/magnetic && sbatch slurm/test_bra08.sbatch
#SBATCH -J mad-test
#SBATCH -p cpu-14
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 00:30:00
#SBATCH -o slurm/test_bra08.out

source "$HOME/miniforge3/etc/profile.d/conda.sh" && conda activate madness
echo "host $(hostname), $(nproc) cores, $(date '+%F %T')"
EXE=$HOME/madness/build/src/madness/chem/magneticatom
D=$HOME/madness/magnetic/slurm/test_bra08 && mkdir -p "$D" && cd "$D"
/usr/bin/time -v "$EXE" B=100 E0=48.1 k=8 thresh=1e-5 L=6 rcut=1.34 rtaper=0.134 maxiter=2 > test.log 2> time.log
E1=$(awk '$1=="iter" && $2=="1" {print $5}' test.log)
echo "iteration 1: E = $E1   (expected: 46.230390629)"
awk -v e="$E1" 'BEGIN{exit !(e != "" && (e-46.230390629)^2 < 1e-10)}' && echo PASSED || echo FAILED
grep -E "Elapsed|Maximum resident" time.log
echo "end $(date '+%F %T')"
