#!/bin/bash
# Build MADNESS + the magnetic-field code on a Linux server.
#
# Usage:   bash install.sh [BUILD_DIR] [NPROC]      (from the folder magnetic/ of this repository)
#            BUILD_DIR    build folder              (default: <repository>/build)
#            NPROC        cores used for compiling  (default: all)
#
# What it does:
#   1. Miniforge (if no conda) + conda environment "madness" (gcc 14, OpenMPI, OpenBLAS, cmake 3) [needs internet]
#      -> skip with USE_CONDA=0 if the server provides compilers/MPI/BLAS itself (see README, section 2b)
#   2. cmake + make magneticatom                                                             [~30-60 min]
#   3. 4-minute test run, prints PASSED or FAILED
#
# Result:  BUILD_DIR/src/madness/chem/magneticatom
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
NPROC=${2:-$(nproc)}
USE_CONDA=${USE_CONDA:-1}
ENV_NAME=${ENV_NAME:-madness}

[ -f "$HERE/../src/madness/chem/magneticatom.cc" ] \
    || { echo "ERROR: install.sh must stay in the folder magnetic/ of the repository"; exit 1; }

echo "=== 1. compilers, MPI, BLAS"
if [ "$USE_CONDA" = 1 ]; then
    if ! command -v conda >/dev/null; then
        if [ ! -d "$HOME/miniforge3" ]; then
            URL=https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
            wget -q "$URL" -O "$HOME/miniforge.sh" 2>/dev/null || curl -sL "$URL" -o "$HOME/miniforge.sh" \
                || { echo "ERROR: cannot download Miniforge (no internet?) -> see README section 2b"; exit 1; }
            bash "$HOME/miniforge.sh" -b -p "$HOME/miniforge3"
            rm "$HOME/miniforge.sh"
        fi
        source "$HOME/miniforge3/etc/profile.d/conda.sh"
    else
        source "$(conda info --base)/etc/profile.d/conda.sh"
    fi
    conda env list | grep -q "^$ENV_NAME " || conda env create -n "$ENV_NAME" -f "$HERE/environment.yml"
    conda activate "$ENV_NAME"
    LAPACK_LIBS="$CONDA_PREFIX/lib/liblapack.so;$CONDA_PREFIX/lib/libblas.so"
else
    # system compilers: mpicc/mpicxx must be in PATH; LAPACK_LIBS must point to LAPACK + BLAS
    : "${LAPACK_LIBS:?USE_CONDA=0 needs LAPACK_LIBS, e.g. LAPACK_LIBS=/usr/lib64/liblapack.so\;/usr/lib64/libopenblas.so}"
fi
command -v mpicxx >/dev/null || { echo "ERROR: mpicxx not found"; exit 1; }
echo "compiler: $(which mpicxx)"

REPO=$(cd "$HERE/.." && pwd)
DIR=${1:-$REPO/build}

echo "=== 2. configure and build (log: $DIR/build.log)"
mkdir -p "$DIR" && cd "$DIR"
cmake "$REPO" -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_C_COMPILER="$(which mpicc)" -DCMAKE_CXX_COMPILER="$(which mpicxx)" \
         -DENABLE_MPI=ON -DMADNESS_TASK_BACKEND=Pthreads -DENABLE_MKL=OFF -DENABLE_PYTHON=OFF \
         -DLAPACK_LIBRARIES="$LAPACK_LIBS" -DBUILD_TESTING=ON > cmake.log 2>&1 \
    || { echo "ERROR: cmake failed, see $DIR/cmake.log"; exit 1; }
make -j"$NPROC" magneticatom > build.log 2>&1 \
    || { echo "ERROR: build failed, see $DIR/build.log"; exit 1; }
EXE="$DIR/src/madness/chem/magneticatom"
echo "built: $EXE"

echo "=== 3. test (H atom, B = 100 B0, 2 iterations, about 4-10 min)"
cd "$DIR"
"$EXE" B=100 E0=48.1 k=8 thresh=1e-5 L=6 rcut=1.34 rtaper=0.134 maxiter=2 > test.log 2>&1 || true
E1=$(awk '$1=="iter" && $2=="1" {print $5}' test.log)
echo "iteration 1: E = $E1   (expected: 46.230390629)"
if awk -v e="$E1" 'BEGIN{exit !(e != "" && (e-46.230390629)^2 < 1e-10)}'; then
    echo "PASSED"
else
    echo "FAILED: see $DIR/test.log"; exit 1
fi
