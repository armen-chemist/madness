#!/bin/bash
# Split a job file into one-calculation job files, so that submit_next.sh can spread them over free nodes.
# Usage:  bash tools/split_jobs.sh <jobs/parent.jobs> <walltime h> [name ...]
#   names given: only those calculations; none given: every calculation of the parent.
#   Writes jobs/split/<parent>__<name>.jobs and prints "<chunk file> <walltime>" lines for batches/order.list.
# A chunk contains "RESULTS <parent>", so its log goes to results/<parent>/ like before (run_queue.sh).
cd "$(dirname "$0")/.."
P=${1:?usage: split_jobs.sh <jobs/parent.jobs> <walltime h> [name ...]}; WALL=${2:?walltime}; shift 2
BASE=$(basename "$P" .jobs)
mkdir -p jobs/split
NAMES=${*:-$(grep -v -e '^#' -e '^COMMON' -e '^[[:space:]]*$' "$P" | awk '{print $1}')}
for N in $NAMES; do
    LINE=$(awk -v n="$N" '$1==n' "$P")
    [ -n "$LINE" ] || { echo "no calculation $N in $P" >&2; exit 1; }
    C=jobs/split/${BASE}__$N.jobs
    { echo "# one calculation of $P, made by tools/split_jobs.sh"
      echo "RESULTS $BASE"
      grep '^COMMON' "$P"
      echo "$LINE"; } > "$C"
    echo "$C $WALL"
done
