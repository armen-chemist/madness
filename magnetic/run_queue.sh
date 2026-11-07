#!/bin/bash
# Run the jobs of one job file, one after another, with a memory watchdog.
#
# Usage:   nohup bash run_queue.sh jobs/01_H_B100_qextra.jobs > queue.out 2>&1 &
#   (nohup: keeps running after you log out.  Stop:  kill <PID of run_queue.sh> and of magneticatom)
#
# Job file format:
#   COMMON <arguments used by every job of the file>
#   <name> <arguments of this job>        (later arguments override COMMON ones)
#   RESULTS <folder>                      (optional: one-calculation file made by tools/split_jobs.sh;
#                                          the log goes to results/<folder>/, which DONE then covers as a whole)
#   # comment lines and empty lines are ignored
# Output: results/<jobfile name>/<name>.log, and DONE in that folder when all jobs are finished.
# A job that ended normally gets the line "# JOB FINISHED" at the end of its log. Those jobs are skipped
# when the queue is started again; all others (stopped, killed, crashed) are run again from the start.

MAGNETICATOM=${MAGNETICATOM:-$(cd "$(dirname "$0")/.." && pwd)/build/src/madness/chem/magneticatom}
MINMEM_GB=${MINMEM_GB:-3}      # kill a job when free memory drops below this (GB)
STOP_AFTER_H=${STOP_AFTER_H:-} # start no new job after this many hours (the rest is left for a resubmission)
T0=$(date +%s)
export OMP_NUM_THREADS=1
# MADNESS threads: default = all cores. Set e.g. MAD_NUM_THREADS=24 to leave cores free.
# CORES=0-11 (taskset list): run every job of this queue on these cores only, MAD_NUM_THREADS = their number.
# Keep each queue inside one NUMA node (lscpu: node0 0-23, node1 24-47): 4 queues x 12 cores were fastest.
PIN=""
if [ -n "$CORES" ]; then
    PIN="taskset -c $CORES"
    NC=$(env -u OMP_NUM_THREADS taskset -c "$CORES" nproc) || { echo "bad CORES=$CORES"; exit 1; }   # nproc obeys OMP_NUM_THREADS
    export MAD_NUM_THREADS=${MAD_NUM_THREADS:-$NC}
fi

JOBS=$(readlink -f "$1")
[ -f "$JOBS" ] || { echo "usage: bash run_queue.sh <jobfile>"; exit 1; }
[ -x "$MAGNETICATOM" ] || { echo "program not found: $MAGNETICATOM (set MAGNETICATOM=...)"; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd)
RESULTS=$(awk '$1=="RESULTS"{print $2}' "$JOBS")
OUT="$HERE/results/${RESULTS:-$(basename "$JOBS" .jobs)}"
mkdir -p "$OUT" && cd "$OUT"
[ -z "$RESULTS" ] && rm -f DONE
COMMON=$(awk '$1=="COMMON"{$1=""; print}' "$JOBS")

grep -v -e '^#' -e '^COMMON' -e '^RESULTS' -e '^AFTER' -e '^RESTART' -e '^[[:space:]]*$' "$JOBS" | while read -r NAME ARGS; do
    LOG="$NAME.log"
    if [ -f "$LOG" ] && grep -q "^# JOB FINISHED" "$LOG"; then echo "skip $NAME (finished)"; continue; fi
    if [ -n "$STOP_AFTER_H" ] && [ $(( $(date +%s) - T0 )) -ge $(( STOP_AFTER_H * 3600 )) ]; then
        echo "$(date '+%F %T')  stop: $STOP_AFTER_H h reached, $NAME and later jobs are left for the next submission"; exit 3
    fi
    echo "$(date '+%F %T')  start $NAME"
    # an unfinished log of an earlier attempt is kept (renamed), never overwritten
    [ -f "$LOG" ] && mv "$LOG" "$LOG.unfinished_$(date +%m%d_%H%M)"
    echo "# $MAGNETICATOM $COMMON $ARGS" > "$LOG"
    echo "# started $(date '+%F %T') on $(hostname), MAD_NUM_THREADS=${MAD_NUM_THREADS:-all}, cores ${CORES:-all}" >> "$LOG"
    eval $PIN stdbuf -oL "$MAGNETICATOM" $COMMON $ARGS >> "$LOG" 2>&1 &
    PID=$!
    while kill -0 $PID 2>/dev/null; do
        AV=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
        if [ "$AV" -lt "$MINMEM_GB" ]; then kill $PID; echo "KILLED: low memory" >> "$LOG"; fi
        sleep 10
    done
    wait $PID 2>/dev/null; RC=$?
    if [ $RC -eq 0 ] && ! grep -q KILLED "$LOG"; then echo "# JOB FINISHED" >> "$LOG"; STATUS=finished
    else STATUS="NOT finished (exit code $RC)"; fi
    echo "$(date '+%F %T')  end   $NAME  $STATUS  $(grep -E '^final|KILLED' "$LOG" | tail -1)"
done
[ "${PIPESTATUS[1]}" = 3 ] && exit 0     # stopped by STOP_AFTER_H: not DONE
if [ -n "$RESULTS" ]; then
    # one-calculation file: the folder is DONE when every calculation of the parent file has finished
    for N in $(grep -v -e '^#' -e '^COMMON' -e '^[[:space:]]*$' "$HERE/jobs/$RESULTS.jobs" | awk '{print $1}'); do
        grep -q "^# JOB FINISHED" "$N.log" 2>/dev/null || exit 0
    done
fi
date > DONE
