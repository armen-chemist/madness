#!/bin/bash
# Start the next unfinished job files of batches/order.list on free nodes, without getting in other users' way:
#   - nothing is submitted while a job of another user waits in the queue (theirs go first);
#   - nothing is submitted unless a node is free right now (ours never wait in the queue ahead of others);
#   - at most CAP of our jobs at the same time (default all 8 nodes; the other rules make room for others);
#   - every job stops starting new calculations after STOP_AFTER_H hours (default 12) and gives its node back;
#     the next call of this script continues it (finished calculations are skipped).
# Usage (on hawaii): cd ~/madness/magnetic && bash submit_next.sh      (DRY=1: only show what it would do)
# Runs every 5 minutes from cron on hawaii (crontab -l); log: batches/submit.log. Pause it: touch batches/PAUSE
cd "$(dirname "$0")"
exec 9> batches/.submit.lock; flock -n 9 || exit 0
# daily read-only snapshot of all logs, job files, scripts and paper sources (raw data safety net; not in git)
SNAP=$HOME/madness_backups/magnetic_$(date +%F).tar.gz
if ! ls "$HOME"/madness_backups/magnetic_$(date +%F)*.tar.gz >/dev/null 2>&1; then
    mkdir -p "$HOME/madness_backups" && tar czf "$SNAP" results jobs batches tools submit_next.sh run_queue.sh run_queue.pbs \
        docs/paper/*.tex docs/paper/*.py docs/paper/*.bib docs/paper/accepted.txt docs/paper/build.sh 2>/dev/null; \
        [ -s "$SNAP" ] && chmod a-w "$SNAP" && echo "$(date '+%F %T')  snapshot $SNAP"
fi
[ -f batches/PAUSE ] && exit 0
CAP=${CAP:-8}
STOP_AFTER_H=${STOP_AFTER_H:-12}
WIDE_RCUT=${WIDE_RCUT:-1.34}   # rcut for an optimization stopped by the off-axis guard (2026-10-24, user request)
ME=$(whoami)

# finished optimizations of an *opt file count for its scan file too (11opt -> 11), so they are not run twice
for D in results/*opt_*; do
    S=results/$(basename "$D" | sed 's/opt_/_/')
    [ -d "$S" ] || continue
    for LOG in "$D"/ropt_*.log; do
        [ -f "$LOG" ] && grep -q "^# JOB FINISHED" "$LOG" && [ ! -f "$S/$(basename "$LOG")" ] && cp "$LOG" "$S/" && echo "copied $LOG to $S/"
    done
done

others_waiting=$(qstat -t 2>/dev/null | awk -v me="$ME" 'NR>2 && $3!=me && $5=="Q"' | wc -l)
mine=$(qstat -t -u "$ME" 2>/dev/null | awk 'NR>5 && ($10=="R" || $10=="Q")' | wc -l)
free_nodes=$(pbsnodes -a 2>/dev/null | grep -c '^ *state = free$')
echo "$(date '+%F %T')  ours $mine (cap $CAP), free nodes $free_nodes, other users' jobs waiting $others_waiting"
[ "$others_waiting" -eq 0 ] || { echo "other users are waiting: submit nothing"; exit 0; }

# job id file -> is that Torque job still queued, running or held?
running() { [ -f "$1" ] && qstat -t "$(cat "$1")" 2>/dev/null | awk 'NR>2{print $5}' | grep -q '[RQH]'; }
# first calculation name of a job file
calcname() { grep -v -e '^#' -e '^COMMON' -e '^RESULTS' -e '^AFTER' -e '^RESTART' -e '^[[:space:]]*$' "$1" | awk 'NR==1{print $1}'; }

slots=$(( CAP - mine )); [ $free_nodes -lt $slots ] && slots=$free_nodes
grep -v -e '^#' -e '^[[:space:]]*$' batches/order.list | while read -r JOBFILE WALL; do
    [ $slots -gt 0 ] || break
    RES=$(awk '$1=="RESULTS"{print $2}' "$JOBFILE")
    if [ -n "$RES" ]; then
        # one-calculation file (tools/split_jobs.sh): log in results/<RES>/, Torque id in results/<RES>/.jobid/
        OUT=results/$RES; BASE=$(basename "$JOBFILE" .jobs); NAME0=$(calcname "$JOBFILE")
        # an optimization killed unfinished is continued in a chain jobs/split/auto/<BASE>_c<K>.jobs -> <NAME0>_c<K>.log;
        # follow the chain to its last member: CUR (job file), CLOG (its log)
        K=0; CUR=$JOBFILE; CLOG=$OUT/$NAME0.log
        while [ -f "jobs/split/auto/${BASE}_c$((K+1)).jobs" ]; do
            K=$((K+1)); CUR=jobs/split/auto/${BASE}_c$K.jobs; CLOG=$OUT/${NAME0}_c$K.log
        done
        JIDF=$OUT/.jobid/$(basename "$CUR" .jobs)
        # an optimization stopped by the off-axis guard (nucleus beyond rcut - 0.70) continues once with rcut WIDE_RCUT
        WIDE=""
        if [ -f "$CLOG" ] && grep -q "^# JOB FINISHED" "$CLOG"; then
            grep -q "^geo: STOP.*off the B axis" "$CLOG" && grep -q 'ropt=1' "$CUR" \
                && awk -v w="$WIDE_RCUT" '/^COMMON/{for(i=1;i<=NF;i++) if($i~/^rcut=/){split($i,a,"=");exit !(a[2]<w)}}' "$CUR" \
                && WIDE=$WIDE_RCUT || continue
        fi
        [ -f "$OUT/.jobid/$BASE.stopped" ] && continue           # chain ended without converging: see that file
        running "$JIDF" && continue
        if [ -f "$CLOG" ] && [ $K -eq 0 ]; then
            # unfinished log of an old whole-file job: wait while that job runs (no JOBID: while written in the last 4 h)
            if [ -f "$OUT/JOBID" ]; then running "$OUT/JOBID" && continue
            else [ -n "$(find "$CLOG" -mmin -240)" ] && continue; fi
        fi
        if [ -f "$CLOG" ] && grep -q 'ropt=1' "$CUR"; then
            # the optimization ended unfinished (walltime) or at the off-axis guard (WIDE): continue it from its last R
            NEXT=jobs/split/auto/${BASE}_c$((K+1)).jobs
            STEPS=$(grep -c "^geo " "$CLOG")   # steps of the last member only: its maxgeo already excludes the earlier ones
            if [ -n "$DRY" ]; then echo "would continue $CUR as $NEXT ($STEPS geo steps done${WIDE:+, rcut $WIDE})"; slots=$(( slots - 1 )); continue; fi
            if MSG=$(bash tools/continue_opt.sh "$CUR" "$CLOG" "$NEXT" "${NAME0}_c$((K+1))" "$STEPS" $WIDE); then
                echo "continue $CUR as $NEXT: $MSG"; CUR=$NEXT; JIDF=$OUT/.jobid/$(basename "$NEXT" .jobs)
            else
                mkdir -p "$OUT/.jobid"; echo "$(date '+%F %T')  $CUR ($CLOG): not continued: $MSG" | tee "$OUT/.jobid/$BASE.stopped"; continue
            fi
        fi
    else
        OUT=results/$(basename "$JOBFILE" .jobs); JIDF=$OUT/JOBID; CUR=$JOBFILE
        [ -f "$OUT/DONE" ] && continue
        # JIDF: the Torque job that runs (or ran) this file; skip it while that job is still queued or running
        running "$JIDF" && continue
    fi
    NAME=q$(basename "$CUR" | cut -d_ -f1)$([ -n "$RES" ] && basename "$CUR" .jobs | sed 's/.*__/_/' | cut -c1-9)
    if [ -n "$DRY" ]; then echo "would submit $CUR (walltime $WALL h)"
    else
        ID=$(qsub -N "$NAME" -l walltime=$WALL:00:00 -v JOBFILE="$CUR",STOP_AFTER_H="$STOP_AFTER_H" run_queue.pbs) || continue
        mkdir -p "$(dirname "$JIDF")" && echo "$ID" > "$JIDF"
        echo "submitted $CUR as $ID (walltime $WALL h)"
    fi
    slots=$(( slots - 1 ))
done
