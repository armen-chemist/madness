#!/bin/bash
# What is running now, memory, and the last lines of every running job.
echo "=== running programs"
ps -u "$USER" -o pid=,etime=,rss=,comm=,args= | awk '$4=="magneticatom" || ($4=="bash" && /run_queue\.sh/) {
    a=""; for (i=5; i<=NF; i++) a=a" "$i; if (length(a) > 110) a=substr(a,1,110)"...";
    printf "pid %-8s time %-11s %5.1f GB %s\n", $1, $2, $3/1048576, a }'
echo "=== memory (GB)"
free -g | head -2
echo "=== queues"
HERE=$(cd "$(dirname "$0")" && pwd)
for D in "$HERE"/results/*/; do
    [ -d "$D" ] || continue
    if [ -f "$D/DONE" ]; then S="DONE $(cat "$D/DONE")"; else S="running or stopped"; fi
    echo "$(basename "$D"):  $S"
    L=$(ls -t "$D"/*.log 2>/dev/null | head -1)
    [ -n "$L" ] && echo "   newest log $(basename "$L"):" && grep -E '^  iter|^  scf|^geo|^final|KILLED|JOB FINISHED' "$L" | tail -2 | cut -c1-150 | sed 's/^/     /'
done
