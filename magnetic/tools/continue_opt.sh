#!/bin/bash
# Write the continuation of an optimization (ropt=1) that ended without "# JOB FINISHED" (walltime kill):
#   bash tools/continue_opt.sh <job file> <its log> <new job file> <new calculation name> <geo steps done in the chain> [<new rcut>]
# With <new rcut> (after a "geo: STOP ... off the B axis"), rcut and rtaper = rcut/10 are replaced, so the nuclei may go
# further off the axis; the first step recomputes the last R with the wider cylinder (rcut check).
# The new file has the same settings, the nuclei scaled about their centroid to the last "geo: new R" of the log
# (as magneticatom.cc does: R = |R_1 - R_0|, shape fixed) and maxgeo = original maxgeo - steps done.
# Exit 1 (no file written, reason on stdout) when there is nothing to continue from or no steps are left.
JOB=$1; LOG=$2; NEW=$3; NEWNAME=$4; DONE=$5; WIDE=$6
LINE0=$(grep -v -e '^#' -e '^COMMON' -e '^RESULTS' -e '^AFTER' -e '^RESTART' -e '^[[:space:]]*$' "$JOB" | head -1)
if echo "$LINE0" | grep -q 'ropt=2'; then
    # ropt=2 (free BFGS on all coordinates, 2026-11-10): restart from the last "geo: new atoms=" (Hessian restarts),
    # steps done = "geo2" lines of this log; the "# PROGRAM" line is copied
    ANEW=$(awk '/^geo: new atoms=/{a=$3} END{print a}' "$LOG" | sed 's/^atoms=//')
    [ -n "$ANEW" ] || { echo "no geometry step finished before the kill: needs a longer walltime"; exit 1; }
    DONE=$(grep -c '^geo2 ' "$LOG")
    MAXGEO=$(echo "$LINE0 $(grep '^COMMON' "$JOB")" | grep -oE '(^| )maxgeo=[0-9]+' | head -1 | cut -d= -f2)
    LEFT=$(( ${MAXGEO:-12} - DONE ))
    [ $LEFT -gt 0 ] || { echo "all maxgeo=${MAXGEO:-12} geometry steps used"; exit 1; }
    mkdir -p "$(dirname "$NEW")"
    {
        echo "# continuation of $JOB (ropt=2), written by tools/continue_opt.sh $(date '+%F %T') after it ended unfinished:"
        echo "# start = last \"geo: new atoms=\" of $LOG; maxgeo $LEFT = ${MAXGEO:-12} - $DONE steps done"
        grep '^# PROGRAM' "$JOB"
        grep -e '^RESULTS' -e '^COMMON' "$JOB"
        echo "$LINE0" | sed -E "s/^[^ ]+/$NEWNAME/; s/maxgeo=[0-9]+/maxgeo=$LEFT/; s/\"atoms=[^\"]*\"/\"atoms=$ANEW\"/"
    } > "$NEW"
    echo "atoms $ANEW, maxgeo $LEFT"; exit 0
fi
RNEW=$(awk '/^geo: new R =/{r=$5} END{if (r!="") printf "%.6f", r}' "$LOG")
[ -n "$RNEW" ] || { echo "no geometry step finished before the kill: needs a longer walltime"; exit 1; }
LINE=$(grep -v -e '^#' -e '^COMMON' -e '^RESULTS' -e '^AFTER' -e '^RESTART' -e '^[[:space:]]*$' "$JOB" | head -1)
MAXGEO=$(echo "$LINE $(grep '^COMMON' "$JOB")" | grep -oE '(^| )maxgeo=[0-9]+' | head -1 | cut -d= -f2)
LEFT=$(( ${MAXGEO:-12} - DONE ))
[ $LEFT -gt 0 ] || { echo "all maxgeo=${MAXGEO:-12} geometry steps used (last new R $RNEW)"; exit 1; }
mkdir -p "$(dirname "$NEW")"
{
    if [ -n "$WIDE" ]; then
        echo "# continuation of $JOB, written by tools/continue_opt.sh $(date '+%F %T') after the off-axis guard stopped it:"
        echo "# rcut widened to $WIDE (rtaper $(awk -v r="$WIDE" 'BEGIN{print r/10}')); first step = rcut check at the last R"
    else
        echo "# continuation of $JOB, written by tools/continue_opt.sh $(date '+%F %T') after it ended unfinished:"
    fi
    echo "# start R $RNEW = last \"geo: new R\" of $LOG; maxgeo $LEFT = ${MAXGEO:-12} - $DONE steps done"
    grep -e '^RESULTS' -e '^COMMON' "$JOB" | if [ -n "$WIDE" ]; then
        sed -E "s/rcut=[0-9.]+/rcut=$WIDE/; s/rtaper=[0-9.]+/rtaper=$(awk -v r="$WIDE" 'BEGIN{print r/10}')/"; else cat; fi
    echo "$LINE" | sed -E "${WIDE:+s/ rcut=[0-9.]+//; s/ rtaper=[0-9.]+//}" | awk -v rn="$RNEW" -v left="$LEFT" -v name="$NEWNAME" '{
        $1 = name
        sub(/maxgeo=[0-9]+/, "maxgeo=" left); if ($0 !~ /maxgeo=/) $0 = $0 " maxgeo=" left
        i = index($0, "\"atoms="); head = substr($0, 1, i + 6); rest = substr($0, i + 7); j = index(rest, "\"")
        n = split(substr(rest, 1, j - 1), at, ";")
        for (k = 1; k <= n; k++) { split(at[k], c, ","); Z[k] = c[1]; for (a = 2; a <= 4; a++) { x[k, a] = c[a]; m[a] += c[a] / n } }
        R0 = sqrt((x[2,2] - x[1,2])^2 + (x[2,3] - x[1,3])^2 + (x[2,4] - x[1,4])^2); f = rn / R0; s = ""
        for (k = 1; k <= n; k++) s = s (k > 1 ? ";" : "") sprintf("%s,%.6f,%.6f,%.6f", Z[k], m[2] + (x[k,2] - m[2]) * f, m[3] + (x[k,3] - m[3]) * f, m[4] + (x[k,4] - m[4]) * f)
        print head s substr(rest, j) }'
} > "$NEW"
echo "R $RNEW, maxgeo $LEFT"
