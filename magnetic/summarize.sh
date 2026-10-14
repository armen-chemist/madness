#!/bin/bash
# Table of all logs in one results folder.
#
# Usage:  bash summarize.sh results/<folder>  [REF]
#
# One-electron jobs (H atom): REF = exact binding (Ha). Columns: binding = B/2 - E, err = binding - REF (uHa).
# Molecules / He atom:        REF = energy of the separated atoms (Ha), e.g. 2 E_phys(H) from file 10.
#                             Columns: E_phys (spin Zeeman included), D = REF - E_phys in Ha and eV
#                             (D > 0: bound with respect to the atoms). For ropt jobs: last R, dE/dR / |F_other|;
#                             R(STOP) = search stopped because a nucleus would leave the rcut region.
# status: ok = converged | NOT-CONVERGED = stopped at maxiter/maxscf | running/stopped | KILLED ...
DIR=${1:?usage: bash summarize.sh <results folder> [REF]}
REF=${2:-}
HA_EV=27.211386
H1=0; H2=0
hdr() {
    if [ "$1" = one ] && [ $H1 = 0 ]; then H1=1
        printf "%-22s %6s %5s %16s %11s %10s %8s  %s\n" job modes iter "binding(Ha)" "<H>-E" "err(uHa)" "wall(s)" status; fi
    if [ "$1" = uhf ] && [ $H2 = 0 ]; then H2=1
        printf "%-22s %6s %5s %16s %12s %9s %10s %20s %8s  %s\n" job modes scf "E_phys(Ha)" "D(Ha)" "D(eV)" "R(ropt)" "dE/dR/F_oth" "wall(s)" status; fi
}
for PASS in one uhf; do
for LOG in $(ls -v "$DIR"/*.log 2>/dev/null); do
    if grep -q -e '^final E' -e 'nelec=1 ' "$LOG" || ! grep -q -e 'nelec=' -e 'UHF' "$LOG"; then T=one; else T=uhf; fi
    [ $T = $PASS ] || continue
    hdr $T
    N=$(basename "$LOG" .log)
    MODES=$(grep -m1 -o 'modes = [0-9]*' "$LOG" | awk '{print $3}')
    ST="running/stopped"
    if [ $T = one ]; then                                      # one electron
        FIN=$(grep '^final E' "$LOG" | tail -1)
        BIND=$(echo "$FIN" | sed -n 's/.*B\/2 - E = \([-0-9.]*\).*/\1/p')
        WALL=$(echo "$FIN" | sed -n 's/.*wall = \([0-9.]*\).*/\1/p')
        IT=$(grep -c '^  iter ' "$LOG")
        HE=$(grep -m1 '<H> abgv' "$LOG" | sed -n 's/.*<H> - E = \([-+0-9.e]*\).*/\1/p')
        ERR=""; [ -n "$REF" ] && ERR=$(awk -v b="$BIND" -v e="$REF" 'BEGIN{printf "%+.2f",(b-e)*1e6}')
        [ -n "$FIN" ] && ST=ok; echo "$FIN" | grep -q "converged = NO" && ST=NOT-CONVERGED
        grep -q KILLED "$LOG" && ST=$(grep -m1 KILLED "$LOG")
        printf "%-22s %6s %5s %16s %11s %10s %8s  %s\n" "$N" "$MODES" "$IT" "$BIND" "$HE" "$ERR" "$WALL" "$ST"
    else                                                       # UHF (He atom, molecules)
        FIN=$(grep '^final UHF' "$LOG" | tail -1)
        IT=$(grep -c '^  scf ' "$LOG")
        EP=""; WALL=""; D=""; DEV=""; RR=""; G=""
        if [ -n "$FIN" ]; then
            EP=$(echo "$FIN" | sed -n 's/.*incl\.) = *\([-0-9.]*\).*/\1/p')
            WALL=$(echo "$FIN" | sed -n 's/.*wall = \([0-9.]*\).*/\1/p')
            ST=ok; echo "$FIN" | grep -q 'converged = NO' && ST=NOT-CONVERGED
        else
            EP=$(grep '^  scf ' "$LOG" | tail -1 | sed -n 's/.*E_phys = \([-0-9.]*\).*/\1/p')
            [ -n "$EP" ] && EP="($EP)"                         # last SCF value of an unfinished job
        fi
        if grep -q '^geo ' "$LOG" && grep -q 'ropt=1' <(head -1 "$LOG"); then
            RR=$(grep '^geo ' "$LOG" | tail -1 | sed -n 's/.*R = *\([0-9.]*\).*/\1/p')
            G=$(grep '^geo ' "$LOG" | tail -1 | sed -n 's/.*dE\/dR (HF) = *\([-+0-9.e]*\).*/\1/p')
            FO=$(grep '^geo ' "$LOG" | tail -1 | sed -n 's/.*|F_other| = *\([-+0-9.e]*\).*/\1/p')
            [ -n "$FO" ] && G="$G/$FO"
            grep -q '^geo: STOP' "$LOG" && RR="$RR(STOP)"
        fi
        E=${EP//[()]/}
        if [ -n "$REF" ] && [ -n "$E" ]; then
            D=$(awk -v r="$REF" -v e="$E" 'BEGIN{printf "%+.6f",r-e}')
            DEV=$(awk -v r="$REF" -v e="$E" -v c=$HA_EV 'BEGIN{printf "%+.3f",(r-e)*c}')
        fi
        grep -q KILLED "$LOG" && ST=$(grep -m1 KILLED "$LOG")
        grep -qE 'MadnessException|terminate called|Segmentation' "$LOG" && ST="ERROR (see log)"
        printf "%-22s %6s %5s %16s %12s %9s %10s %20s %8s  %s\n" "$N" "$MODES" "$IT" "$EP" "$D" "$DEV" "$RR" "$G" "$WALL" "$ST"
    fi
done
done
