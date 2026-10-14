#!/bin/bash
# All results at once. Molecule folders (11-20) get their dissociation energy D = E(atoms) - E(molecule)
# automatically, with E(atoms) taken from results/10_atoms_B425 at the SAME rcut as the molecule.
# Usage:  bash summarize_all.sh            (or:  bash summarize_all.sh > all_results.txt)
HERE=$(cd "$(dirname "$0")" && pwd)
R="$HERE/results"
eatom() {   # eatom <H|He> <rcut>  -> E_phys of the separated atom (empty if not finished)
    local L="$R/10_atoms_B425/$1_atom_rcut$2.log"
    [ -f "$L" ] || return
    if [ "$1" = H ]; then grep '^final E' "$L" | sed -n 's/.*B\/2 - E = \([-0-9.]*\).*/-\1/p'
    else grep '^final UHF' "$L" | tail -1 | sed -n 's/.*incl\.) = *\([-0-9.]*\).*/\1/p'; fi
}
declare -A EXACT=([01]=3.789804236305 [02]=3.789804236305 [05]=6.257087674681 [06]=0.831168896733)
for D in $(ls -d "$R"/*/ 2>/dev/null | sort); do
    F=$(basename "$D"); NUM=${F:0:2}; REF=${EXACT[$NUM]:-}
    case $NUM in
        11|15|19) RC=0.75 ;;
        *)        RC=1.0 ;;
    esac
    case $NUM in
        11|12|13|14) A=$(eatom H $RC);  [ -n "$A" ] && REF=$(awk -v a="$A" 'BEGIN{printf "%.10f",2*a}') ;;
        15|16|17|18) A=$(eatom H $RC);  [ -n "$A" ] && REF=$(awk -v a="$A" 'BEGIN{printf "%.10f",3*a}') ;;
        19|20)       A=$(eatom He $RC); [ -n "$A" ] && REF=$(awk -v a="$A" 'BEGIN{printf "%.10f",2*a}') ;;
    esac
    echo "=================== $F   $( [ -f "$D/DONE" ] && echo "(DONE)" || echo "(not finished)")   REF = ${REF:-none}"
    bash "$HERE/summarize.sh" "$D" $REF
    echo
done
