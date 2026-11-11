#!/bin/bash
# First use case: every product for one North Pacific winter (November to March), 0-50 m.
#
#   scripts/run_winter.sh fetch 2014    # download sources for November 2014 - March 2015
#   scripts/run_winter.sh build 2014    # build U1-U6 and the checks (after fetch)
#
# Normally run through scripts/fetch.sbatch and scripts/build.sbatch (docs/SERVER_PLAN.md), with the
# ocean-uncertainty env active. Storage locations: scripts/paths.sh. Every step skips outputs that already exist, so
# an interrupted task is resumed by submitting it again.
# Settings and the evaluation that chose the U6 error model: docs/EVALUATION_2014-15.md.
set -eo pipefail   # no -u: conda's shell functions reference unset variables
cd "$(dirname "$0")/.."
source scripts/paths.sh

STAGE=$1
YEAR=$2
WINTER=$YEAR-$(printf %02d $(( (YEAR + 1) % 100 )))   # e.g. 2014-15
START=$YEAR-11-01
END=$((YEAR + 1))-03-31
REGION=(140 240 20 60)
MAX_DEPTH=50
DATA=$DATA_ROOT/$WINTER
OUT=$PRODUCT_ROOT/$WINTER
WOA=$DATA_ROOT/woa23

if [[ $STAGE == fetch ]]; then
  F=(python -m ocean_uncertainty.fetch)
  # Climatology months either side of November-March, for interpolation between mid-months. Shared by all winters:
  # the lock lets one array task download it while the others wait, then find the files and skip them.
  mkdir -p "$DATA_ROOT"
  flock "$DATA_ROOT/woa23.lock" "${F[@]}" woa23 --out "$WOA" --region "${REGION[@]}" --max-depth $MAX_DEPTH \
      --months 10 11 12 1 2 3 4
  "${F[@]}" argo --out "$DATA/argo" --region "${REGION[@]}" --start $START --end $END --max-pressure 100
  # Representativeness (U5) from one 1/12-degree day per month: the 15th.
  for month in $YEAR-11 $YEAR-12 $((YEAR + 1))-01 $((YEAR + 1))-02 $((YEAR + 1))-03; do
    "${F[@]}" glorys12 --out "$DATA/glorys12" --region "${REGION[@]}" --start $month-15 --end $month-15 --daily \
        --max-depth $MAX_DEPTH
  done
  "${F[@]}" ostia --out "$DATA/ostia" --region "${REGION[@]}" --start $START --end $END
  "${F[@]}" grep --out "$DATA/grep" --region "${REGION[@]}" --start $START --end $END --max-depth $MAX_DEPTH \
      --variables thetao mlotst
  touch "$DATA/FETCHED"

elif [[ $STAGE == build ]]; then
  [[ -e $DATA/FETCHED ]] || { echo "fetch $WINTER first" >&2; exit 1; }
  P=(python -m ocean_uncertainty.products)
  C=(python -m ocean_uncertainty.checks)
  STATS=(--max-distance-km 1000 --bin-km 25 --max-days-apart 3 --min-pairs 100 --map-bin-deg 2)
  GRID=$(ls "$DATA"/grep/*.nc | head -1)

  "${P[@]}" u1 --grep "$DATA"/grep/*.nc --out-dir "$OUT/u1"   # skips months already built
  [[ -e $OUT/u2.nc ]] || "${P[@]}" innovations --argo "$DATA"/argo/*.nc --background grep_mean \
      --background-files "$DATA"/grep/*.nc --out "$OUT/u2.nc" "${STATS[@]}"
  [[ -e $OUT/u3.nc ]] || "${P[@]}" innovations --argo "$DATA"/argo/*.nc --background woa23 \
      --background-files "$WOA"/*.nc --out "$OUT/u3.nc" "${STATS[@]}"
  [[ -e $OUT/u4_mixed_layer.nc ]] || "${P[@]}" u4 --ostia "$DATA"/ostia/*.nc --grep "$DATA"/grep/*.nc \
      --out "$OUT/u4.nc"
  [[ -e $OUT/u5.nc ]] || "${P[@]}" u5 --glorys12 "$DATA"/glorys12/*.nc --grid "$GRID" --out "$OUT/u5.nc"

  [[ -e $OUT/checks.md ]] || "${C[@]}" report --u1 "$OUT"/u1/*.nc --u2 "$OUT/u2.nc" --u3 "$OUT/u3.nc" \
      --u4 "$OUT/u4.nc" --u5 "$OUT/u5.nc" --pattern-depths 10 40 --min-count 10 --out "$OUT/checks.md"
  # Non-overlapping 7-day windows, centred from 4 November: the same windows as the 2014-15 evaluation.
  [[ -e $OUT/holdout.md ]] || "${C[@]}" holdout --u1 "$OUT"/u1/*.nc --u2 "$OUT/u2.nc" --u5 "$OUT/u5.nc" \
      --start $YEAR-11-04 --end $((YEAR + 1))-03-28 --window-days 7 --folds 5 --distance-edges-km 0 100 5000 \
      --sigma-classes 3 --out "$OUT/holdout.md"

  centres=$(python -c "import pandas as pd; print(*pd.date_range('$YEAR-11-04', '$((YEAR + 1))-03-28', freq='7D').date)")
  for centre in $centres; do
    month=${centre:0:4}${centre:5:2}
    u1=$(ls "$OUT"/u1/*_${month}01_*.nc)
    target=$OUT/u6/u6_${centre}.nc
    [[ -e $target ]] || "${P[@]}" u6 --u2 "$OUT/u2.nc" --u1 "$u1" --u5 "$OUT/u5.nc" --date $centre --window-days 7 \
        --background-model separable --obs-model u5_scaled --block-rows 4000 --out "$target"
  done
  touch "$OUT/BUILT"
else
  echo "usage: $0 fetch|build YEAR" >&2
  exit 2
fi
