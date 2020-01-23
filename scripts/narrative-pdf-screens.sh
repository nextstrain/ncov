#!/usr/bin/env bash
set -e
ADDRESS=$1
PARTIAL_FNAME=$2

echo "Make sure that you have auspice / nextstrain running at so that ${ADDRESS} is valid!"
echo "(e.g. from the 'ncov' directory run 'auspice view --datasetDir auspice --narrativeDir narratives'"
echo "This script will save PDFs starting with the prefix ${PARTIAL_FNAME}"
echo ""

# https://gs.statcounter.com/screen-resolution-stats/desktop/worldwide
RESOLUTIONS=(3200x1350 1920x1080 1600x900 366x768 )
#             james'     iphone+            iphone

for RES in ${RESOLUTIONS[@]}; do
  F="${PARTIAL_FNAME}.${RES}.pdf"
  echo ""
  echo "-------------------------------------"
  echo "Making ${F}"
  echo ""
  decktape generic --load-pause 3000 --key ArrowDown --size ${RES} ${ADDRESS} ${F}
done

exit 0
