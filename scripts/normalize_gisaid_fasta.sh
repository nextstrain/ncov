#!/usr/bin/env bash
set -e
GISAID_SARSCOV2_IN=$1
GISAID_SARSCOV2_OUT=$2
MIN_LENGTH=$3

if [[ ! -r "$GISAID_SARSCOV2_IN" ]]
then
	echo "$0: input $GISAID_SARSCOV2_IN not found"
	exit 1
fi

if [[ -z "$MIN_LENGTH" ]]
then
	echo "Using default minimum length of 25000"
	MIN_LENGTH=25000
fi

echo "Normalizing GISAID file $GISAID_SARSCOV2_IN to $GISAID_SARSCOV2_OUT (min length $MIN_LENGTH)"

# Remove leading virus name prefix from sequence names
# Remove embedded spaces in sequence names (Hong Kong sequences)
# Remove trailing |EPI_ISL_id|datestamp from sequence names
# Remove sequences shorter than minimum length
# Eliminate duplicate sequences (keep only the first seen)

#cat $GISAID_SARSCOV2_IN |
	sed 's/^>[hn]Co[Vv]-19\//>/g' $GISAID_SARSCOV2_IN |	# remove leading prefix
	sed 's/ //g' |					# remove embedded spaces
	sed 's/|.*$//' | 				# remove trailing metadata
	awk "BEGIN{RS=\">\";FS=\"\n\"}length>$MIN_LENGTH{print \">\"\$0}" |	# remove short seqs
	awk 'BEGIN{RS=">";FS="\n"}!x[$1]++{print ">"$0}' | 	# remove duplicates
	grep -v '^>*$' > $GISAID_SARSCOV2_OUT

exit 0
