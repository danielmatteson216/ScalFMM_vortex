#!/bin/bash

convert_to_csv()
{
	#stop_mark=`grep "stop_profiling" $1 | cut -f 2`
	#start_mark=`grep "start_profiling" $1 | cut -f 2`
	grep -e '^\(\(%\)\|\(\(0\|1\|2\|3\|4\|5\|6\|7\|9\)\>\)\)' $1 > start.trace
	grep -e '^\(\(%\)\|\(\(0\|1\|2\|3\|4\|5\|6\|7\|9\|18\|19\)\>\)\)' -v  $1 > end.trace
	sort -s -V --key=2,2 end.trace > endSorted.trace
	cat start.trace endSorted.trace > outputSorted.trace
	egrep -v "[[:digit:]]_t(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24)[[:space:]]" outputSorted.trace > sorted.trace
	#echo "$start_mark - $stop_mark"
	# Converting to .csv
	#pj_dump --start=$start_mark --end=$stop_mark -n -u outputSorted.trace > $2
	pj_dump -n -u sorted.trace > $2
	# Keeping only the states
	#perl -i -ne 'print if /^State/' $2
	# Delete temporary files
	rm outputSorted.trace start.trace end.trace endSorted.trace sorted.trace
	# Adding a header manually
	sed -i '1s/^/Nature, ResourceId, Type, Start, End, Duration, Depth, Value, Footprint, JobId, Params, Size, Tag\n/' $2 
}
directory=$1

cd $directory
list_fxt_file=`ls fxt`
args=""
gantt_args=""

for j in $list_fxt_file; do
	args="$args fxt/$j"
done

#Aggregate trace
starpu_fxt_tool -i $args
echo "Generate CSV ..."
convert_to_csv paje.trace paje.csv
cd ..

mkdir only_gantt 2>/dev/null
mkdir only_gantt/output 2>/dev/null
cd only_gantt

#Run a python script
echo "Run python ..."
../loutre.py -t ../$directory/trace.rec -i ../$directory/stdout -g gantt_only.db -o tmp.db
count=$(($count+1))
rm -f tmp.db

Rscript ../only_gantt.R
rm gantt_only.db
