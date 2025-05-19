#!/bin/bash

# parameters are in the order they are passed to vapor (see snakemake)
input_fw=$1
input_rev=$2
input_ref=$3
output=$4
benchmark_out=$5

if [ $# -ne 4]; then
	echo "Provide 4 arguments for vapor and one for the benchmark file"
	exit 1
fi

START_TIME=$(date +%s)

vapor.py -fq "$input_fw" "$input_rev" -fa "$input_ref" -o "$output" &
PID=$!

MEM_PEAK="0"
while ps | grep " $PID " > /dev/null; do
	MEM_PEAK=$(grep -e "VmHWM" /proc/$PID/status)
	sleep 1
done

END_TIME=$(date +%s)
ELLAPSED_TIME=$((END_TIME - START_TIME))

echo "$input_ref : $MEM_PEAK , $ELLAPSED_TIME" >> "$benchmark_out"