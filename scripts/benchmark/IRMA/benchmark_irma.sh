#!/bin/bash
# creates a subfolder for the sample name, provided as an cl argument. 
# a "bench_out.txt" file is created, including 1. the highest RAM peak (measured every second) and 2. the total runtime.
# the used IRMA version was stopped after read gathering phase, with a simple exit statement in the code, as described in the master thesis.
# for reuse, adapt folder and read names.

sample_name=$1

START_TIME=$(date +%s)

mkdir "$sample_name" &
IRMA FLU /home/galaxy/akolbecher/workdir/data/benchmark_samples_TEST/"$sample_name"/forward.fastqsanger.gz /home/galaxy/akolbecher/workdir/data/benchmark_samples_TEST/"$sample_name"/reverse.fastqsanger.gz "$sample_name"/irmaOut &
PID=$!

MEM_PEAK="0"
while kill -0 $PID 2>/dev/null; do
	ALL_PIDS=$(pstree -p $PID | grep -o '[0-9]\+')

	CURR_MEM=0
	for pid in $ALL_PIDS; do
		RSS=$(ps -o rss= -p $pid 2>/dev/null)
		CURR_MEM=$((CURR_MEM + RSS))
	done
	if [ "$CURR_MEM" -gt "$MEM_PEAK" ]; then
		MEM_PEAK=$CURR_MEM
	fi

	sleep 1
done

END_TIME=$(date +%s)
ELLAPSED_TIME=$((END_TIME - START_TIME))

echo "$MEM_PEAK, $ELLAPSED_TIME" > "$sample_name"/bench_out.txt

