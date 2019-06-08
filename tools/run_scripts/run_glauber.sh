#!/bin/bash

# Sample script for running glauber simulations

numCpus=4
beginRun=1
endRun=10
#tsteps=40
tsteps=10

t=0.7915589182064
tc=`echo "1/l(1+sqrt(5))" | bc -l`
#tmin=`echo "$tc-0.01" | bc -l`
tmin=`echo "$tc-0.002" | bc -l`
tmax=`echo "$tc+0.01" | bc -l`
t=$tmin

outputDir="potts-run"
mkdir -p $outputDir
for r in $(seq $beginRun 1 $endRun)
do
	outfile=out_sim_$r\_10samples.txt
	./potts-entropy glauber -L 32 -U 100000 --samples 10 --seed 0 --run-ID $r --init-mode 5 -S 10000 --T-count $tsteps --T-min $tmin --T-max $tmax --threads $numCpus --out-dir $outputDir > $outfile 2>&1
	mv $outfile $outputDir
done
