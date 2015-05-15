#!/bin/bash
INFILE=$1

for n in 1 3 5; do
	for d in 0.001 0.002 0.005 0.01 0.02; do
		for l in 1000 5000 10000; do
			for s in 10 20 30 40 50; do
				python recombSim.py -i ${INFILE} -d${d} -l${l} -n${n} -o${INFILE}_d${d}_l${l}_n${n}_s${s} -s${s}
			done
		done
	done
done
