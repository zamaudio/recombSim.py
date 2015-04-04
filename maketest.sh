#!/bin/bash
INFILE=$1

for n in 1 5 10; do
	for d in 0.0005 0.001 0.002 0.005 0.01 0.02; do
		for l in 100 1000 5000 10000; do
			for s in 1 2 3 4 5 6 7 8 9 10; do
				python recombTest.py -i ${INFILE} -d${d} -l${l} -n${n} -o${INFILE}_d${d}_l${l}_n${n}_s${s} -s${s}
			done
		done
	done
done
