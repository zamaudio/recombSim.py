#!/bin/bash
GENOMESIZE=20000
STEPS=50

for maxmut in 0.01 0.1 0.5 1
do
	dir1=sim${maxmut}_${STEPS}_${GENOMESIZE}
	mkdir ${dir1}
	cd ${dir1}
	for rep in 0 1 2 3
	do
		mkdir ${rep}
		cd ${rep}
		recombsim -g ${GENOMESIZE} -m ${maxmut} -s ${STEPS}
		cd ..
	done
	cd ..
done
