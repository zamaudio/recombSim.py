#!/bin/bash
GENOMESIZE=20000
STEPS=50

for maxmut in 1 10
do
	dir1=sim${maxmut}_${STEPS}_${GENOMESIZE}
	mkdir ${dir1}
	cd ${dir1}
	for rep in 0 1 2 3 4 5 6 7 8 9
	do
		mkdir ${rep}
		cd ${rep}
		recombsim -g ${GENOMESIZE} -m ${maxmut} -s ${STEPS}
		cd ..
	done
	cd ..
done
