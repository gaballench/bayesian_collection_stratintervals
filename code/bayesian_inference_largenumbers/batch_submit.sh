#!/usr/bin/bash

mkdir batchsims
cd batchsims

SIGMA=(1 20)
NDATA=(10 50 100 200 500)
#NDATA=(5 10)
NSIMS=1000
#NSIMS=5
MCMCITERS=10000
#MCMCITERS=100
SEED=1507
SAMPLING=("theta1" "theta2" "both")

for sampling in ${SAMPLING[@]}; do
    for ndata in ${NDATA[@]}; do
	for sigma in ${SIGMA[@]}; do
	    #echo $sampling
	    #echo $ndata
	    #echo $sigma
	    cp ../template.pbs ${sampling}_${ndata}_${sigma}.pbs
	    sed -i -e s/_SAMPLING/$sampling/g ${sampling}_${ndata}_${sigma}.pbs
	    sed -i -e s/_NDATA/$ndata/g ${sampling}_${ndata}_${sigma}.pbs
	    sed -i -e s/_SIGMA/$sigma/g ${sampling}_${ndata}_${sigma}.pbs
	    sed -i -e s/_SEED/$SEED/g ${sampling}_${ndata}_${sigma}.pbs
	    sed -i -e s/_MCMCITERS/$MCMCITERS/g ${sampling}_${ndata}_${sigma}.pbs
	    sed -i -e s/_NSIMS/$NSIMS/g ${sampling}_${ndata}_${sigma}.pbs
	    qsub ${sampling}_${ndata}_${sigma}.pbs
	    let "SEED++"
	done
    done
done
