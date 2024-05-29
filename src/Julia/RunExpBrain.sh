#!/bin/bash

MAXPROCS=1

for HOMCONSTR in {2000,8000,100000000} 	
#100000000
do
	for HOMOGENCONSTR_PERVOX in 1000000 #{2,5,10,1000000}
	do
		for OARCONSTRINIT in {2000,10000000}
		do
		    	echo  $HOMCONSTR $HOMOGENCONSTR_PERVOX $OARCONSTRINIT
			echo output_${HOMCONSTR}_${HOMOGENCONSTR_PERVOX}_${OARCONSTRINIT}.txt
			julia maxmin_test_runtime.jl $HOMCONSTR $HOMOGENCONSTR_PERVOX $OARCONSTRINIT &> ./output_${HOMCONSTR}_${HOMOGENCONSTR_PERVOX}_${OARCONSTRINIT}.txt &
			sleep 1
			if [ $(jobs -r | wc -l) -ge $MAXPROCS ]; then
				wait $(jobs -r -p | head -1)
			fi
		done
	done	
done	
