#!/usr/bin/bash
MAXPROCS=9

for mu in 1.1
do
	for gamma in {0,1,2,3,4,5}
	do
		gamma_step=0.01
		gamma_new=$(echo $gamma*$gamma_step | bc)
		for delta in {0,1,2,3,4,5,6,7,8,9,10}
		do
			delta_step=0.01
			delta_new=$(echo $delta*$delta_step | bc)
		    echo  $mu $gamma_new $delta_new
    		#run_julia_single_cpu $seed $N $T &
			echo output_${mu}_${gamma_new}_${delta_new}.txt
			julia maxmin_test.jl $mu $gamma_new $delta_new &> ./Logs/output_${mu}_${gamma_new}_${delta_new}.txt &
			sleep 1
			if [ $(jobs -r | wc -l) -ge $MAXPROCS ]; then
      	wait $(jobs -r -p | head -1)
			fi
		done
	done
done
