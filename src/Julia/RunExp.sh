#!/usr/bin/bash

for mu in {1.05,1.1,1.15,1.2,1.25} 
do	
	for gamma in {0..10}
	do
		let gamma_new=$gamma*0.01
		for delta in {0}
		#{0..10}
			do
				delta_new=$delta*0.01
		        echo  $mu $gamma_new $delta_new
    			#run_julia_single_cpu $seed $N $T &
				julia maxmin_test.jl $mu $gamma_new $delta_new &>/dev/null &
			done
		done	
	done	
done
