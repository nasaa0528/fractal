#!/bin/bash 

binary_file="./Fractal"
noimage="--noimage"
for np in {1..4}
do
	for dim in $(seq 256 256 5120)
	do
		for iter in $(seq 100 100 400)
		do 
			$(echo "nproc: ${np} dim: ${dim} iter: ${iter}" >> result.txt)
			$(mpirun -np ${np} ${binary_file} $dim $dim $iter >> result.txt)  
		done
	done
	#echo $np
done
#$(mpirun -np ${np} ./Fractal ${dim} ${dim} ${iter} ${noimage} >> result.txt)

