#!/bin/bash

if [ "$1" != "" ]; then
	m="$1"
fi

i=512
m=2048
nb_test=5
np=2
n=16
# while [ "$i" -le "$m" ]; do
	sum_time1=0
	# sequential
	for j in `seq 1 ${nb_test}`
	do
		ret=`make -s exec n=${np} m=${i} seq=1`
		sum_time1=$(echo "${ret} + ${sum_time1}" | bc)
	done
	mean_time1=$(echo "scale=6;${sum_time1}/${nb_test}" | bc)
	while [ "$np" -le "$n" ]; do
		sum_time2=0
		# parallel $np procs
		for j in `seq 1 ${nb_test}`
		do
			ret=`make -s exec n=${np} m=${i} seq=0`
			sum_time2=$(echo "${ret} + ${sum_time2}" | bc)
		done
		mean_time2=$(echo "scale=6;${sum_time2}/${nb_test}" | bc)
		# echo "$mean_time2"
		mean_sp=$(echo "scale=6;${sum_time1}/${sum_time2}" | bc)
		echo "$np $mean_sp"
	# i=$(($i+256))
	np=$(($np+2))
done
# done
