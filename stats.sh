#!/bin/bash

m=2048
if [ "$1" != "" ]; then
	m="$1"
fi

i=1024
nb_test=3
np=2
n=4
# while [ "$i" -le "$m" ]; do
	sum_time1=0
	echo -e "\n$i sequential"
	for j in `seq 1 ${nb_test}`
	do
		ret=`make -s exec n=${np} m=${i} seq=1`
		sum_time1=$(echo "${ret} + ${sum_time1}" | bc)
		echo "${ret}"
	done
	mean_time1=$(echo "scale=6;${sum_time1}/${nb_test}" | bc)
	echo "$mean_time1"
	while [ "$np" -le "$n" ]; do
		sum_time2=0
		echo -e "\n$i parallel $np procs"
		for j in `seq 1 ${nb_test}`
		do
			ret=`make -s exec n=${np} m=${i} seq=0`
			sum_time2=$(echo "${ret} + ${sum_time2}" | bc)
			echo "${ret}"
		done
		mean_time2=$(echo "scale=6;${sum_time2}/${nb_test}" | bc)
		echo "$mean_time2"
		mean_sp=$(echo "scale=6;${sum_time1}/${sum_time2}" | bc)
	# i=$(($i+128))
	np=$(($np+1))
done
# done
