#!/bin/bash

m=256
if [ "$1" != "" ]; then
	m="$1"
fi

i=128
nb_test=10
while [ "$i" -lt "$m" ]; do
	sum_time1=0
	sum_time2=0
	np=4
	echo -e "\n sequential"
	for j in `seq 1 ${nb_test}`
	do
		ret=`make -s exec n=${np} m=${i} seq=1`
		sum_time1=$(echo "${ret} + ${sum_time1}" | bc)
		echo "${ret}"
	done
	mean_time1=$(echo "scale=6;${sum_time1}/${nb_test}" | bc)
	echo "$mean_time1"
	echo -e "\n parallel"
	for j in `seq 1 ${nb_test}`
	do
		ret=`make -s exec n=${np} m=${i} seq=0`
		sum_time2=$(echo "${ret} + ${sum_time2}" | bc)
		echo "${ret}"
	done
	mean_time2=$(echo "scale=6;${sum_time2}/${nb_test}" | bc)
	echo "$mean_time2"
	mean_sp=$(echo "scale=6;${sum_time1}/${sum_time2}" | bc)
	i=$(($i+16))
done
