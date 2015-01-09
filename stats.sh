#!/bin/bash

m=256
if [ "$1" != "" ]; then
	m="$1"
fi

i=16
while [ "$i" -lt "$m" ]; do
	sum_time1=0
	sum_time2=0
	np=4
	echo -e "\n sequential"
	for j in 1 2 3 4
	do
		ret=`make exec n=$np m=$i seq=1`
		sum_time1=$(echo "${ret} + ${sum_time1}" | bc)
	done
	echo -e "\n parallel"
	for j in 1 2 3 4
	do
		ret=`make exec n=${np} m=${i} seq=0`
		echo "${ret}"
		sum_time2=$(echo "${ret} + ${sum_time2}" | bc)
	done
	mean_sp=$(echo "scale=6;${sum_time1}/${sum_time2}" | bc)
	echo "${i} ${mean_sp}"
	i=$(($i+16))
done
