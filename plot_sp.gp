set yrange [0:]
set xlabel xname
set ylabel "Speedup"
set terminal png
set output output
plot data u 1:2 title name with linespoints lt rgb "#0099cc";
