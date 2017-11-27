#!/bin/bash


infile="$1"
outfile="$2"

echo "
set datafile separator ','
set term png size 800, 600
set xlabel \"The dimension of matrix M during iteration\"
set ylabel \"The number of nonzero elements in matrix M during iteration\"
set output \"$outfile\"  
plot \"$infile\" using 1:2 title \"The number of nonzero elements in matrix M VS The dimension of matrix M during iteration\"
" > sample.gp && gnuplot sample.gp
