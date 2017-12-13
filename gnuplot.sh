#!/bin/bash


infile="$1"
outfile="$2"

echo "
set key font \",20\"
set datafile separator ','
set term png size 1000, 800  font 'times.ttf,13'
set xlabel \"k\" font 'times.ttf,20'
set ylabel \"|Mk|\" font 'times.ttf,20'
set zlabel \"N(Mk)\" font 'times.ttf,20'
set output \"$outfile\"  
splot \"$infile\" using 1:2:3 title \" N(Mk) VS  |Mk| during iteration\" 
" > sample.gp && gnuplot sample.gp
