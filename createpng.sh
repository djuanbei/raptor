#!/bin/bash
allfiles=`ls  *saturate*.csv`


for ff in  $allfiles
do

    bbname=${ff%.csv}
    echo $ff  $bbname
    ./gnuplot.sh $ff $bbname.png
done

