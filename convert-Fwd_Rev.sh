#!/bin/bash

dir=`pwd`
Fwd=$dir/T3_5-Fwd

for chip in `ls $dir | grep Chip`
do

    for barcode in `ls $dir/$chip | grep barcode`
    do

	cd $Fwd/$chip/$barcode


	python3 $dir/Scripts/convert-Fwd_Rev.py --inp inserts.out --ncol 4 --icol 3 --pcol 4 --acol 2 --bcol 1 

	echo $chip $barcode "Done"
	cd $dir

    done
done
