#!/bin/bash

dir=`pwd`

for chip in `ls $dir/Combined | grep Chip`
do
    
    cd $dir/Combined/$chip
    if [ -f lengths.dat ]
    then
	rm -f lengths.dat
    fi
    touch lengths.dat
    
    for barcode in `ls $dir/Combined/$chip | grep barcode`
    do

	cd $dir/Combined/$chip/$barcode
	python3 $dir/Scripts/lengths_inserts-bulk.py --inp inserts.out --ncol 4 --icol 2 --Nsample 10000 > lengths_inserts.dat
	
	cd $dir/Combined/$chip
	echo $barcode > lengths.temp
	cat $dir/Combined/$chip/$barcode/lengths_inserts.dat >> lengths.temp
	
	paste lengths.dat lengths.temp > lengths.tmp
	mv lengths.tmp lengths.dat
	rm lengths.temp

	echo $chip $barcode "done"

    done
    
    if [ -f lengths.dat ]
    then
	python3 ../../Scripts/lengths_to_violin.py --inp lengths.dat
    fi

done
