#!/bin/bash

dir=`pwd`
Fwd=$dir/T3_5-Fwd
Rev=$dir/T3_5-Rev

if [ ! -d $dir/Combined ]
then
    mkdir $dir/Combined
fi

for chip in `ls $Fwd | grep Chip`
do

    if [ ! -d $dir/Combined/$chip ]
    then
	mkdir $dir/Combined/$chip
    fi
    
    for barcode in `ls $Fwd/$chip | grep barcode`
    do

	if [ ! -d $dir/Combined/$chip/$barcode ]
	then
	    mkdir $dir/Combined/$chip/$barcode
	fi
	
	cd $dir/Combined/$chip/$barcode

	
	cat $Fwd/$chip/$barcode/inserts_Rev.out $Rev/$chip/$barcode/inserts.out > $dir/Combined/$chip/$barcode/inserts.out
	python3 $dir/Scripts/score_inserts-bulk.py --inp inserts.out --ncol 4 --icol 2
	
	numis=`python3 $dir/Scripts/inserts_to_per_template.py --inp inserts.out --ncol 4 --icol 2 --pcol 1 --acol 3 --bcol 4 --spol Rev --ncut 10`
	
	echo $chip $barcode $numis

	if [ -f inserts_per_template.out ]
	then

	    python3 $dir/Scripts/score_inserts.py --inp inserts_per_template.out
	    
	fi
	
	cd $dir

    done
done
