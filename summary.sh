#!/bin/bash

dir=`pwd`
Fwd=$dir/T3_5-Fwd
Rev=$dir/T3_5-Rev

for chip in `ls $Fwd | grep Chip`
do

    for barcode in `ls $Fwd/$chip | grep barcode`
    do



	nrev=`grep "T31p," $Rev/$chip/$barcode/sequential_summary.csv | sed "s|\,|\ |g" | awk '{print $2}'`
	nfwd=`grep "T35," $Fwd/$chip/$barcode/sequential_summary.csv |  sed "s|\,|\ |g" | awk '{print $2}'`

	prev=`grep "T31p," $Rev/$chip/$barcode/sequential_summary.csv | sed "s|\,|\ |g" | awk '{print $3}'`
	pfwd=`grep "T35," $Fwd/$chip/$barcode/sequential_summary.csv |  sed "s|\,|\ |g" | awk '{print $3}'`


	ntot=`echo $nfwd $nrev | awk '{print $1+$2}'`
	ptot=`echo $pfwd $prev | awk '{print $1+$2}'`
	
	echo $chip $barcode $ntot $ptot
    done
done
