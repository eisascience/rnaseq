#!/bin/bash

set -e
set -x
BASEDIR=/home/exacloud/lustre1/ONPRCPGP/bimber/rnaseq-singlecell

function doWork(){
	#directory for analysis
	if [ ! -e $4 ];then
		mkdir $4
	fi

	cd $4 

	#directory for read source
	if [ ! -e $2 ];then
		mkdir $2
	fi

	cd $2

	#cp $BASEDIR/ImmuneGeneList.txt ./
	touch de.R
	truncate -s 0 de.R
	
	cp $BASEDIR/${3} ./
	cp ../../Metadata.txt ./
	if [ -e ../../filter.R ];then
		cp ../../filter.R ./	
		cat filter.R > de.R
		echo "" >> de.R
	fi

	cat $BASEDIR/de.R >> de.R
	
	WD=`pwd`
	docker run --rm=true -v "${WD}:/work" -t bbimber/rnaseq Rscript /work/de.R ${3} ${1}

	count=`ls -1 *.png 2>/dev/null | wc -l`
	if [ $count != 0 ]
	then 
		rm *.png
	fi 
	#docker run --rm=true -v "${WD}:/work" -t bbimber/rnaseq Rscript /work/makePlots.R ${1}

	#rm ImmuneGeneList.txt
	
	TARGET=/home/groups/prime-seq/production/Internal/Bimber/142/\@files/rnaseq/${5}/${4}/${2}
	if [ ! -e $TARGET ];then
		mkdir -p $TARGET
	fi

	rm -Rf ${TARGET}/*
	cp -R ./ $TARGET

	cd ../../
}

#ensure up to date
docker pull bbimber/rnaseq

WD=`pwd`
WORKDIR=`basename $WD`

#First using STAR's counts:
doWork $1 'starCounts' 'TCR_1263.txt' $1 $WORKDIR
#doWork $1 'subreadCounts' 'geneCounts.txt' $1 $WORKDIR
