#!/bin/bash

set -e
set -x

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

	cp /home/exacloud/lustre1/ONPRCPGP/bimber/rnaseq-singlecell/makePlots2.R ./
	cp /home/exacloud/lustre1/ONPRCPGP/bimber/rnaseq-singlecell/RNASeq.rmd ./
	cp /home/exacloud/lustre1/ONPRCPGP/bimber/rnaseq-singlecell/ImmuneGeneTargets.txt ./
	cp /home/exacloud/lustre1/ONPRCPGP/bimber/rnaseq-singlecell/${3} ./
	cp ../../Metadata.txt ./
	if [ -e ../../filter.R ];then
		cp ../../filter.R ./	
	fi

	WD=`pwd`
	docker run --rm=true -v "${WD}:/work" -t bbimber/rnaseq Rscript -e "geneCountTableFile=\"${3}\"; groupColName = \"${1}\"; rmarkdown::render('/work/RNASeq.rmd', clean=TRUE)"

	count=`ls -1 *.png 2>/dev/null | wc -l`
	if [ $count != 0 ]
	then 
		rm *.png
	fi 
	docker run --rm=true -v "${WD}:/work" -t bbimber/rnaseq Rscript /work/makePlots2.R ${1}

	#rm Metadata.txt 
	rm makePlots2.R 
	rm RNASeq.rmd
	rm ImmuneGeneTargets.txt
	
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
BASEDIR=`basename $WD`

#First using STAR's counts:
doWork $1 'starCounts' 'TCR_1263.txt' $1 $BASEDIR
#doWork $1 'subreadCounts' 'geneCounts.txt' $1 $BASEDIR
