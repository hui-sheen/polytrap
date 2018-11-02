#!/bin/bash

if [ ! -d "log" ]; then
    mkdir log
fi
echo "newtracts.R is creating tracts/breaks for the following genomes:"
gns=("$@") # all genomes supplied in arguments
for gn in "${gns[@]}"; do 
	echo -e "\t$gn"
	Rscript newtracts.R $gn &>>log/newtracts.Rout #rn6 galGal5 danRer10 dm6 sacCer3 &>log/newtracts.Rout
done
mapfile -t gnGffs < ./genomes.meta
#gnGffs=("${gnGffs0[@]:1}")
for gn in "${gns[@]}"; do
	for ((i=0; i<${#gnGffs[@]}; ++i)); do
		gffi=${gnGffs[i]} 
		if [[ "$gffi" == $gn,* ]]; then
			simGff=${gffi/*,/}
			lower_gn=`echo $gn | tr '[:upper:]' '[:lower:]'`
			if [[ "$simGff" != "NA" ]]; then
				echo "for genome $gn mapped to $lower_gn, newtracts.pl annotates all CSVs with GTF $simGff"
				./newtracts.pl $lower_gn ./genomes/$simGff &>log/newtracts.$gn.pOut
				#sleep 5
			fi
			break
		fi
	done
done
echo "calc2stats.R is calculating tract statistics for the following genomes:"
for gn in "${gns[@]}"; do
	echo $gn
	Rscript calc2stats.R $gn
done
gns=("$@")
gzip -k newtracts/*.csv
for gn in "${gns[@]}"; do
	lower_gn=`echo $gn | tr '[:upper:]' '[:lower:]'`
	mv newtracts/*$lower_gn*.csv.gz ../tracts/
	cp newtracts/*$lower_gn*stat* ../tracts/
done
echo "JOB DONE!"
