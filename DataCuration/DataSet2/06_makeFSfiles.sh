#!/bin/bash

cd Outputs/fineStructureOnlyNat

mkdir Inputs/
a=$(ls  shapeITphased/withCM.6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr*_alignedRef.bim | wc -l)
if [ $a != 22 ]
then
	Rscript ../../6a_putGenMap_inMapFile.R \
	T \
	shapeITphased/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr \
	T \
	T \
	/Volumes/MARIOLOLO/Data/genetMap_1KG/genetic_map_chr \
	6 \
	shapeITphased/withCM.6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr
else
	echo cM ok!
fi
for chr in {1..22}
do
	echo $chr
	echo "start.pos recom.rate.perbp" > Inputs/chr$chr.recomb
	awk 'BEGIN{prevCM=0;prevBP=0} {rate=($3-prevCM)/($4-prevBP); print $4,rate; prevCM=$3; prevBP=$4}' shapeITphased/withCM.6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef.bim" >> Inputs/chr$chr.recomb
	if [  ! -s Inputs/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef_phased.phase" ]
        then
                perl ~/Documents/PostDoc/scripts/Tools/fineSTRUCTURE_4.0.0/impute2chromopainter.pl \
			shapeITphased/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef_phased.haps" \
			Inputs/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef_phased"
        else
                echo Inputs/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr$chr"_alignedRef_phased.phase" already generated
        fi

done

awk '{if(NR>2)print $2,$1,1}' shapeITphased/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr22_alignedRef_phased.sample > Inputs/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.ids 

test=$(awk '{print $2}' Inputs/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.ids | uniq | wc -l)
test2=$(awk '{print $2}' Inputs/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.ids | sort | uniq | wc -l)
if [ $test != $test2 ]
then
	echo "your .phase file must be orderd by individuals... need to fix this"
	exit
fi
awk '{print $2}' Inputs/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.ids | uniq -c | awk '{print $2,$1}' > Inputs/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.initialDonorList


