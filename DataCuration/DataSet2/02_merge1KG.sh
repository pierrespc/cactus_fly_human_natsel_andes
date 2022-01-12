#!/bin/sh


mkdir Outputs/1KG
cd Outputs/
rm 1KG/1KG.keep
grep AFR /Volumes/MARIOLOLO/PolymorphismDataFromLitterature/1KG_Phase3/integrated_call_samples_v3.20130502.ALL.panel  | awk '{print $1"\t"$1}' > 1KG/1KG.keep
grep EUR /Volumes/MARIOLOLO/PolymorphismDataFromLitterature/1KG_Phase3/integrated_call_samples_v3.20130502.ALL.panel  | awk '{print $1"\t"$1}' >> 1KG/1KG.keep
grep AMR /Volumes/MARIOLOLO/PolymorphismDataFromLitterature/1KG_Phase3/integrated_call_samples_v3.20130502.ALL.panel  | awk '{print $1"\t"$1}' >> 1KG/1KG.keep

awk '{print $2}' SouthAmerica_dataset_QCed.KinshipFilter.bim > 1KG/SNPs.extract

rm 1KG/mergeList
for chr in {1..22}
do
	if [[ ! -e  1KG/1KG.chr$chr.bed ]]
	then
		~/src/plink1.9/plink --vcf /Volumes/MARIOLOLO/PolymorphismDataFromLitterature/1KG_Phase3/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep 1KG/1KG.keep --extract 1KG/SNPs.extract --make-bed --out  1KG/1KG.chr$chr
	fi
	if [[ $chr != 1 ]]
	then
		echo 1KG/1KG.chr$chr.{bed,bim,fam} >> 1KG/mergeList
	fi
done

First=T
if [[ ! -e 1KG/1KG.GW.bed ]]
then
	~/src/plink1.9/plink --bfile 1KG/1KG.chr1 --merge-list 1KG/mergeList --make-bed --out 1KG/1KG.GW
else
	echo  1KG/1KG.GW.bed already generated
fi

if [[ ! -e 1KG/1KG.GW.bed ]]
then
	~/src/plink1.9/plink --bfile 1KG/1KG.chr1 --merge-list 1KG/mergeList --exclude 1KG/1KG.GW-missnp --make-bed --out 1KG/1KG.GW

fi
if [[ ! -e 1KG/1KG.GW.bed ]]
then
	echo 1KG/1KG.GW.bed not generated 
	exit
fi

if [[ ! -e 1KG/1KG.GW.Flipped.bed ]]
then
	Rscript 02a_ArrangeStrandWithRefBim.R 1KG/1KG.GW SouthAmerica_dataset_QCed.KinshipFilter 1KG/1KG.GW.Flipped ~/src/plink1.9/plink
fi
 
awk '{print $2}' 1KG/1KG.GW.Flipped.bim > 1KG/SNPsCommon.bim

~/src/plink1.9/plink --bfile SouthAmerica_dataset_QCed.KinshipFilter --extract 1KG/SNPsCommon.bim --bmerge 1KG/1KG.GW.Flipped.{bed,bim,fam} --make-bed --allow-no-sex --out 1KG/SouthAmerica_dataset_QCed.KinshipFilter.1KG 

