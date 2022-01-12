#!/bin/sh


mkdir Outputs/1KG-CHB/
cd Outputs/1KG-CHB/

grep CHB /Volumes/MARIOLOLO/PolymorphismDataFromLitterature/1KG_Phase3/integrated_call_samples_v3.20130502.ALL.panel  | awk '{print $1"\t"$1}' > 1KG-CHB.keep

awk '{print $2}' ../SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.bim > 1KG-CHB.extract

rm 1KG-CHB.mergeList
for chr in {22..1}
do
	if [[ ! -e  1KG-CHB.chr$chr.bed ]]
	then
		~/src/plink1.9/plink --vcf /Volumes/MARIOLOLO/PolymorphismDataFromLitterature/1KG_Phase3/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep 1KG-CHB.keep --extract 1KG-CHB.extract --make-bed --out  1KG-CHB.chr$chr
	fi
	if [[ $chr != 1 ]]
	then
		echo 1KG-CHB.chr$chr.{bed,bim,fam} >> 1KG-CHB.mergeList
	fi
done

if [[ ! -e 1KG-CHB.GW.bed ]]
then
	~/src/plink1.9/plink --bfile 1KG-CHB.chr1 --merge-list 1KG-CHB.mergeList --make-bed --out 1KG-CHB.GW
else
	echo  1KG-CHB.bed already generated
fi

if [[ ! -e 1KG-CHB.GW.bed ]]
then
	for chr in {1..22}
	do
		~/src/plink1.9/plink --bfile 1KG-CHB.chr$chr  --exclude 1KG-CHB.GW-merge.missnp --make-bed --out 1KG-CHB.chr$chr
	done
        ~/src/plink1.9/plink --bfile 1KG-CHB.chr1 --merge-list 1KG-CHB.mergeList --exclude 1KG-CHB.GW-merge.missnp --make-bed --out 1KG-CHB.GW

fi

if [[ ! -e 1KG-CHB.GW.bed ]]
then
	echo 1KG-CHB.GW.bed not generated 
	exit
fi

if [[ ! -e 1KG-CHB.GW.Flipped.bed ]]
then
	Rscript ../../02a_ArrangeStrandWithRefBim.R 1KG-CHB.GW ../SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS 1KG-CHB.GW.Flipped ~/src/plink1.9/plink
fi
 
awk '{print $2}' 1KG-CHB.GW.Flipped.bim > 1KG-CHB.SNPsCommon.bim

if [[ ! -e SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.1KG_CHB ]]
then
	~/src/plink1.9/plink --bfile ../SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS --extract 1KG-CHB.SNPsCommon.bim --bmerge 1KG-CHB.GW.Flipped.{bed,bim,fam} --make-bed --allow-no-sex --out SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.1KG_CHB 
fi

awk '{if($1!="Wichi" && $1 !="Yanesha" && $1 != "Andean"){$1="Han"}; $6=$1; print $0}' SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.1KG_CHB.fam > SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.1KG_CHB.fam2
mv SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.1KG_CHB.fam2 SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.1KG_CHB.fam
