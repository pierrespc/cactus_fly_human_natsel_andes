#!/bin/bash


mkdir Outputs/TreeSelect/
for i in Yukpa Bari
do
	grep $i Outputs/1KG-CHB/6pops.KinshipFilter.AdmFilter.CPFS.1KG_CHB.fam > Outputs/TreeSelect/AndeanHan$i.KEEP
	for j in Han Andean 
	do 
		grep $j Outputs/1KG-CHB/6pops.KinshipFilter.AdmFilter.CPFS.1KG_CHB.fam >> Outputs/TreeSelect/AndeanHan$i.KEEP
	done 
	~/src/plink1.9/plink --bfile Outputs/1KG-CHB/6pops.KinshipFilter.AdmFilter.CPFS.1KG_CHB --keep Outputs/TreeSelect/AndeanHan$i.KEEP --make-bed --out Outputs/TreeSelect/AndeanHan$i
	 ~/src/plink1.9/plink --bfile Outputs/TreeSelect/AndeanHan$i --maf 0.01 --geno 0.02 --mind 0.05 --make-bed --out Outputs/TreeSelect/AndeanHan$i.SNPFiltered

done
	


