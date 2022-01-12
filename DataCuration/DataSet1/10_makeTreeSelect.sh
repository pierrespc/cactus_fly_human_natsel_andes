#!/bin/bash


mkdir Outputs/TreeSelect/

grep Titicaca_Uro  ../../Gnecchio/Reformate/Outputs/SouthAmerica_dataset_QCed.fam | awk '{print "Andean",$2,$3,$4,$5,"Andean"}' > Outputs/RemoveUros.RM 

for i in Yanesha Wichi
do
	if [[ $i == "Yanesha" ]]
	then
		grep Wichi ../../Gnecchio/Reformate/Outputs/SouthAmerica_dataset_QCed.fam  > Outputs/TreeSelect/ForYanesha.RM
		cat  Outputs/RemoveUros.RM >> Outputs/TreeSelect/ForYanesha.RM
	else
		grep Yanesha Outputs/TreeSelect/AndeanHanYanesha.SNPFiltered.fam  > Outputs/TreeSelect/ForWichi.RM
		cat  Outputs/RemoveUros.RM >> Outputs/TreeSelect/ForWichi.RM
	fi
	~/src/plink1.9/plink --bfile ../../Gnecchio/Reformate/Outputs/1KG-CHB/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.1KG_CHB --remove Outputs/TreeSelect/For$i.RM --make-bed --out Outputs/TreeSelect/AndeanHan$i
	 ~/src/plink1.9/plink --bfile Outputs/TreeSelect/AndeanHan$i --maf 0.01 --geno 0.02 --mind 0.05 --make-bed --out Outputs/TreeSelect/AndeanHan$i.SNPFiltered

done
	


