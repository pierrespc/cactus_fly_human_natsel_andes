#!/bin/bash


for i in bed bim fam
do
	cp ../MaoMoreno/Reformate/Outputs/6pops.KinshipFilter.AdmFilter.CPFS.$i  MaoMoreno/MaoMoreno.$i
	cp ../Gnecchio/Reformate/Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.$i Gnecchio/Gnecchio.$i
done

Rscript putOriLab.R  MaoMoreno/MaoMoreno.fam ../MaoMoreno/Reformate/Outputs/6pops.KinshipFilter.fam
Rscript putOriLab.R Gnecchio/Gnecchio.fam ../Gnecchio/Reformate/Outputs/SouthAmerica_dataset_QCed.KinshipFilter.fam

mkdir Merged
echo "geno1: MaoMoreno/MaoMoreno.bed
snp1:  MaoMoreno/MaoMoreno.bim
ind1:  MaoMoreno/MaoMoreno.fam 
geno2: Gnecchio/Gnecchio.bed
snp2:  Gnecchio/Gnecchio.bim
ind2:  Gnecchio/Gnecchio.fam
genooutfilename:   Merged/Merged.bed
snpoutfilename:    Merged/Merged.bim
indoutfilename:    Merged/Merged.fam
outputformat: PACKEDPED " > Merged/Merging.param

 ~/Documents/PostDoc/scripts/Tools/EIG7.2.0/bin/mergeit -p Merged/Merging.param

awk '{print $2}' Merged/Merged.fam | awk -F ":" '{print $1,$2,0,0,0,$1}' > tmp
mv tmp  Merged/Merged.fam 


mkdir Merged_noUros
grep Uros  Merged/Merged.fam  > Merged_noUros/Merged_noUros.RM
~/src/plink1.9/plink --bfile Merged/Merged --remove Merged_noUros/Merged_noUros.RM --make-bed --out Merged_noUros/Merged_noUros

for set in Merged_noUros Merged Gnecchio MaoMoreno
do

	~/src/plink1.9/plink --bfile  $set/$set --maf 0.01 --geno 0.02 --mind 0.05 --make-bed --out $set/$set.MAF0.01.GENO0.02.MIND0.05
	awk '{$6=$1;print $0}' $set/$set.MAF0.01.GENO0.02.MIND0.05.fam > tmp
	mv tmp  $set/$set.MAF0.01.GENO0.02.MIND0.05.fam
	mkdir $set/FSTpost/
	root=$(pwd)/$set/
	pref=$set.MAF0.01.GENO0.02.MIND0.05
	###make FST and PCA
echo "genotypename:    $root/$pref.bed
snpname:         $root/$pref.bim
indivname:         $root/$pref.fam
evecoutname:     $root/FSTpost/$pref.evec
evaloutname:     $root/FSTpost/$pref.eval
numoutevec:    16
numoutlieriter:    0 
fstonly:        yes" > $root/FSTpost/FST.param

	if [ ! -s $root/FSTpost/$pref.evec ]
	then
		~/Documents/PostDoc/scripts/Tools/EIG7.2.0/bin/smartpca -p $root/FSTpost/FST.param > $root/FSTpost/FST.log
	else
		echo $root/FSTpost/$pref.evec already exists 
	fi
	Rscript plotPCA_FINAL.R \
		$root/FSTpost/$pref.evec \
		$root/FSTpost/$pref.eval \
		ColorManifest.tsv \
		$root/FSTpost/$pref.pdf \
		16
	#retrieve FST matrix
	grep "population:" $root/FSTpost/FST.log | awk '{print $2,$3}' > $root/FSTpost/FST.PopIndexes

	awk 'BEGIN {pri=0} {if($0=="fst *1000:"){pri=1;next}; if(pri==1){if($0!=""){o=$1; for(i=2;i<=NF;i=i+1){o=o" "$i};print o}else{pri=0}}}'  $root/FSTpost/FST.log > $root/FSTpost/FST.Matrix
awk 'BEGIN {pri=0} {if($0=="s.dev * 1000000:"){pri=1;next}; if(pri==1){if($0!=""){o=$1; for(i=2;i<=NF;i=i+1){o=o" "$i};print o}else{pri=0}}}'  $root/FSTpost/FST.log > $root/FSTpost/FST.Std.Matrix


	Rscript plotFST.R $root/FSTpost/
	Rscript MakeMDS.R $root/FSTpost/ $pref
	exit

done


