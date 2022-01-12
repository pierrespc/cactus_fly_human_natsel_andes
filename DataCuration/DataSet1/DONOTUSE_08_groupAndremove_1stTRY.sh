#!/bin/bash

###this is done looking at the tree (generrated with itol online for example)
##we remove individuals not clustering correctly and rename some pops 
rm Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.Remove
rm Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
grep -w YASH41 Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.fam >> Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.Remove
grep -w SHI6 Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.fam >> Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.Remove
grep -w HUA01 Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.fam >> Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.Remove
grep -w TQ19 Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.fam >> Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.Remove
grep -w TQ29 Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.fam >> Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.Remove

~/src/plink1.9/plink --bfile Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter \
	--remove Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.Remove \
	--make-bed --out Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS


grep Yanesha_IntermediateSelva Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.fam | cut -f1,2 | awk '{print $2,$1,"Yanesha"}' >> Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
grep Yanesha_HighSelva Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.fam | cut -f1,2 | awk '{print $2,$1,"Yanesha"}' >> Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS

echo "B28 Bolivia_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
echo "B12 Bolivia_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
echo "B19 Bolivia_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
echo "B29 Bolivia_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
echo "B24 Bolivia_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
echo "TA2 Titicaca_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
echo "TA14 Titicaca_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
echo "TA15 Titicaca_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS
echo "TA16 Titicaca_Aymara Aymara2" >>Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS

while read line 
do
	a=($line)
	awk -v ind=${a[0]} -v old=${a[1]} -v new=${a[2]} '{if($2==ind){$1=new};$6=$1;print $0}' Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.fam > tmp
	mv tmp Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.fam
done < Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.NEWLABELS

cut -f1 -d' ' Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.fam  | sort | uniq -c

~/src/plink1.9/plink --bfile  Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS --maf 0.01 --geno 0.02 --mind 0.05 --make-bed --out Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.MAF0.01.GENO0.02.MIND0.05

awk '{$6=$1;print $0}' Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.MAF0.01.GENO0.02.MIND0.05.fam > tmp
mv tmp  Outputs/SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.MAF0.01.GENO0.02.MIND0.05.fam

mkdir Outputs/fineStructureOnlyNat/FSTpost/

root=$(pwd)/Outputs/
pref=SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.CPFS.MAF0.01.GENO0.02.MIND0.05


###make FST and PCA
echo "genotypename:    $root/$pref.bed
snpname:         $root/$pref.bim
indivname:         $root/$pref.fam
evecoutname:     $root/fineStructureOnlyNat/FSTpost/$pref.evec
evaloutname:     $root/fineStructureOnlyNat/FSTpost/$pref.eval
numoutevec:    16
numoutlieriter:    0 
fstonly:        yes" > $root/fineStructureOnlyNat/FSTpost/FST.param

if [ ! -s $root/fineStructureOnlyNat/FSTpost/$pref.evec ]
then
	~/Documents/PostDoc/scripts/Tools/EIG7.2.0/bin/smartpca -p $root/fineStructureOnlyNat/FSTpost/FST.param > $root/fineStructureOnlyNat/FSTpost/FST.log
else
	echo $root/fineStructureOnlyNat/FSTpost/$pref.evec already exists 
fi
Rscript Outputs/1KG/PCA/plotPCA_FINAL.R \
	$root/fineStructureOnlyNat/FSTpost/$pref.evec \
	$root/fineStructureOnlyNat/FSTpost/$pref.eval \
	ColorManifest.tsv \
	$root/fineStructureOnlyNat/FSTpost/$pref.pdf \
	16

#retrieve FST matrix
grep "population:" $root/fineStructureOnlyNat/FSTpost/FST.log | awk '{print $2,$3}' > $root/fineStructureOnlyNat/FSTpost/FST.PopIndexes

awk 'BEGIN {pri=0} {if($0=="fst *1000:"){pri=1;next}; if(pri==1){if($0!=""){o=$1; for(i=2;i<=NF;i=i+1){o=o" "$i};print o}else{pri=0}}}'  $root/fineStructureOnlyNat/FSTpost/FST.log > $root/fineStructureOnlyNat/FSTpost/FST.Matrix
awk 'BEGIN {pri=0} {if($0=="s.dev * 1000000:"){pri=1;next}; if(pri==1){if($0!=""){o=$1; for(i=2;i<=NF;i=i+1){o=o" "$i};print o}else{pri=0}}}'  $root/fineStructureOnlyNat/FSTpost/FST.log > $root/fineStructureOnlyNat/FSTpost/FST.Std.Matrix


Rscript 8a_plotFST.R




