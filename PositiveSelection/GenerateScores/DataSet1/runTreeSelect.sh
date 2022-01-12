#!/bin/bash


outFold=/Users/pierrespc/Documents/PostDoc/DetoxFly/LAST/Gnecchio_noUros/Treeselect/
inFold=/Users/pierrespc/Documents/PostDoc/DetoxFly/LAST/Gnecchio_noUros/Reformate/Outputs/TreeSelect/
for i in Yanesha Wichi
do
	mkdir $outFold/ref$i/
	cd $outFold/ref$i/
	for pop in $i Han Andean
	do
		grep -w $pop $inFold/AndeanHan$i.fam > $pop.KEEP
		~/src/plink1.9/plink --bfile $inFold/AndeanHan$i.SNPfiltered --keep $pop.KEEP --freq --out $pop
	done
	perl ~/Documents/These/scripts/Tools/TreeSelect/RunLRT/Make_Input_Run_TreeSelect.pl \
		$(pwd) \
		F \
		$inFold/AndeanHan$i.SNPfiltered.bim \
		0.01 \
		$(pwd)/Results/ \
		$(pwd)/Temp/ \
		~/Documents/These/scripts/Tools/TreeSelect/ \
		3 \
		Andean Han $i
	exit
done
		

		

	
	
	
