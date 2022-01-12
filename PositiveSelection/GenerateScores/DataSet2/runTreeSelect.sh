#!/bin/bash


outFold=/Users/pierrespc/Documents/PostDoc/DetoxFly/LAST/MaoMoreno/Treeselect/
inFold=/Users/pierrespc/Documents/PostDoc/DetoxFly/LAST/MaoMoreno/Reformate/Outputs/TreeSelect/

for i in Bari Yukpa
do
	mkdir $outFold/ref$i/
	cd $outFold/ref$i/
	mkdir Frq/
	for pop in $i Han Andean
	do
		grep -w $pop $inFold/AndeanHan$i.fam > Frq/$pop.KEEP
		~/src/plink1.9/plink --bfile $inFold/AndeanHan$i.SNPfiltered --keep Frq/$pop.KEEP --freq --out Frq/$pop
		
	done
	pwd
	ls
	perl ~/Documents/These/scripts/Tools/TreeSelect/RunLRT/Make_Input_Run_TreeSelect.pl \
		$(pwd)/Frq/ \
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
		

		

	
	
	
