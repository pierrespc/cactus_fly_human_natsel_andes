#!/bin/sh

folder=$1
folderBED=$folder
prefBED=$2
rootOut=$folder/Admixture/
mkdir $rootOut
cd $rootOut
###Cross-Validation
if [[ ! -e $prefBED.pruned.bed ]]
then
	cp $folderBED/$prefBED.pruned.bed .
	cp $folderBED/$prefBED.pruned.bim .
	cp $folderBED/$prefBED.pruned.fam .
	echo "copied"
fi
if [[ ! -e $prefBED.pruned.bed ]]
then
	~/src/plink1.9/plink --bfile $folderBED/$prefBED --indep-pairwise 50 5 0.5 --out $prefBED.pruned
	~/src/plink1.9/plink --bfile $folderBED/$prefBED --extract $prefBED.pruned.prune.in --make-bed --out $prefBED.pruned
else 
	echo already in folder...
fi
for run in {1..10}
#for run in 10
do
	mkdir RUN$run
	mv $prefBED.pruned.bed RUN$run/
	mv $prefBED.pruned.fam RUN$run/
	mv $prefBED.pruned.bim RUN$run/
	cd RUN$run
	echo K CVscore > $prefBED.pruned.CV

	for K in {3..10}
	#for K in 8 
	do
		if [[ ! -e $prefBED.pruned.${K}.Q ]]
		then

			/Users/pierrespc/Documents/PostDoc/scripts/Tools/admixture_macosx-1.3.0/admixture --seed $run --cv $prefBED.pruned.bed $K | tee $prefBED.pruned.K${K}.out
			#echo /Users/pierrespc/Documents/PostDoc/scripts/Tools/admixture_macosx-1.3.0/admixture $prefBED.pruned.bed $K
		else 
			echo run $run $prefBED.pruned.K${K}.out already exists!
		fi
		C=$(grep -h CV $prefBED.pruned.K$K.out | awk '{print $4}')
		echo $C
		echo $K $C >> $prefBED.pruned.CV
	done
	mv $prefBED.pruned.bed ../
	mv $prefBED.pruned.fam ../
	mv $prefBED.pruned.bim ../
	cd ..
done


 
