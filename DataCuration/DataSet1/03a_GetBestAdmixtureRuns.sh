#!/bin/bash

folder=$1
folderBED=$folder
prefBED=$2
#prefED=FinalSetPierre.1KG.GENO.MIND.MAF
rootOut=$folder/Admixture/
mkdir $rootOut/BestRUNperK/
cd $rootOut

echo K RUN log CVscore > BestRUNperK/ListBESTruns
for k in {3..10}
do
	rm BestRUNperK/K$k.LogLik
	for run in {1..10}
	do 
		if [ ! -e  RUN$run/$prefBED.pruned.$k.P ]
		then
			continue
		fi
		#log=$(grep Loglik RUN$run/$prefBED.pruned.K$k.out | grep -v delta | awk '{print $2}')
		log=$(grep CV RUN$run/$prefBED.pruned.K$k.out | cut -d' ' -f4)
		echo $run $log >> BestRUNperK/K$k.LogLik
	done
	#RUNbest=$(sort -k2 -n BestRUNperK/K$k.LogLik | tail -1 | awk '{print $1}')
	RUNbest=$(sort -k2 -n BestRUNperK/K$k.LogLik | head -1 | awk '{print $1}')
	#Logbest=$(sort -k2 -n BestRUNperK/K$k.LogLik | tail -1 | awk '{print $2}')
	Logbest=$(sort -k2 -n BestRUNperK/K$k.LogLik | head -1 | awk '{print $2}')
	CVbest=$(grep "CV error" RUN$RUNbest/$prefBED.pruned.K$k.out | awk '{print $4}')
	echo $k $RUNbest $Logbest $CVbest >> BestRUNperK/ListBESTruns
	cp RUN$RUNbest/$prefBED.pruned.$k.P BestRUNperK/
	cp RUN$RUNbest/$prefBED.pruned.$k.Q BestRUNperK/
	cp RUN$RUNbest/$prefBED.pruned.K$k.out BestRUNperK/
done
	
	
		

