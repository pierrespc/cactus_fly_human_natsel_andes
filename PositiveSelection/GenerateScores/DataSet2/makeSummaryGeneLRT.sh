#!/bin/sh


root=~/Documents/PostDoc/DetoxFly/LAST/MaoMoreno/

for stat in mean median
do
	for g in 1
	do
		for pop in Bari Yukpa
		do
			for fl in 0 5 10
			do
				echo $pop $fl KB $g 
				if [ -e $root/TreeSelect/ref$pop/Results/Andean.BKGDGenes.$fl"KB.$stat."$g"Groups.tsv" ]
				then
					echo already performed
					continue
				fi
				Rscript $root/AnalysesScripts/makeSummaryScorePerGene.R \
					$root/TreeSelect/ref$pop/Results/Andean.output \
					$root/../../OrthoGenes.txt \
					$root/../../GeneStartEnd.BED \
	        			2 \
					3 \
					8 \
					$stat \
					$root/TreeSelect/ref$pop/Results/Andean.DetoxGenes.$fl"KB.$stat."$g"Groups.tsv" \
		        		$root/TreeSelect/ref$pop/Results/Andean.BKGDGenes.$fl"KB.$stat."$g"Groups.tsv" \
					$g \
					$fl
			done
		done
	done
done

