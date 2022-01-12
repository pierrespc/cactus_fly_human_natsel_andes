#!/bin/sh


root=~/Documents/PostDoc/DetoxFly/LAST/Gnecchio_noUros/

for stat in mean median
do
	for g in 1
	do
		for fl in 0 5 10
		do
			echo KB $fl $g
			Rscript $root/AnalysesScripts/makeSummaryScorePerGene.R \
				$root/IHS/Andean.absIHS \
				$root/../../OrthoGenes.txt \
				$root/../../GeneStartEnd.BED \
	        		2 \
				3 \
				4 \
				$stat \
				$root/IHS/Andean.DetoxGenes.$fl"KB.absIHS.$stat."$g"Groups.tsv" \
	        		$root/IHS/Andean.BKGDGenes.$fl"KB.absIHS.$stat."$g"Groups.tsv" \
				$g \
				$fl
		done
	done
done

