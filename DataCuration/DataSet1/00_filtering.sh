#!/bin/bash


awk '{if($1>=1 && $1<=22)print $2}' /Volumes/MARIOLOLO/PolymorphismDataFromLitterature/GnnuchioRuscone_SouthAmerica/SouthAmerica_dataset_pre_QC.bim > Outputs/Autosomal.SNPs
~/src/plink1.9/plink --bfile /Volumes/MARIOLOLO/PolymorphismDataFromLitterature/GnnuchioRuscone_SouthAmerica/SouthAmerica_dataset_pre_QC --extract Outputs/Autosomal.SNPs --maf 0.01 --geno 0.02 --mind 0.05 --make-bed --out Outputs/SouthAmerica_dataset_QCed
