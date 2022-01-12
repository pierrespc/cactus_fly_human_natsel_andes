#!/bin/bash

perl 01a_getListIndToRemove_from_RelatedPairs.pl Outputs/SouthAmerica_dataset_QCed Outputs/SouthAmerica_dataset_QCed.fam 0.0884 0.01 Outputs/SouthAmerica_dataset_QCed
~/src/plink1.9/plink --bfile Outputs/SouthAmerica_dataset_QCed --remove Outputs/SouthAmerica_dataset_QCed.RemoveKin0.0884 --make-bed --out Outputs/SouthAmerica_dataset_QCed.KinshipFilter
