#!/bin/Rscript

fam<-read.table("1KG/SouthAmerica_dataset_QCed.KinshipFilter.1KG.fam",stringsAsFactors=F,header=F)

pop<-read.table("/Volumes/MARIOLOLO/PolymorphismDataFromLitterature/1KG_Phase3/integrated_call_samples_v3.20130502.ALL.panel",stringsAsFactors=F,header=T)

for(i in c(1:nrow(fam))){
	lab<-pop$pop[ pop$sample==fam$V2[i]]
	if(length(lab)==0){
		fam$V6[i]=fam$V1[i]
	}else{
		fam$V6[i]=lab
		fam$V1[i]=lab
	}
}

table(fam$V6)
write.table(fam,"1KG/SouthAmerica_dataset_QCed.KinshipFilter.1KG.fam",col.names=F,row.names=F,quote=F)
