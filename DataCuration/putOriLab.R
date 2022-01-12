#!/bin/Rscript

params<-commandArgs(trailingOnly=T)


fam<-read.table(params[1],stringsAsFactors=F,header=F)
print(table(fam$V1))

ori<-read.table(params[2],stringsAsFactors=F,header=F)

fam$V6<-fam$V1
whichAndean<-which(fam$V1=="Andean")
for(i in whichAndean){
	fam$V1[i]<-fam$V6[i]<-ori$V1[ ori$V2==fam$V2[i]]
}

print(table(fam$V1))

write.table(fam,params[1],col.names=F,row.names=F,sep=" ",quote=F)


