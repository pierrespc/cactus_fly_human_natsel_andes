#!/bin/Rscript



params<-commandArgs(trailingOnly = T)

#params=c("~/Documents/PostDoc/DetoxFly/RunIHS//AYMARAN.absIHS",
#         "~/Documents/PostDoc/DetoxFly/OrthoGenes.txt",
#         "~/Documents/PostDoc/DetoxFly/GeneStartEnd.BED",
#         2,3,4,
#         "median",
#         "~/Documents/PostDoc/DetoxFly/RunIHS//AYMARAN.DetoxGenes.Median.tsv",
#         "~/Documents/PostDoc/DetoxFly/RunIHS//AYMARAN.BKGDGenes.Median.tsv",
#         1)

if(length(params)!=11){
  stop("
param1: <inputScorePerSNP>
param2: <orthoGenes> with ID | ENSID | chr | start | end
param3: <backgroundGenes>  with ID | ENSID | chr | start | end
param4: column with chr of snps in <inputScorePerSNP
param5: column with pos of snps in <inputScorePerSNP
param6: column with score of snps in <inputScorePerSNP
param7: summary stat (mean,median,proportion)
param8: output file for orthogenes
param9: output file for background genes
param10: <number of gene classes according to NSNPs to compute Pvalues
param11: number of KB for flanking regions")
}

scoreFile=read.table(params[1],stringsAsFactors = F,header=T)[,as.numeric(params[c(4:6)])]
genes=read.table(params[2],stringsAsFactors = F,header=F)
bkg=read.table(params[3],stringsAsFactors = F,header=F)
names(bkg)<-c("HGNC","ENSID","chr","Start","End")
bkg$Start=bkg$Start-as.numeric(params[11])*1000
bkg$End=bkg$End+as.numeric(params[11])*1000
names(genes)<-c("HGNC","ENSID","chr","Start","End")
genes$Start=genes$Start-as.numeric(params[11])*1000
genes$End=genes$End+as.numeric(params[11])*1000

bkg<-cbind(bkg,NA)
names(bkg)[length(bkg)]<-"NSNPs"
bkg<-cbind(bkg,NA)
names(bkg)[length(bkg)]<-params[7]
bkg<-cbind(bkg,NA)
names(bkg)[length(bkg)]<-paste("P_",params[7],sep="")


for(i in c(1:nrow(bkg))){
  #print(i)
  tmp<-scoreFile[scoreFile[,1]==bkg[i,3] & scoreFile[,2]>=bkg[i,4] & scoreFile[,2]<=bkg[i,5],3]
  bkg[i,"NSNPs"]<-length(tmp)
  bkg[i,params[7]]<-eval(parse(text=paste(params[7],"(tmp)",sep="")))
}


seq<-quantile(bkg$NSNPs[ bkg$NSNPs>0],probs=seq(0,1,length.out=as.numeric(params[10])+1))
print(seq)
if(as.numeric(params[10]>1)){
  for(i in c(2:(as.numeric(params[10])))){
    bkg[bkg$NSNPs>=seq[i-1] & bkg$NSNPs<seq[i],paste("P_",params[7],sep="")]<-(rank(1/bkg[bkg$NSNPs>=seq[i-1] & bkg$NSNPs<seq[i],params[7]],ties.method = "max")+1)/(sum(bkg$NSNPs>=seq[i-1] & bkg$NSNPs<seq[i])+1)
  }
  i=i+1
}else{
  i=2
}
bkg[bkg$NSNPs>=seq[i-1] & bkg$NSNPs<=seq[i],paste("P_",params[7],sep="")]<-(rank(1/bkg[bkg$NSNPs>=seq[i-1] & bkg$NSNPs<=seq[i],params[7]],ties.method = "max")+1)/(sum(bkg$NSNPs>=seq[i-1] & bkg$NSNPs<seq[i])+1)


genes<-merge(genes,bkg[,c("ENSID","NSNPs",params[7],paste("P_",params[7],sep=""))],by="ENSID",all.x=T)


write.table(bkg,params[9],sep="\t",col.names=T,row.names=F,quote=F)
write.table(genes,params[8],sep="\t",col.names=T,row.names=F,quote=F)
