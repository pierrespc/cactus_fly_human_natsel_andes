#!/bin/Rscript


setwd("/Users/pierrespc/Documents/PostDoc/DetoxFly/LAST/SummaryTables/MaoMoreno_Gnecchio_noUros")
out<-read.table("../../../OrthoGenes.txt",stringsAsFactors=F,header=F)
names(out)<-c("HGNC","ENSGID","chr","start","end")
stat="median"
kb="10"

numSET=0
for(set in c("MaoMoreno","Gnecchio_noUros")){
	numSET=numSET+1
	iHS<-read.table(paste("../../",set,"/IHS/Andean.DetoxGenes.",kb,"KB.absIHS.",stat,".1Groups.tsv",sep=""),stringsAsFactors=F,header=T)[,c("HGNC","NSNPs",stat,paste("P_",stat,sep=""))]
	names(iHS)[c(2:4)]<-paste("IHS",paste("set",numSET,sep=""),names(iHS)[c(2:4)],sep="_")
        out<-merge(out,iHS,by="HGNC")
	if(set=="MaoMoreno"){
		pops=c("Bari","Yukpa")
	}else{
		pops=c("Yanesha","Wichi")
	}
	for(pop in pops){
		LRT<-read.table(paste("../../",set,"/TreeSelect/ref",pop,"/Results/Andean.DetoxGenes.",kb,"KB.",stat,".1Groups.tsv",sep=""),stringsAsFactors=F,header=T)[,c("HGNC","NSNPs",stat,paste("P_",stat,sep=""))]
		names(LRT)[c(2:4)]<-paste("LRT",paste("set",numSET,sep=""),pop,names(LRT)[c(2:4)],sep="_")

		out<-merge(out,LRT,by="HGNC")
		Fisher<-read.table(paste("../../",set,"/Fisher/Andean.ref",pop,".DetoxGenes.",kb,"KB.Fisher.",stat,".1Groups.tsv",sep=""),stringsAsFactors=F,header=T)[,c("HGNC","Zf","P_Zf")]
		names(Fisher)[c(2:3)]<-paste(names(Fisher)[c(2:3)],paste("set",numSET,sep=""),pop,sep="_")
		 out<-merge(out,Fisher,by="HGNC")
	}
}
write.table(out,paste("SummaryPerGene.",stat,".",kb,"KB.tsv",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
	

