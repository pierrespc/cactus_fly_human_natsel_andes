#!/bin/Rscript

require(distances)
require(cluster)



change<-function(string,split,pos){
  return(strsplit(string,split=split)[[1]][pos])
}

params<-commandArgs(trailingOnly=T)
setwd(params[1])
eigenvector<-read.table(paste(params[2],".evec",sep=""),stringsAsFactors = F,header=F,comment.char = "#",skip=1)
eigenval<-read.table(paste(params[2],".eval",sep=""),stringsAsFactors = F,header=F,sep="\t")


colInfo<-read.table("../../ColorManifest.tsv",stringsAsFactors = F,header=T,comment.char = "@",sep="\t")

numPC=7
val<-eigenval[c(1:numPC),1]
weig=val/sum(val)
vec=eigenvector[,c(2:(numPC+1))]
rownames(vec)<-eigenvector$V1
dm<-data.frame(as.matrix(daisy(vec,metric="gower",weights=weig)))
row.names(dm)<-names(dm)<-rownames(vec)
mds<-cmdscale(dm)
mds<-data.frame(mds)
mds$Ind<-sapply(row.names(mds),change,split=":",pos=2)
mds$Population<-sapply(row.names(mds),change,split=":",pos=1)
mds<-merge(mds,colInfo,by="Population")
pdf(paste("MDS-",numPC,"PCs.pdf",sep=""))
plot(mds$X1[mds$cex==min(mds$cex)],mds$X2[mds$cex==min(mds$cex)],col="black",bg=mds$Color[mds$cex==min(mds$cex)],cex=mds$cex[mds$cex==min(mds$cex)],pch=mds$Point,xlim=c(min(mds$X1),max(mds$X1)),ylim=c(min(mds$X2),max(mds$X2)),main=paste("MDS for ",numPC,"-first PC based weighted euclidean distance",sep=""),xlab="1st Dimension",ylab="2nd Dimension")
for(cex in unique(mds$cex[mds$cex!=min(mds$cex)])){
    points(mds$X1[mds$cex==cex],mds$X2[mds$cex==cex],col=mds$Color[mds$cex==cex],pch=mds$Point[mds$cex==cex],cex=mds$cex[mds$cex==cex])
  }

  
#plot(0,0,"n",ann=F,axes=F)
colInfo<-colInfo[ colInfo$Population %in% mds$Population,]
colInfo<-colInfo[ order(as.numeric(colInfo$cex),colInfo$Region,colInfo$Population),]
legend("bottomleft",pch=colInfo$Point,col="black",pt.bg=colInfo$Color,legend=paste(colInfo$Population,": ",colInfo$Region,sep=""),ncol = 1,cex=0.8)

val<-eigenval[,1]
val<-100*val/sum(val)
plot(val[c(1:20)],ylab="% variance explained",xlab="PC",pch=3)
points(val[c(1:20)],type="l")
points(numPC,val[numPC],cex=2,lwd=2,col="red")



dev.off()


