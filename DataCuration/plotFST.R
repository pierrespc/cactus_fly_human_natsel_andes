#!/bin/Rscript

require(gplots)

param<-commandArgs(trailingOnly=T)
setwd(param)


orderReg<-c(1:6)
names(orderReg)<-c("Andes (set 1)","Andes (set 2)",
            "Equatorial-Tucanoan (set 2)",
            "Ge-Pano-Carib (set 1)","Ge-Pano-Carib (set 2)",
            "Chibchan-Paezan (set 1)")


ColorManifest<-read.table("../../ColorManifest.tsv",stringsAsFactors = F,header=T,
                          sep="\t",comment.char = "@")

pop<-read.table("FST.PopIndexes",stringsAsFactors = F,header=F)
pop$V1=pop$V1+1
Fst<-read.table("FST.Matrix",header = F,skip=1,stringsAsFactors = F)
Fst<-Fst[,-1]
Fst<-Fst/1000
names(Fst)<-row.names(Fst)<-pop$V2

SD<-read.table("FST.Std.Matrix",header = F,skip=1,stringsAsFactors = F)
SD<-SD[,-1]
SD<-SD/1e6
names(SD)<-row.names(SD)<-pop$V2

outTab<-c()

for(p1 in c(1:(nrow(Fst)-1))){
  for(p2 in c((p1+1):nrow(Fst))){
      outTab<-rbind(outTab,cbind("p1"=names(Fst)[p1],
                               "r1"=ColorManifest$Region[ ColorManifest$Population==names(Fst)[p1]],
                               "p2"=names(Fst)[p2],
                               "r2"=ColorManifest$Region[ ColorManifest$Population==names(Fst)[p2]],
                               "Fst"=Fst[p1,p2],
                               "SD"=SD[p1,p2]))
  }
}

outTab<-as.data.frame(outTab,stringsAsFactors = F)
outTab$Fst<-as.numeric(outTab$Fst)
outTab$SD<-as.numeric(outTab$SD)
#outTab<-outTab[order(outTab$r1,outTab$r2,outTab$Fst),]
outTab<-outTab[order(outTab$Fst),]
color<-list()
for(i in unique(c(outTab$p1,outTab$p2))){
  color[i]<-ColorManifest$Color[ColorManifest$Population==i]
}


labelName<-list()
count=0

orderCol=c()
for(i in c(1:nrow(outTab))){
  if(!outTab$p1[i] %in% names(labelName)){
    count=count+1
    labelName[[outTab$p1[i]]]<-count
    orderCol<-c(orderCol,color[[outTab$p1[i]]])

  }
  if(!outTab$p2[i] %in% names(labelName)){
    count=count+1
    labelName[[outTab$p2[i]]]<-count
    orderCol<-c(orderCol,color[[outTab$p2[i]]])
  }
}


for(i in c(1:nrow(outTab))){
  if(orderReg[outTab$r1[i]] > orderReg[outTab$r2[i]]){
    tmpP<-outTab$p1[i]
    tmpR<-outTab$r1[i]
    outTab$p1[i]<-outTab$p2[i]
    outTab$r1[i]<-outTab$r2[i]
    outTab$p2[i]<-tmpP
    outTab$r2[i]<-tmpR
  }
}

    
  

pdf("Fst.pairwise.pdf")
#heatmap.2(as.matrix(Fst),col=gray.colors(200),density.info = "none",trace="none",main = "Fst heatmap")
#plot(0,0,"n",axes=F,ann=F)
#text(-1,seq(1,-1,length.out = length(orderCol)),paste(unlist(labelName),":",names(labelName)),col = orderCol,pos=4)

plot(0,0,xlim=c(0,nrow(outTab)),ylim=c(min(c(0,outTab$Fst-3*outTab$SD)),max(outTab$Fst+9*outTab$SD)),ann=F,axes=F,"n")
points(c(1:nrow(outTab)),outTab$Fst,pch=15)
segments(x0 = c(1:nrow(outTab)),y0=outTab$Fst-3*outTab$SD,y1=outTab$Fst+3*outTab$SD,pch=15)
abline(h=0)
xseq<-seq(0.5,nrow(outTab)+0.5,1)

for(i in c(1:nrow(outTab))){
  #rect(xleft = xseq[i],
  #   xright = xseq[i+1],
  #   ybottom = max(outTab$Fst+3*outTab$SD),
  #   ytop =  max(outTab$Fst+6*outTab$SD),
  #   col = color[[outTab$p1[i]]])
  #rect(xleft = xseq[i],
  #     xright = xseq[i+1],
  #     ybottom = max(outTab$Fst+6*outTab$SD),
  #     ytop =  max(outTab$Fst+9*outTab$SD),
  #     col = color[[outTab$p2[i]]])
  #text(x=i,y=max(outTab$Fst+4.5*outTab$SD),labelName[[outTab$p1[i]]],cex=0.6)
  #text(x=i,y=max(outTab$Fst+7.5*outTab$SD),labelName[[outTab$p2[i]]],cex=0.6)
  points(x=i,y=max(outTab$Fst+4.5*outTab$SD),pch=ColorManifest$Point[ColorManifest$Population==outTab$p1[i]],col="black",bg=ColorManifest$Color[ColorManifest$Population==outTab$p1[i]])
  points(x=i,y=max(outTab$Fst+7.5*outTab$SD),pch=ColorManifest$Point[ColorManifest$Population==outTab$p2[i]],col="black",bg=ColorManifest$Color[ColorManifest$Population==outTab$p2[i]])
  

}


ColorManifest<-ColorManifest[ ColorManifest$Population %in% c(outTab$p1,outTab$p2),]
ColorManifest<-ColorManifest[order(ColorManifest$Region),]
segments(y0 = 0,y1=max(outTab$Fst+3*outTab$SD),x0 = xseq,lwd=0.5)
axis(2)
title(xlab = "Pairwise comparison",ylab="Fst",main="Fst")
legend(x = nrow(outTab),y = 0,xjust = 1,yjust=0,legend=paste(ColorManifest$Population,"; ",ColorManifest$Region,sep=""),
       pch=ColorManifest$Point,pt.bg = ColorManifest$Color,col = "black",bg="white",cex=0.7)
dev.off()

