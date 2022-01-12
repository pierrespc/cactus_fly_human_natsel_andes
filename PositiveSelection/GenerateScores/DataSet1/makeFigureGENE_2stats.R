#!/bin/Rscript

library(stringr)
library(ggrepel)
library(ggplot2)
library(dplyr)
require(maps)
require(mapdata)
library(ggplot2)
setwd("~/Documents/PostDoc/DetoxFly/LAST/Gnecchio_noUros/Fisher/")


source("~/Documents/PostDoc/scripts/R/MultiGGPLOT2.R")
up<-read.table("../../../UpGenes.txt",stringsAsFactors = F,header=F)
names(up)<-"HGNC"
up$Up<-T
down<-read.table("../../../DownGenes.txt",stringsAsFactors = F,header=F)
names(down)<-"HGNC"
down$Down=T

sens<-merge(up,down,by="HGNC",all=T)
sens$Down[ is.na(sens$Down)]<-F
sens$Up[ is.na(sens$Up)]<-F

sens$pch[ sens$Up]<-2
sens$pch[ sens$Down]<-6
sens$pch[ sens$Down & sens$Up ]<-5

BASE=F


colors=c("Wichi"="deeppink","Yanesha"="#9751A1")
system("mkdir GeneSets")
for(pop in c("Yanesha","Wichi")){
	for(grou in c(1,4)[1]){
  		for(fl in c(0,5,10)){
  		  
    			#pdf(paste("2pops.Scatter.",stat,".",fl,"KB.",grou,"Groups.pdf",sep=""))
  		    postscript(paste("2Stats.Scatter.ref",pop,".",fl,"KB.",grou,"Groups.eps",sep=""),width = 7,height=7)
  		    listGenes<-c()
  		    
      		refmean<-read.table(paste("Andean.ref",pop,".DetoxGenes.",fl,"KB.Fisher.mean.",grou,"Groups.tsv",sep=""),stringsAsFactors = F,header=T)
      		refmean<-refmean[,c("HGNC","P_Zf")]
      		names(refmean)<-c("HGNC","Pmean")
   				refmean$Pmean<--log10(refmean$Pmean)
      		refmean$Pmean[is.na(refmean$Pmean)]<--0.1
      		
      		
      		refmedian<-read.table(paste("Andean.ref",pop,".DetoxGenes.",fl,"KB.Fisher.median.",grou,"Groups.tsv",sep=""),stringsAsFactors = F,header=T)
      		
      		refmedian<-refmedian[,c("HGNC","P_Zf")]
      		names(refmedian)<-c("HGNC","Pmedian")
      		refmedian$Pmedian<--log10(refmedian$Pmedian)
      		refmedian$Pmedian[is.na(refmedian$Pmedian)]<--0.1
      		
		merged<-merge(refmean,refmedian,by="HGNC")
		merged<-merge(merged,sens[,c("HGNC","pch")],by="HGNC")
  
      		merged$Color="black"
      		merged$Color[ merged$Pmean<0 | merged$Pmedian<0]<-"white"
      		merged$Color[ merged$Pmean>-log10(0.05) | merged$Pmedian>-log10(0.05)]<-colors[pop]
 		merged$cex=ifelse(merged$Color==colors[pop],4,2)
 		print(table(merged$Color)) 
      		min=min(c(merged$Pmean,merged$Pmedian))
      		max=max(c(merged$Pmean,merged$Pmedian))
      		mergedText=merged[ merged$Pmean>-log10(0.06)| merged$Pmedian>-log10(0.06),]
        	plo<-ggplot()+coord_cartesian(xlim=c(min,max),ylim=c(min,max))+
          		 geom_rect(mapping=aes(xmin=-100,xmax=0,ymin = -100,ymax=100),fill="grey10",color=NA)+
          		 geom_rect(mapping=aes(xmin=-100,xmax=100,ymin = -100,ymax=0),fill="grey10",color=NA)+
          		 geom_rect(mapping=aes(xmin=0,xmax=-log10(0.05),ymin =0,ymax=-log10(0.05)),fill="grey75",color = NA) +
          		 #geom_rect(mapping=aes(xmin=-log10(0.05),xmax=-log10(0.01),ymin =0,ymax=-log10(0.01)),fill="grey85",color = NA) +
          		 #eom_rect(mapping=aes(xmin=0,xmax=-log10(0.01),ymin =-log10(0.05),ymax=-log10(0.01)),fill="grey85",color = NA) +
          		 #geom_rect(mapping=aes(xmin=-log10(0.01),xmax=-log10(0.001),ymin =0,ymax=-log10(0.001)),fill="grey95",color = NA) +
          		 #geom_rect(mapping=aes(xmin=0,xmax=-log10(0.001),ymin =-log10(0.01),ymax=-log10(0.001)),fill="grey95",color = NA) +
          		 geom_hline(yintercept =-log10(c(0.05)),lwd=1.3,lty=2)+
          		 geom_vline(xintercept =-log10(c(0.05)),lwd=1.3,lty=2)+
        	         #geom_hline(yintercept =-log10(c(0.05,0.01,0.001)),lwd=1.3,lty=2)+
        	         #geom_vline(xintercept =-log10(c(0.05,0.01,0.001)),lwd=1.3,lty=2)+
        			  
        	         geom_abline(intercept=0,slope=1,lwd=1.3,lty=2)+  
          		 theme_bw(base_size=20,base_line_size = NA)+
          		 labs(x="-log10(P-value) from mean")+labs(y="-log10(P-value) from median",title = ifelse(pop=="Yanesha","C. with Yanesha as reference in LRT","D. with Wichi as reference in LRT"))+
         	 	 #geom_point(aes(x=merged$Pmean,y=merged$Pmedian),colour=merged$Color,shape=merged$pch,size=merged$cex,stroke=2)+
	         	 geom_point(aes(x=merged$Pmean,y=merged$Pmedian),colour="black",shape=merged$pch,size=merged$cex,stroke=2)+
          		 geom_text_repel(data=mergedText,x =mergedText$Pmean, y = mergedText$Pmedian,label=mergedText$HGNC,
			        #colour=mergedText$Color,
				colour="black",
			        #size=mergedText$cex,
			        segment.color=NA,
                        	#colour=colors[pop],
                        	size=5,
                        	fontface = "italic",
                        	box.padding = 0.5,
                        	point.padding = 0.1,
                        	seed = 1
			        #stroke=2,					
                        	#segment.color = colors[pop]
			)
        	print(plo)
    		dev.off()
    			
  		}
	}
}

