#!/bin/Rscript

library(stringr)
library(ggrepel)
library(ggplot2)
library(dplyr)
require(maps)
require(mapdata)
library(ggplot2)
setwd("~/Documents/PostDoc/DetoxFly/LAST/MaoMoreno/Fisher/")


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

system("mkdir GeneSets")
for(stat in c("mean","median")){
	for(grou in c(1,4)[1]){
  		for(fl in c(0,5,10)){
  		  
    			#pdf(paste("2pops.Scatter.",stat,".",fl,"KB.",grou,"Groups.pdf",sep=""))
  		   postscript(paste("2pops.Scatter.",stat,".",fl,"KB.",grou,"Groups.eps",sep=""),width = 7,height=7)
  		    listGenes<-c()
  		    
      		refBari<-read.table(paste("Andean.refBari.DetoxGenes.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),stringsAsFactors = F,header=T)
      		
      		write.table(refBari$ENSID[!is.na(refBari$P_Zf)],paste("GeneSets/Studied.Andean.refBari.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),col.names = F,row.names = F,quote=F)
      		write.table(refBari$ENSID[!is.na(refBari$P_Zf) & refBari$P_Zf<0.05],paste("GeneSets/Signif.Andean.refBari.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),col.names = F,row.names = F,quote=F)
      		refBari<-refBari[,c("HGNC","P_Zf")]
      		names(refBari)<-c("HGNC","PrefBari")
   				refBari$PrefBari<--log10(refBari$PrefBari)
      		refBari$PrefBari[is.na(refBari$PrefBari)]<--0.1
      		
      		
      		refYukpa<-read.table(paste("Andean.refYukpa.DetoxGenes.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),stringsAsFactors = F,header=T)
      		write.table(refYukpa$ENSID[!is.na(refYukpa$P_Zf)],paste("GeneSets/Studied.Andean.refYukpa.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),col.names = F,row.names = F,quote=F)
      		write.table(refYukpa$ENSID[!is.na(refYukpa$P_Zf) & refYukpa$P_Zf<0.05],paste("GeneSets/Signif.Andean.refYukpa.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),col.names = F,row.names = F,quote=F)
      		
      		refYukpa<-refYukpa[,c("HGNC","P_Zf")]
      		names(refYukpa)<-c("HGNC","PrefYukpa")
      		refYukpa$PrefYukpa<--log10(refYukpa$PrefYukpa)
      		refYukpa$PrefYukpa[is.na(refYukpa$PrefYukpa)]<--0.1
      		
		merged<-merge(refBari,refYukpa,by="HGNC")
		merged<-merge(merged,sens[,c("HGNC","pch")],by="HGNC")
  
      		merged$Color="grey50"
      		merged$Color[ merged$PrefBari<0 | merged$PrefYukpa<0]<-"white"
      		merged$Color[ merged$PrefBari>-log10(0.05) | merged$PrefYukpa>-log10(0.05)]<-"black"
 		merged$cex=ifelse(merged$Color=="black",4,2)
 		print(table(merged$Color)) 
      		min=min(c(merged$PrefBari,merged$PrefYukpa))
      		max=max(c(merged$PrefBari,merged$PrefYukpa))
      		mergedText=merged[ merged$PrefBari>-log10(0.06)| merged$PrefYukpa>-log10(0.06),]
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
          		 labs(x="-log10(P-value) with Bari as reference")+labs(y="-log10(P-value) with Yukpa as reference",title = "")+
         	 	 labs(title="A. Genotype data set #1")+
			 geom_point(aes(x=merged$PrefBari,y=merged$PrefYukpa),colour=merged$Color,shape=merged$pch,size=merged$cex,stroke=2)+
          		 geom_text_repel(data=mergedText,x =mergedText$PrefBari, y = mergedText$PrefYukpa,label=mergedText$HGNC,
			        colour=mergedText$Color,
			        #size=mergedText$cex,
			        segment.color=NA,
                        	#colour="lightseagreen",
                        	size=5,
                        	fontface = "italic",
                        	box.padding = 0.5,
                        	point.padding = 0.1,
                        	seed = 1
			        #stroke=2,					
                        	#segment.color = "lightseagreen"
			)
        	print(plo)
    		dev.off()
    			
  		}
	}
}

