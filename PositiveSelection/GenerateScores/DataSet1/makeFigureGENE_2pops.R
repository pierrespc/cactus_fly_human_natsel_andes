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

system("mkdir GeneSets")
for(stat in c("mean","median")){
	for(grou in c(1,4)[1]){
  		for(fl in c(0,5,10)){
  		  
    			#pdf(paste("2pops.Scatter.",stat,".",fl,"KB.",grou,"Groups.pdf",sep=""))
  		   postscript(paste("2pops.Scatter.",stat,".",fl,"KB.",grou,"Groups.eps",sep=""),width = 7,height=7)
  		    listGenes<-c()
  		    
      		refWichi<-read.table(paste("Andean.refWichi.DetoxGenes.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),stringsAsFactors = F,header=T)
      		
      		write.table(refWichi$ENSID[!is.na(refWichi$P_Zf)],paste("GeneSets/Studied.Andean.refWichi.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),col.names = F,row.names = F,quote=F)
      		write.table(refWichi$ENSID[!is.na(refWichi$P_Zf) & refWichi$P_Zf<0.05],paste("GeneSets/Signif.Andean.refWichi.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),col.names = F,row.names = F,quote=F)
      		refWichi<-refWichi[,c("HGNC","P_Zf")]
      		names(refWichi)<-c("HGNC","PrefWichi")
   				refWichi$PrefWichi<--log10(refWichi$PrefWichi)
      		refWichi$PrefWichi[is.na(refWichi$PrefWichi)]<--0.1
      		
      		
      		refYanesha<-read.table(paste("Andean.refYanesha.DetoxGenes.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),stringsAsFactors = F,header=T)
      		write.table(refYanesha$ENSID[!is.na(refYanesha$P_Zf)],paste("GeneSets/Studied.Andean.refYanesha.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),col.names = F,row.names = F,quote=F)
      		write.table(refYanesha$ENSID[!is.na(refYanesha$P_Zf) & refYanesha$P_Zf<0.05],paste("GeneSets/Signif.Andean.refYanesha.",fl,"KB.Fisher.",stat,".",grou,"Groups.tsv",sep=""),col.names = F,row.names = F,quote=F)
      		
      		refYanesha<-refYanesha[,c("HGNC","P_Zf")]
      		names(refYanesha)<-c("HGNC","PrefYanesha")
      		refYanesha$PrefYanesha<--log10(refYanesha$PrefYanesha)
      		refYanesha$PrefYanesha[is.na(refYanesha$PrefYanesha)]<--0.1
      		
		merged<-merge(refWichi,refYanesha,by="HGNC")
		merged<-merge(merged,sens[,c("HGNC","pch")],by="HGNC")
  
      		merged$Color="grey50"
      		merged$Color[ merged$PrefWichi<0 | merged$PrefYanesha<0]<-"white"
      		merged$Color[ merged$PrefWichi>-log10(0.05) | merged$PrefYanesha>-log10(0.05)]<-"black"
 		merged$cex=ifelse(merged$Color=="black",4,2)
 		print(table(merged$Color)) 
      		min=min(c(merged$PrefWichi,merged$PrefYanesha))
      		max=max(c(merged$PrefWichi,merged$PrefYanesha))
      		mergedText=merged[ merged$PrefWichi>-log10(0.06)| merged$PrefYanesha>-log10(0.06),]
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
          		 labs(x="-log10(P-value) with Wichi as reference")+labs(y="-log10(P-value) with Yanesha as reference",title = "")+
			 labs(title="B. Genotype data set #2")+
         	 	 geom_point(aes(x=merged$PrefWichi,y=merged$PrefYanesha),colour=merged$Color,shape=merged$pch,size=merged$cex,stroke=2)+
          		 geom_text_repel(data=mergedText,x =mergedText$PrefWichi, y = mergedText$PrefYanesha,label=mergedText$HGNC,
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

