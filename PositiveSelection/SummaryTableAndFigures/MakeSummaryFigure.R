#!/bin/Rscript


require(stringr)
setwd("~/Documents/PostDoc/DetoxFly/LAST/SummaryTables/MaoMoreno_Gnecchio_noUros/")
mean<-read.table("SummaryPerGene.mean.10KB.tsv",stringsAsFactors = F,header=T)
mean$NumSignal<-apply(mean[,grepl("P_Zf",names(mean)) ]<0.05,1,sum,na.rm=T)
mean$sumNA<-apply(is.na(mean[,grepl("P_Zf",names(mean)) ]),1,sum)
row.names(mean)<-mean$HGNC
names(mean)<-str_replace_all(names(mean),"_set1","")
names(mean)<-str_replace_all(names(mean),"_set2","")

#mean<-mean[order(mean$NumSignal,decreasing = T),]
median<-read.table("SummaryPerGene.median.10KB.tsv",stringsAsFactors = F,header=T)
median$NumSignal<-apply(median[,grepl("P_Zf",names(median)) ]<0.05,1,sum,na.rm=T)
median$sumNA<-apply(is.na(median[,grepl("P_Zf",names(median)) ]),1,sum)

row.names(median)<-median$HGNC

median<-median[ row.names(mean),]
names(median)<-str_replace_all(names(median),"_set1","")
names(median)<-str_replace_all(names(median),"_set2","")




meanBU<-mean
medianBU<-median
ColorManifest<-read.table("../../CheckStructure/ColorManifest.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
color<-ColorManifest$Color[ ColorManifest$Population %in% c("Wichi","Bari","Yukpa","Yanesha")]
names(color)<-ColorManifest$Population[ ColorManifest$Population %in% c("Wichi","Bari","Yukpa","Yanesha")]
color["Wichi"]<-"deeppink"


genes<-unique(c(row.names(meanBU)[ meanBU$NumSignal>0],row.names(medianBU)[ medianBU$NumSignal>0]))
for(set in c("signif","noSignif")){
  pdf(paste("summary.",set,".svg",sep=""))
  if(set=="signif"){
    mean<-meanBU[genes,]
    median<-medianBU[genes,]
  }else{
    mean<-meanBU[ ! row.names(meanBU) %in% genes,]
    median<-medianBU[ ! row.names(medianBU) %in% genes,]
  }
  plot(0,0,"n",axes=F,
     xlim=c(0,nrow(mean))+0.5,
     ylim=c(ifelse(sum(c(mean$sumNA,median$sumNA))>0,-0.1,0),
            max(cbind(mean[,grepl("Zf",names(mean)) & !grepl("P_Zf",names(mean))],
                      median[,grepl("Zf",names(median)) & !grepl("P_Zf",names(median))],
                      qchisq(0.95,4)),na.rm=T)),
     main="",xlab="",ylab="Fisher's combination score")
  rect(xleft = 0,xright=nrow(mean)+1,ybottom=0,ytop = qchisq(0.95,4),border=NA,col="lightgrey")
  abline(h=0)
  axis(2,seq(0,20,2))
  x=0
  for(i in row.names(median)){
    abline(v=x+0.5,col="darkgrey")
    x=x+1
  
    axis(1,at=x,labels = i,las=ifelse(set=="signif",2,2),tick = F,col = ifelse(x%%2==1,"grey","black"))
    buf<-0
    print(i)
    for(pop in names(color)){
      buf=buf+0.05
      points(x-buf , ifelse(is.na(mean[i,paste("Zf_",pop,sep="")]),-0.1,mean[i,paste("Zf_",pop,sep="")]),bg=color[pop],pch=21,col="black",cex=1.3)
      points(x+buf , ifelse(is.na(median[i,paste("Zf_",pop,sep="")]),-0.1,median[i,paste("Zf_",pop,sep="")]),bg=color[pop],pch=22,col="black",cex=1.3)
    }
  }
  abline(h=qchisq(0.95,4),lty=2,lwd=2)
  abline(v=x+0.5,col="darkgrey")
  if(set=="signif"){
    legend("topright",pch=c(15,16,rep(NA,6)),text.col=c(rep(c("black",NA),each=2),color),legend = c("Median","Mean",NA,NA,names(color)),ncol=2,bg="white")
  }
  dev.off()
}

