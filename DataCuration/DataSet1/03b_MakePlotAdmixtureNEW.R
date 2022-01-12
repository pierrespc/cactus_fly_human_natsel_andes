#!/usr/bin/Rscript

library(RColorBrewer)
library(stringr)
setwd("/Users/pierrespc/Documents/PostDoc/DetoxFly/Gnecchio/Reformate/Outputs/1KG/Admixture/BestRUNperK//")
orderRegions<-c("Africa","AfroAmericans","Europe","Admixed",
                "Northern-Amerind","Equatorial-Tucanoan","Ge-Pano-Carib","Andes")

#orderINDS<-read.table("ORDERK8.IND.txt",stringsAsFactors = F,header=F)$V1
orderINDS<-NULL
fam<-read.table("../../SouthAmerica_dataset_QCed.KinshipFilter.1KG.fam",stringsAsFactors=F,header=F)
regions<-read.table("../../../../ColorManifest.tsv",stringsAsFactors=F,header=T,comment.char = "@",sep="\t")

adjustTextRight=300
pdf("AllAdmixure.pdf", width=180,height=20)
for(KNUM in c(8)){
  a<-read.table(paste("SouthAmerica_dataset_QCed.KinshipFilter.1KG.pruned.",KNUM,".Q",sep=""),stringsAsFactors=F,header=F)
	outfile<-paste("SouthAmerica_dataset_QCed.KinshipFilter.1KG.pruned.K",KNUM,sep="")
  numberK=dim(a)[2]
  numberInd=dim(a)[1]

  #MyColors=c("green4","darkolivegreen3","green3","mediumseagreen")
  MyColors=c(brewer.pal(12,"Paired"),"darkolivegreen","limegreen","turquoise")
  if(numberK==2){
	#  MyColors=MyColors[c(3,2)]
  }
  if(numberK==3){
   MyColors=MyColors[c(6,4,2)]
  }
  if(numberK==4){
    MyColors=MyColors[c(6,2,4,7)]
  }
  if(numberK==5){
    MyColors=MyColors[c(4,7,6,1,2)]
  }
  if(numberK==6){
    MyColors=MyColors[c(3,7,6,2,1,4)]
  }
  if(numberK==7){
    MyColors=MyColors[c(7,3,2,4,8,6,13)]
  }
  if(KNUM==8){
   MyColors=MyColors[c(4,1,6,7,2,13,3,8)]
  }
  if(KNUM==9){
    MyColors=MyColors[c(13,6,1,3,12,7,2,4,8)]

  }
  if(KNUM==10){
    MyColors=MyColors[c(12,2,13,6,4,14,8,3,1,7)]
  }
  if(KNUM==11){
    MyColors=MyColors[c(12,2,8,5,3,14,7,1,4,6,13)]
   # MyColors=c(MyColors[c(2,3,7,8)],"paleturquoise2",MyColors[6],"darkolivegreen",MyColors[1],"goldenrod",MyColors[c(5,4)])
  }
  if(KNUM==12){
    MyColors=MyColors[c(12,15,3,4,8,14,6,13,2,1,7,5)]
  }
  if(KNUM==13){
    MyColors=MyColors[c(8,13,7,3,6,12,1,11,4,5,2,15,14)]
  }
  a$Ind<-fam$V2
  a$Population<-fam$V1
  
  a<-merge(a,regions,by="Population")

  #a$Region[ is.na(a$Region)]<-a$

  numPops<-length(unique(a$Population))
  numRegions<-length(unique(a$Region))

  if(dim(a)[1] != numberInd){
	  stop("your regions file and admixture output do not coincide: do not have the same number of Pops")
  }
  a<-a[ ! a$Population %in% c("MAYAN","NAHUAN"),]
 
  out<-c()
  for(reg in orderRegions){
    print(reg)
    temp<-a[ a$Region == reg,]
    if(dim(temp)[1]==0){
      print(paste(reg,"skipped"))
      
      next
    }	
	  meanOverRegion<-vector(length=numberK)
	  names(meanOverRegion)<-paste("V",c(1:numberK),sep="")
	
	
	  meanOverPop<-data.frame(matrix(NA,length(unique(temp$Population)),numberK))
	  names(meanOverPop)<-paste("V",c(1:numberK),sep="")
	  rownames(meanOverPop)<-unique(temp$Population)
	
	  for(K in c(1:numberK)){
		  meanOverRegion[paste("V",K,sep="")]<-mean(temp[,paste("V",K,sep="")])
		  for(pop in unique(temp$Population)){
			  meanOverPop[pop,paste("V",K,sep="")]<-mean(temp[temp$Population==pop,paste("V",K,sep="")])
		  }
	  }
	
	  meanOverRegion<-meanOverRegion[order(as.numeric(meanOverRegion),decreasing=T)]
	  meanOverPop<-meanOverPop[,names(meanOverRegion)]
	  meanOverPop<-meanOverPop[order(meanOverPop[,1],decreasing=T),]
	  for(pop in rownames(meanOverPop)){
		  temp2<-a[ a$Population == pop,]
		  temp2<-temp2[order(temp2[,names(meanOverRegion)[1]],decreasing=T),]
		  out<-rbind(out,temp2)
	  }
  }
  row.names(out)<-out$Ind
  if(! is.null(orderINDS)){
    out<-out[orderINDS,]
    out<-out[ ! is.na(out$Ind),]
    if(nrow(out)!=numberInd){
      stop("pb making out")
    }
  }
  
  
  separ<-ceiling(numberInd/500)
  #separPop<-ceiling(numberInd/1000)
  separPop=1

  #pdf(paste(outfile,".pdf",sep=""), width=180,height=20)
  par(mar=c(15, 2, 40, 2) + 0.1)

  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(dim(a)[1]+separPop*(numPops-numRegions)+numRegions*(separ))),main=paste("K =", numberK),ann=F,axes=F)
  title(paste("K =",numberK,sep=""))
  xleft=0
  atPop<-c()


  dimPrevRegion<- 0
  atRegion<-c()


  ITER = 0
  for(reg in orderRegions){
	  temp<-out[ out$Region == reg,]
	  Population=temp$Population[1]
	  ITER=ITER+1
	  #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
	  axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,tick=T)
	  #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
	  for(ind in unique(temp$Ind)){
		  ybottom=0
		  if(temp$Population[ temp$Ind==ind] != Population){
			  Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
			  xleft=xleft+separPop
			  ITER=ITER+1
			
	  		#axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
			  axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,tick=T)
		  	#axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
			  dimPrevRegion=dimPrevRegion+separPop
		  }
		  for(k in c(1:numberK)){
			  ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
			  rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
			  ybottom=ytop
		  }
		  xleft=xleft+1
	  }
	  xleft=xleft+separ
	  atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
	  dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }

  axis(1,at=atRegion,label=names(atRegion),cex.axis=5,tick=F,line=4)
  text(x=-30,y=0.5,paste("K = ",KNUM,sep=""),cex=12)
  forLeg<-unique(out[,c("Population","Color")])
  #plot(0,0,"n",pch=0,ann=F,axes=F)
  #legend("topleft",legend=paste(c(1:dim(forLeg)[1]),": ",forLeg$Population,sep=""),ncol=10,cex=8,text.col=forLeg$Color)

  
  #dev.off()
  
  meanByPop=data.frame(matrix(NA,0,KNUM))
  names(meanByPop)=paste(c(1:KNUM),sep="")
  iterP=0
  for(region in unique(out$Region)){
    iterP=iterP+1
    meanByPop[iterP,]<-apply(out[out$Region==region,paste("V",c(1:KNUM),sep="")],2,mean)
    rownames(meanByPop)[iterP]<-region
    for(pop in unique(out$Population[out$Region==region])){
      iterP=iterP+1
      meanByPop[iterP,]<-apply(out[out$Population==pop,paste("V",c(1:KNUM),sep="")],2,mean)
      rownames(meanByPop)[iterP]<-pop
    }
  }
  iterP=iterP+1
  #write.table(meanByPop,paste("chr1.1KG.PopArg.pruned.",KNUM,".MeanByGroup.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  write.table(meanByPop,paste(outfile,".MeanByGroup.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  #  write.table(out,paste("chr1.1KG.PopArg.pruned.",KNUM,".AncestryComponentByIndividual.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)  
  write.table(out,paste(outfile,".AncestryComponentByIndividual.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  
}
dev.off()

CV<-list()
a<-read.table("ListBESTruns",stringsAsFactors = F,header=T)
a<-a[ a$K >2,]
for(run in c(1:10)){
	cv<-read.table(paste("../RUN",run,"/SouthAmerica_dataset_QCed.KinshipFilter.1KG.pruned.CV",sep=""),stringsAsFactors = F,header=T)
	for(k in c(3:10)){
	  if(run==1){
	    CV[[k]]<-cv$CVscore[cv$K==k]
	  }else{
	    CV[[k]]<-c(CV[[k]],cv$CVscore[cv$K==k])
	  }
	}
}
pdf("CrossValidation.pdf")
boxplot(CV,xlim=c(2.5,10.5),main="Cross-Validation Score",axes=F)
points(CVscore~K,data=a,type="l",col="red")
axis(2,at=seq(round(min(a$CVscore),digits = 4),round(max(a$CVscore),digits = 4),by = 2e-4))
axis(1,at=a$K)
points(CVscore~K,data=a,pch=4,col="red")
#points(a$K[ a$CVscore==min(a$CVscore)],a$CVscore[ a$CVscore==min(a$CVscore)],col="red",pch=1,cex=3,lwd=2)
dev.off()
