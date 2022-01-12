#!/bin/Rscript


system("mkdir ~/Documents/PostDoc/DetoxFly/LAST/Gnecchio_noUros/Fisher/")
setwd("~/Documents/PostDoc/DetoxFly/LAST/Gnecchio_noUros/Fisher/")
ortho<-read.table("../../../OrthoGenes.txt",stringsAsFactors = F,header=F)
for(stat in c("mean","median")){
	for(numG in c(1)){
		for(pop in c("Wichi","Yanesha")){
	  		#for(set in c(paste(c("DetoxGenes","BKGDGenes"),rep(c("0KB","5KB","10KB"),each=2),sep="."),"BKGDGenes.eQTLs","DetoxGenes.eQTLs")){
			#for(set in c(paste(c("DetoxGenes","BKGDGenes"),rep(c("0KB"),each=2),sep="."),"BKGDGenes.eQTLs","DetoxGenes.eQTLs")){
			for(set in paste(c("DetoxGenes","BKGDGenes"),rep(c("0KB","5KB","10KB"),each=2),sep=".")){
	    			print(c(pop,set))
	    			ihs<-read.table(paste("../IHS/Andean.",set,".absIHS.",stat,".",numG,"Groups.tsv",sep=""),stringsAsFactors = F,header=T)
	    			lrt<-read.table(paste("../TreeSelect/ref",pop,"/Results/Andean.",set,".",stat,".",numG,"Groups.tsv",sep=""),stringsAsFactors = F,header=T)
	    			if(! grepl("eQTLs",set)){
		    			names(ihs)[c(6:8)]<-c("Nsnps_absIHS","absIHS","P_absIHS")
		    			names(lrt)[c(6:8)]<-c("Nsnps_LRT","LRT","P_LRT")
		    			out<-merge(ihs,lrt,by=names(ihs)[c(1:5)])
		    			out<-out[ order(as.numeric(out$chr),as.numeric(out$Start)),]
	    			}else{
		    			names(ihs)[c(2:4)+grepl("Detox",set)]<-c("Nsnps_absIHS","absIHS","P_absIHS")
		    			names(lrt)[c(2:4)+grepl("Detox",set)]<-c("Nsnps_LRT","LRT","P_LRT")	
		    			out<-merge(ihs,lrt,by=names(ihs)[c(1)])
	    			}
	    			out$Zf=-2*(log(out$P_absIHS)+log(out$P_LRT))
	    			out$P_Zf=pchisq(out$Zf,df=4,lower.tail = F)
	    			write.table(paste("Andean.ref",pop,".",set,".Fisher.",stat,".",numG,"Groups.tsv",sep=""),col.names=F,row.names = F,quote=F,sep="\t")
	    			if(grepl("BKGDGenes",set)){
	      				png(paste("Andean.ref",pop,".VerifDistribution.Fisher.",stat,".",set,".",numG,"Groups.png",sep=""))
	      				hist1<-hist(rchisq(100000,df = 4),nclass=25,plot=F)
	      				hist1$counts=hist1$counts/sum(hist1$counts)
	      				hist2<-hist(out$Zf,nclass=25,plot=F)
	      				hist2$counts=hist2$counts/sum(hist2$counts)
      
			      		plot(hist1$mids,hist1$counts,type="l",col="black",lwd=2,
        			   		xlab="Z-score",ylab="Frequency",
						main=paste("In",pop,ifelse(grepl("eQTLs",set),"for eQTLs","for coding region")),
	           		   		xlim=c(0,max(c(hist1$mids,hist2$mids))),
	           		   		ylim=c(0,max(c(hist1$counts,hist2$counts))))
	      		      		points(hist2$mids,hist2$counts,type="l",col="red",lwd="2")
	      		      		legend("topright",lty=1,col=c("red","black"),lwd=2,legend = c("Observed","Expected"))

	    			}else{
	      				print(na.omit(out[ out$P_Zf<0.05,]))
	    			}
	    
	    			write.table(out,paste("Andean.ref",pop,".",set,".Fisher.",stat,".",numG,"Groups.tsv",sep=""),col.names = T,row.names=F,sep="\t",quote=F)
	  		}
		}
	}
}
