#!/usr/bin/Rscript
 
usage <- paste('usage: /aplic/R-2.11.1/bin/R CMD BATCH PutPopInfinPed.R <evecFile> <evalFile> <FileCodePCA> <outFile> <numEvec>\nWith\n
1. <evecFile> the .eval file from PCA\n
2. <evalFile> the .evec file from PCA\n
3. <FileCodePCA> name of file with parameters for plot\n
4. <outFile> name of output pdf file\n
5. <numEvec> number of PCs to plot to do (if N or N+1, the plots will be 1 VS 2; 3 VS 4 .... N VS N-1\n')


if (!exists("parameters")){
	parameters <- commandArgs(trailingOnly=T)
}
if (length(parameters) == 5){
	evecFile <- parameters[1]
	evalFile <- parameters[2]
	FileCodePCA <- parameters[3]
	outFile <- parameters[4]
	numEvec <- as.numeric(parameters[5])


} else{
	stop(usage)
}



numPlot=floor(numEvec/2)

colg <- read.table(FileCodePCA,header=T,stringsAsFactors=F,comment.char="@",sep="\t")
print("colg read")

pc<- read.table(evecFile,header=FALSE,stringsAsFactors=F)
print("evec read")
evals<-read.table(evalFile,header=FALSE,stringsAsFactors=F)
print("eval read")
colg<-colg[colg$Population %in% pc[,numEvec+2],]

pc$couleur=NA
pc$pointtype=NA
pc$cex=NA
pc$text=NA
evals<-evals/sum(evals)
evals<-round(evals*100000)/1000
#print(pc)

if(! "cex" %in% names(colg)){
    colg$cex=1
}
for (i in 1:length(pc$V2)) {
	print(pc[i,numEvec+2])
	pc$couleur[i] = colg$Color[colg$Population==pc[i,numEvec+2]]
	pc$pointtype[i] = colg$Point[colg$Population==pc[i,numEvec+2]]
    if("cex" %in% names(colg)){
        pc$cex[i] = colg$cex[colg$Population==pc[i,numEvec+2]]
    }
  if("Text" %in% names(colg)){
        pc$text[i] = colg$Text[colg$Population==pc[i,numEvec+2]]
    }
	#print(colg[colg$Population==pc[i,numEvec+2],])
}
colg$cex=as.numeric(colg$cex)
colg<-colg[order(colg$Region,colg$Population),]
minCex=min(pc$cex)
pc<-pc[order(as.numeric(pc$cex)),]
print(pc)
pdf(outFile)
par(mar=c(0, 0, 0, 0)+0.5)
plot(0,0,type="n",axes=FALSE,ann=FALSE)
colg<-colg[ order(colg$cex,colg$Region),]
tail(colg)
legend("topright",
	       legend=paste(colg$Text[! is.na(colg$Text)]," : ",colg$Population[! is.na(colg$Text)]," (",colg$Region[! is.na(colg$Text)],")",sep=""),
	       text.col=as.character(colg$Color[! is.na(colg$Text)]),
	       #pch=colg$Point[! is.na(colg$Text)],
	       cex=0.4,ncol=2
	     )
legend("topleft",
	       legend=paste(colg$Population[ is.na(colg$Text)]," (",colg$Region[ is.na(colg$Text)],")",sep=""),
	       text.col=as.character(colg$Color[ is.na(colg$Text)]),
	       col=ifelse(colg$Point[ is.na(colg$Text)]>=21,1,colg$Color[ is.na(colg$Text)]),
         pt.bg=as.character(colg$Color[ is.na(colg$Text)]),
	       pch=colg$Point[ is.na(colg$Text)],
	       cex=0.4,ncol=2
	)
par(mar=c(5, 4, 4, 2)+0.1)
pc<-pc[ order(as.numeric(pc$cex)),]
pc[ pc$V18=="CLM",]
for(iter in c(1:numPlot)){
	col=iter*2+1
    
	plot(pc[,col-1],pc[,col],"n",
		     col=ifelse(pc$pointtype<21,as.character(pc$couleur),1),
	       bg=as.character(pc$couleur),
		     pch=ifelse(is.na(pc$text),pc$pointtype,NA),
		     cex=pc$cex,
		     xlab=paste("PC",iter*2-1,sep=""),
		     ylab=paste("PC",iter*2,sep=""),
		     main=paste("PC",iter*2," (",evals[iter*2,1],"%)\nVS\nPC",iter*2-1," (",evals[iter*2-1,1],"%)",sep=""),
		     xlim=c(min(pc[,col-1]),max(pc[,col-1])),
		     ylim=c(min(pc[,col]),max(pc[,col])))
	points(pc[,col-1],pc[,col],
               col=ifelse(pc$pointtype<21,as.character(pc$couleur),1),
               bg=as.character(pc$couleur),
               pch=ifelse(is.na(pc$text),pc$pointtype,NA),
	       cex=pc$cex
        )
        text(y=pc[,col],x=pc[,col-1],
                     col=as.character(pc$couleur),
                      labels=pc$text,
                     cex=pc$cex)
		
      
}
print(warnings())
dev.off()

