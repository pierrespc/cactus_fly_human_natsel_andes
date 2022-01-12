#!/bin/Rscript

flipAllele<-function(s){
  return(ifelse(s=="A","T",ifelse(s=="T","A",ifelse(s=="C","G",ifelse(s=="G","C",ifelse(s=="0","0",stop(s)))))))
}

flipBim<-function(vec){
  A1toFlip=vec[1]
  A2toFlip=vec[2]
  A1Ref=vec[3]
  A2Ref=vec[4]
  ### ambiguous strand
  if(A1toFlip == flipAllele(A2toFlip) | A1Ref == flipAllele(A2Ref)){
    return(c(NA,NA))
  }
    
  ###alleles match
  if((A1toFlip %in% c(A1Ref,A2Ref,"0") & A2toFlip %in% c(A1Ref,A2Ref)) | (A1Ref %in% c(A1toFlip,A2toFlip,"0") & A2Ref %in% c(A1toFlip,A2toFlip))){
    return(c(A1toFlip,A2toFlip))
  }
  if((flipAllele(A1toFlip) %in% c(A1Ref,A2Ref,"0") & flipAllele(A2toFlip) %in% c(A1Ref,A2Ref)) | (flipAllele(A1Ref) %in% c(A1toFlip,A2toFlip,"0") & flipAllele(A2Ref) %in% c(A1toFlip,A2toFlip))){
      return(c(flipAllele(A1toFlip),flipAllele(A2toFlip)))
  }

  return(c(NA,NA))
  
}

params=commandArgs(trailingOnly=T)
if(length(params)!=4){
	stop("call with <prefFileToMerge> <prefFileStart> <prefOut> <plinkPath>")
}
prefFileToMerge<-params[1]
prefFileStart<-params[2]
prefOut<-params[3]
plinkPath=params[4]

bimStart<-read.table(paste(prefFileStart,".bim",sep=""),stringsAsFactors = F,header=F)
bimToMerge<-read.table(paste(prefFileToMerge,".bim",sep=""),stringsAsFactors = F,header=F)
bimToMerge<-unique(bimToMerge)

bimToMerge<-bimToMerge[bimToMerge$V5 %in% c("A","T","G","C","0") & bimToMerge$V6 %in% c("A","T","G","C","0") , ]
###merge by ID
merged1<-merge(bimToMerge,bimStart,by="V2")
if(dim(merged1)[1]>0){
  try<-t(apply(merged1[,c(5,6,10,11)],1,flipBim))
  merged1$NewA1=try[,1]
  merged1$NewA2=try[,2]
  #out<-merged1[! is.na(merged1$NewA1),c(2,1,8,4,12,13)]
  dict<-merged1[ ! is.na(merged1$NewA1),c(2,1,3,4,5,6,7,1,8,9,10,11,12,13)]
}else{
  #out<-data.frame(matrix(NA,0,6))
  dict<-data.frame(matrix(NA,0,14))
}
diffPos=merged1[ merged1$V4.x!=merged1$V4.y,]
diffPos<-diffPos[,c(2,1,3,4,5,6,7,1,8,9,10,11)]
names(diffPos)<-paste(rep(c("chr","snpID","cM","pos","A1","A2"),2),rep(c("toMerge","Ref"),each=6),sep="_")
diffPos$diffPos<-abs(diffPos$pos_toMerge - diffPos$pos_Ref)

#names(out)<-c("chr","snpID","cM","pos","A1","A2")
names(dict)<-c(paste(rep(c("chr","snpID","cM","pos","A1","A2"),2),rep(c("toMerge","Ref"),each=6),sep="_"),"NewA1","NewA2")
pbStrand1<-merged1[ is.na(merged1$NewA1),c(2,1,3,4,5,6,7,1,8,9,10,11)]
names(pbStrand1)<-paste(rep(c("chr","snpID","cM","pos","A1","A2"),2),rep(c("toMerge","Ref"),each=6),sep="_")

bimToMerge<-bimToMerge[ ! bimToMerge$V2 %in% merged1$V2,]
bimStart<-bimStart[ ! bimStart$V2 %in% merged1$V2,]

###merge by pos

merged2<-merge(bimToMerge,bimStart,by=c("V1","V4"))
if(dim(merged2)[1]>0){
  try<-t(apply(merged2[,c(5,6,9,10)],1,flipBim))
  merged2$NewA1=try[,1]
  merged2$NewA2=try[,2]
  #tmp<-merged2[! is.na(merged2$NewA1),c(1,7,8,2,11,12)]
  #names(tmp)<-c("chr","snpID","cM","pos","A1","A2")
  #out<-rbind(out,tmp)
  #remove(tmp)
  
  tmp2<-merged2[! is.na(merged2$NewA1),c(1,3,4,2,5,6,1,7,8,2,9,10,11,12)]
  names(tmp2)<-c(paste(rep(c("chr","snpID","cM","pos","A1","A2"),2),rep(c("toMerge","Ref"),each=6),sep="_"),"NewA1","NewA2")
  dict<-rbind(dict,tmp2)
  remove(tmp2)
  
}
diffID=merged2[ merged2$V2.x!=merged2$V2.y,]
bimToMerge<-bimToMerge[ ! bimToMerge$V2 %in% merged2$V2.x,]


pbStrand2<-merged2[ is.na(merged2$NewA1),c(1,3,4,2,5,6,1,7,8,2,9,10)]
names(pbStrand2)<-paste(rep(c("chr","snpID","cM","pos","A1","A2"),2),rep(c("toMerge","Ref"),each=6),sep="_")
pbStrand<-rbind(pbStrand1,pbStrand2)
remove(pbStrand1)
remove(pbStrand2)

if(dim(pbStrand)[1]>0){
  print(paste("some SNVs removed while flipping... see file ",prefOut,".pbStrand",sep=""))
  write.table(pbStrand,paste(prefOut,".pbStrand",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
}
if(dim(bimToMerge)[1]>0){
  print(paste(" SNVs not found by ID nor position... see file ",prefOut,".notFound",sep=""))
  write.table(bimToMerge,paste(prefOut,".notFound",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
}

if(sum(duplicated(dict$snpID_Ref))>0 | sum(duplicated(dict$snpID_toMerge))>0){
  message(paste("duplicated snpID...see file ",prefOut,".duplicatedSNP"))
  write.table(dict[ dict$snpID_Ref %in% dict$snpID_Ref[duplicated(dict$snpID_Ref)] | dict$snpID_toMerge %in% dict$snpID_toMerge[duplicated(dict$snpID_toMerge)],],paste(prefOut,".duplicatedSNP",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
}

if(dim(diffPos)[1]>0){
  print(paste("pb position see file ",prefOut,".diffPos",sep=""))
  write.table(diffPos,paste(prefOut,".diffPos",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
}

write.table(dict,paste(prefOut,".Changes",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
write.table(dict$snpID_toMerge,paste(prefOut,".ExtractList",sep=""),col.names=F,row.names=F,quote=F)
row.names(dict)<-dict$snpID_toMerge
system(paste(plinkPath," --bfile ",prefFileToMerge," --extract ",prefOut,".ExtractList"," --make-bed --out ",prefOut,sep=""))
if(file.exists(paste(prefOut,".duplicatedSNP",sep=""))){
  system(paste(plinkPath," --bfile ",prefOut," --exclude ",prefOut,".duplicatedSNP"," --make-bed --out ",prefOut,sep=""))
}
bim<-read.table(paste(prefOut,".bim",sep=""),header=F,stringsAsFactors = F)
out<-dict[ bim$V2,c(paste(c("chr","snpID","cM","pos"),rep("Ref",4),sep="_"),"NewA1","NewA2")]
if(dim(out)[1]!=dim(bim)[1]){
  system(paste("rm ",prefOut,".bim ",prefOut,".bed ",prefOut,".fam ",sep=""))
  stop("pb dimension byw dictionnary and bim")
}

write.table(out,paste(prefOut,".bim",sep=""),col.names=F,row.names=F,quote=F)
