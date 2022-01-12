#!/bin/Rscript


listCHROM<-c(22:15)
setwd("~/Documents/PostDoc/DetoxFly/Gnecchio/Reformate/Outputs/fineStructureOnlyNat/fineStructure/")


for(chr in listCHROM){
  print(chr)
  chuck<-read.csv(paste("stage2.chr",chr,".regionchunkcounts.out",sep=""),stringsAsFactors = F,header=T,sep=" ")
  row.names(chuck)<-chuck$Recipient
  
  if(chr==listCHROM[1]){
    numReg<-chuck$num.regions
    names(numReg)<-chuck$Recipient
    chuck<-
    out<-chuck[,-c(1,2)]
    listINDS<-row.names(chuck)
    
  }else{
    chuck<-chuck[listINDS,c("Recipient","num.regions",listINDS)]
    numReg<-numReg+chuck$num.regions
    names(numReg)<-chuck$Recipient
    out<-out+chuck[,-c(1,2)]
  }
}

out$Recipient<-row.names(out)
out$num.regions<-numReg
out<-out[,c("Recipient","num.regions",listINDS)]
write.table(out,"stage2.Combined",col.names = T,row.names = F,sep=" ",quote=F)


