#!/bin/Rscript


K=8
setwd("~/Documents/PostDoc/DetoxFly//LAST/MaoMoreno/Reformate/Outputs/")

adm<-read.table(paste("1KG/Admixture/BestRUNperK/6pops.KinshipFilter.1KG.MAF0.0001.GENO0.02.MIND0.05.pruned.K",K,".AncestryComponentByIndividual.txt",sep=""),stringsAsFactors = F,header=T,comment.char = "@")


andean<-apply(adm[ adm$Population %in% c("QUECHUAN","AYMARAN"),paste("V",c(1:K),sep="")],2,mean)
andean<-which(andean==max(andean))

mex<-apply(adm[ adm$Population %in% c("MAYAN","NAHUAN"),paste("V",c(1:K),sep="")],2,mean)
mex<-which(mex==max(mex))

bar<-apply(adm[ adm$Population %in% c("BAR"),paste("V",c(1:K),sep="")],2,mean)
bar<-which(bar==max(bar))

yuk<-apply(adm[ adm$Population %in% c("YUK"),paste("V",c(1:K),sep="")],2,mean)
yuk<-which(yuk==max(yuk))

clm<-apply(adm[ adm$Population %in% c("CLM"),paste("V",c(1:K),sep="")],2,mean)
clm<-which(clm==max(clm))

mxl<-apply(adm[ adm$Population %in% c("MXL"),paste("V",c(1:K),sep="")],2,mean)
mxl<-which(mxl==max(mxl))

pur<-apply(adm[ adm$Population %in% c("PUR"),paste("V",c(1:K),sep="")],2,mean)
pur<-which(pur==max(pur))

#adm$NAT<-apply(adm[,paste("V",unique(c(clm,mex,andean,bar,mxl,pur,yuk)),sep="")],1,sum)
adm$NAT<-apply(adm[,paste("V",unique(c(mex,andean,bar,yuk)),sep="")],1,sum)

print(summary(adm$NAT[ adm$Population == "MAYAN"]))
###NO NAT
esn<-apply(adm[ adm$Population %in% c("ESN"),paste("V",c(1:K),sep="")],2,mean)
esn<-which(esn==max(esn))

lwk<-apply(adm[ adm$Population %in% c("LWK"),paste("V",c(1:K),sep="")],2,mean)
lwk<-which(lwk==max(lwk))

msl<-apply(adm[ adm$Population %in% c("MSL"),paste("V",c(1:K),sep="")],2,mean)
msl<-which(msl==max(msl))

gwd<-apply(adm[ adm$Population %in% c("GWD"),paste("V",c(1:K),sep="")],2,mean)
gwd<-which(gwd==max(gwd))


tsi<-apply(adm[ adm$Population %in% c("TSI"),paste("V",c(1:K),sep="")],2,mean)
tsi<-which(tsi==max(tsi))

fin<-apply(adm[ adm$Population %in% c("FIN"),paste("V",c(1:K),sep="")],2,mean)
fin<-which(fin==max(fin))


adm$NONAT<-apply(adm[,paste("V",unique(c(msl,lwk,esn,gwd,tsi,fin)),sep="")],1,sum)
#last<-c(1:K)[ ! c(1:K) %in% c(msl,lwk,esn,gwd,tsi,fin,mxl,yuk,bar,pur,clm,andean,mex)]
#adm$NAT<-apply(adm[,paste("V",unique(c(mxl,yuk,bar,last,andean,mex,pur,clm)),sep="")],1,sum)

adm<-adm[ order(as.numeric(adm$NAT),decreasing = T),]
keep<-adm[adm$NAT>0.95 & adm$Population %in% c("YUK","BAR","NAHUAN","MAYAN","AYMARAN","QUECHUAN"),c("Population","Ind")]
keep$Code<-as.numeric(as.factor(keep$Population))
keep$CodP<-keep$Population


write.table(keep,"6pops.KinshipFilter.AdmFilter.KEEP",quote=F,sep="\t",col.names = F,row.names = F)
system("~/src/plink1.9/plink --bfile 6pops.KinshipFilter --keep 6pops.KinshipFilter.AdmFilter.KEEP --make-bed --out 6pops.KinshipFilter.AdmFilter")
system("~/src/plink1.9/plink --bfile 6pops.KinshipFilter.AdmFilter --geno 0.02 --make-bed --out 6pops.KinshipFilter.AdmFilter")
system("~/src/plink1.9/plink --bfile 6pops.KinshipFilter.AdmFilter --mind 0.05 --maf 0.0001 --make-bed --out 6pops.KinshipFilter.AdmFilter.MAF0.0001.GENO0.02.MIND0.05")

print(table(keep$Population))
print("went well")

