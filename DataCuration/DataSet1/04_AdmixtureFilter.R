#!/bin/Rscript


K=8
setwd("~/Documents/PostDoc/DetoxFly/LAST/Gnecchio/Reformate/Outputs/")

adm<-read.table(paste("1KG/Admixture/BestRUNperK/SouthAmerica_dataset_QCed.KinshipFilter.1KG.pruned.K",K,".AncestryComponentByIndividual.txt",sep=""),stringsAsFactors = F,header=T,comment.char = "@")


andean<-apply(adm[ adm$Population %in% c("Titicaca_Aymara"),paste("V",c(1:K),sep="")],2,mean)
andean<-which(andean==max(andean))

ET<-apply(adm[ adm$Population %in% c("Ashaninka"),paste("V",c(1:K),sep="")],2,mean)
ET<-which(ET==max(ET))


GPC<-apply(adm[ adm$Population %in% c("Wichi"),paste("V",c(1:K),sep="")],2,mean)
GPC<-which(GPC==max(GPC))


#adm$NAT<-apply(adm[,paste("V",unique(c(clm,mex,andean,bar,mxl,pur,yuk)),sep="")],1,sum)
adm$NAT<-apply(adm[,paste("V",unique(c(GPC,andean,ET)),sep="")],1,sum)


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

adm<-adm[ order(as.numeric(adm$NAT),decreasing = T),]
#stop(paste(adm$NAT[ adm$Population=="Tzotzil"],collapse=" "))
keep<-adm[adm$NAT>0.99 & adm$Region %in% c("Northern-Amerind","Andes","Ge-Pano-Carib","Equatorial-Tucanoan"),c("Population","Ind")]
keep$Code<-as.numeric(as.factor(keep$Population))
keep$CodP<-keep$Population

write.table(keep,"SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.KEEP",quote=F,sep="\t",col.names = F,row.names = F)
system("~/src/plink1.9/plink --bfile SouthAmerica_dataset_QCed.KinshipFilter --keep SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.KEEP --make-bed --out SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter")
system("~/src/plink1.9/plink --bfile SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter --geno 0.02 --make-bed --out SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter")
system("~/src/plink1.9/plink --bfile SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter --mind 0.05 --maf 0.0001 --make-bed --out SouthAmerica_dataset_QCed.KinshipFilter.AdmFilter.MAF0.0001.GENO0.02.MIND0.05")

print("went well")

