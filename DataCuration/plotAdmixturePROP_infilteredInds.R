#!/bin/Rscript
setwd("~/Documents/OtherCollaborations/DetoxFly/LAST/CheckStructure/Merged_noUros/")
famFINAL<-read.table("Merged_noUros.fam",stringsAsFactors = F,header=F)


COLORS<-c("Bari"="orangered",
           "Yukpa"="#FD7F23",
          "Wichi"="deeppink",
          "Yanesha"="#9751A1",
          "NoNatAm"="grey",
          "Andes"="#A7D65D")
admGnecchio<-read.table("../../Gnecchio/Reformate/Outputs/1KG/Admixture/BestRUNperK/SouthAmerica_dataset_QCed.KinshipFilter.1KG.pruned.K8.AncestryComponentByIndividual.txt",stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
admGnecchio<-admGnecchio[ row.names(admGnecchio) %in% famFINAL$V2,]
admGnecchioPOP<-read.table("../../Gnecchio/Reformate/Outputs/1KG/Admixture/BestRUNperK/SouthAmerica_dataset_QCed.KinshipFilter.1KG.pruned.K8.MeanByGroup.txt",stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
yanesha=which(admGnecchioPOP["Yanesha_IntermediateSelva",]==max(admGnecchioPOP["Yanesha_IntermediateSelva",]))
wichi<-which(admGnecchioPOP["Wichi",]==max(admGnecchioPOP["Wichi",]))
andes<-which(admGnecchioPOP["Andes",]==max(admGnecchioPOP["Andes",]))

admGnecchio$Yanesha=admGnecchio[,paste("V",yanesha,sep="")]
admGnecchio$Wichi=admGnecchio[,paste("V",wichi,sep="")]
admGnecchio$Andes=admGnecchio[,paste("V",andes,sep="")]
admGnecchio$NoNatAm=apply(admGnecchio[,paste("V",c(1:8)[-c(yanesha,wichi,andes)],sep="")],1,sum)


svg
plot(0,0,"n",xlim=c(1,2+nrow(admGnecchio)),ylim=c(0,1),axes=F,ann=F)
x=-1
for(pop in c("Yanesha","Wichi","Andes")){
  x=x+1
  if(pop=="Yanesha" | pop=="Wichi"){
    tmp=admGnecchio[ grepl(pop,admGnecchio$Population),]
  }else{
    if(pop=="Andes"){
      tmp=admGnecchio[ admGnecchio$Region==pop,]
    }else{
      stop("HU")
    }
  }
  
  tmp<-tmp[order(tmp$Yanesha,tmp$Wichi,tmp$Andes),]
    for(ind in row.names(tmp)){
      x=x+1
      ybottom=0
      for(ele in c("Yanesha","Wichi","Andes","NoNatAm")){
        ytop=ybottom+tmp[ind,ele]
        rect(xleft = x,xright = x+1,ybottom = ybottom,ytop=ytop,border=NA,
             col = COLORS[ele])
        ybottom=ytop
      }
    }
  }    
}

admMao<-read.table("../../MaoMoreno//Reformate/Outputs/1KG/Admixture/BestRUNperK/6pops.KinshipFilter.1KG.MAF0.0001.GENO0.02.MIND0.05.pruned.K8.AncestryComponentByIndividual.txt",stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
admMao<-admMao[ row.names(admMao) %in% famFINAL$V2,]
admMaoPOP<-read.table("../../MaoMoreno//Reformate/Outputs/1KG/Admixture/BestRUNperK/6pops.KinshipFilter.1KG.MAF0.0001.GENO0.02.MIND0.05.pruned.K8.MeanByGroup.txt",stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
bari=which(admMaoPOP["BAR",]==max(admMaoPOP["BAR",]))
yukpa<-which(admMaoPOP["YUK",]==max(admMaoPOP["YUK",]))
andes<-which(admMaoPOP["Andes",]==max(admMaoPOP["Andes",]))

admMao$Yukpa=admMao[,paste("V",yukpa,sep="")]
admMao$Bari=admMao[,paste("V",bari,sep="")]
admMao$Andes=admMao[,paste("V",andes,sep="")]
admMao$NoNatAm=apply(admMao[,paste("V",c(1:8)[-c(yukpa,bari,andes)],sep="")],1,sum)


plot(0,0,"n",xlim=c(1,2+nrow(admMao)),ylim=c(0,1),axes=F,ann=F)
x=-1
for(pop in c("YUK","BAR","Andes")){
  x=x+1
  if(pop=="YUK" | pop=="BAR"){
    tmp=admMao[ admMao$Population==pop,]
  }else{
    if(pop=="Andes"){
      tmp=admMao[ admMao$Region==pop,]
    }else{
      stop("HU")
    }
  }
  
  tmp<-tmp[order(tmp$Yukpa,tmp$Bari,tmp$Andes,tmp$NoNatAm),]
  for(ind in row.names(tmp)){
    x=x+1
    ybottom=0
    for(ele in c("Yukpa","Bari","Andes","NoNatAm")){
      ytop=ybottom+tmp[ind,ele]
      rect(xleft = x,xright = x+1,ybottom = ybottom,ytop=ytop,border=NA,
           col = COLORS[ele])
      ybottom=ytop
    }
  }
}    


