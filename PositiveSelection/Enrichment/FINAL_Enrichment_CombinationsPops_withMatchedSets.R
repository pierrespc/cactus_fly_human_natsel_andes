#!/bin/Rscript
require(plotrix)
require(stringr)
setwd("~/Documents/OtherCollaborations/DetoxFly/LAST/Permutations_Intersection2sets/MatchedGeneSets/")

####
checkNA<-function(vector){
    return(sum(is.na(vector))>0)
}


allSignif<-function(vector,th){
  return(sum(vector<th)==length(vector))
}

atleastOne<-function(vector){
  return(sum(vector)>0)
}

signifLevel=0.05
listPops<-c("Yanesha","Wichi","Bari","Yukpa")
listPopsCombination<-c("Bari","Yukpa","Wichi")
listOrtho<-read.table("../../../OrthoGenes.txt",stringsAsFactors = F,header=F)
names(listOrtho)<-c("HGNC","ENSID","chr","Start","End")


numGenes=1000
buf=10
outtable<-data.frame(matrix(NA,0,4))
names(outtable)<-c("N","Nsignif","Group.proportion","P.proportion")

for(pop in listPops){
  outtable<-rbind(outtable,outtable[1,])
  row.names(outtable)[nrow(outtable)]<-paste(pop,"mean",sep = "_")
  outtable<-rbind(outtable,outtable[1,])
  row.names(outtable)[nrow(outtable)]<-paste(pop,"median",sep = "_")
}

for(n in c(1:length(listPopsCombination))){
      outtable<-rbind(outtable,outtable[1,])
      row.names(outtable)[nrow(outtable)]<-paste("InAtLeast",n,"Comparisons_mean",sep="")
      outtable<-rbind(outtable,outtable[1,])
      row.names(outtable)[nrow(outtable)]<-paste("InAtLeast",n,"Comparisons_median",sep="")
}

proporSets<-list()
for(sumStat in c("mean","median")){
  print(paste(sumStat,buf,"KB"))
  print("reading and merging")
  #####FROM MAO MORENO DATASET
  bari<-read.table(paste("../../MaoMoreno/Fisher/Andean.refBari.BKGDGenes.",buf,"KB.Fisher.",sumStat,".1Groups.tsv",sep=""),
  stringsAsFactors = F,header=T)
  bari<-bari[ bari$Nsnps_LRT>0 & bari$Nsnps_absIHS>0,c("HGNC","ENSID","chr","Start","End","Nsnps_absIHS","Nsnps_LRT","Zf","P_Zf")]
  names(bari)<-c("HGNC","ENSID","chr","Start","End","Nihs_Moreno","N_Bari","Zf_Bari","P_Zf_Bari")
  yuk<-read.table(paste("../../MaoMoreno/Fisher/Andean.refYukpa.BKGDGenes.",buf,"KB.Fisher.",sumStat,".1Groups.tsv",sep=""),
                  stringsAsFactors = F,header=T)
  yuk<-yuk[ yuk$Nsnps_LRT>0 & yuk$Nsnps_absIHS>0,c("HGNC","ENSID","chr","Start","End","Nsnps_LRT","Zf","P_Zf")]
      names(yuk)<-c("HGNC","ENSID","chr","Start","End","N_Yukpa","Zf_Yukpa","P_Zf_Yukpa")

  out1<-merge(bari,yuk,by=c("HGNC","ENSID","chr","Start","End"),all=T)

  #####FROM GNECCHI_RUSCONE DATASET
  yanesha<-read.table(paste("../../Gnecchio_noUros//Fisher/Andean.refYanesha.BKGDGenes.",buf,"KB.Fisher.",sumStat,".1Groups.tsv",sep=""),
                  stringsAsFactors = F,header=T)
  yanesha<-yanesha[ yanesha$Nsnps_LRT>0 & yanesha$Nsnps_absIHS>0,c("HGNC","ENSID","chr","Start","End","Nsnps_absIHS","Nsnps_LRT","Zf","P_Zf")]
  names(yanesha)<-c("HGNC","ENSID","chr","Start","End","Nihs_Gnecchi","N_Yanesha","Zf_Yanesha","P_Zf_Yanesha")

  wichi<-read.table(paste("../../Gnecchio_noUros//Fisher/Andean.refWichi.BKGDGenes.",buf,"KB.Fisher.",sumStat,".1Groups.tsv",sep=""),
                stringsAsFactors = F,header=T)

  wichi<-wichi[ wichi$Nsnps_LRT>0 & wichi$Nsnps_absIHS>0,c("HGNC","ENSID","chr","Start","End","Nsnps_LRT","Zf","P_Zf")]
  names(wichi)<-c("HGNC","ENSID","chr","Start","End","N_Wichi","Zf_Wichi","P_Zf_Wichi")

  out2<-merge(yanesha,wichi,by=c("HGNC","ENSID","chr","Start","End"),all=T)

  ####all
  out<-merge(out1,out2,by=c("HGNC","ENSID","chr","Start","End"),all=T)

  remove(out1)
  remove(out2)
  remove(yuk)
  remove(bari)
  remove(yanesha)
  remove(wichi)
  ###now make average recomb rate...
  out$recomb=NA
  outfinal<-out[0,]
  for(chr in seq(1,22)){
    print(paste("average recomb rate chr",chr))
    recomb<-read.table(paste("/Volumes/MARIOLOLO/Data/1KG_Haplotypes/Phase3/genetic_map_chr",chr,"_combined_b37.txt",sep=""),
    stringsAsFactors = F,header=T)
    tmp<-out[ out$chr==chr,]
    for(i in c(1:nrow(tmp))){
      tmp$recomb[i]<-mean(recomb$COMBINED_rate.cM.Mb.[recomb$position>=tmp$Start[i] & recomb$position<=tmp$End[i]])
    }
    outfinal<-rbind(outfinal,tmp)
  }

      
  print(paste("making ",numGenes,"matched sets"))
  outfinal<-outfinal[ ! is.na(outfinal$recomb),]
  outfinal$ISNA<-apply(outfinal[ ,paste("Zf_",listPops,sep="")],1,checkNA)
  matchedSets<-list()
  for(i in c(1:nrow(listOrtho))){
    line=outfinal[ outfinal$ENSID==listOrtho$ENSID[i],]
    if(nrow(line)==0){
      print(listOrtho$HGNC[i])
      next
    }
    if(nrow(line)!=1){
      print(listOrtho$HGNC[i])
      stop(i)
    }
    outfinal$distRecomb<-rank(abs(outfinal$recomb-line$recomb),na.last = T,ties.method = "min")
    outfinal$distIHSmoreno<-rank(abs(outfinal$Nihs_Moreno-line$Nihs_Moreno),na.last = T,ties.method = "min")
    outfinal$distIHSgnecchi<-rank(abs(outfinal$Nihs_Gnecchi-line$Nihs_Gnecchi),na.last = T,ties.method = "min")
    outfinal$distIHS<-rank(abs(outfinal$Nihs_Moreno-line$Nihs_Moreno)+abs(outfinal$Nihs_Gnecchi-line$Nihs_Gnecchi),na.last = T,ties.method = "min")
    outfinal$distTot=outfinal$distRecomb+outfinal$distIHSmoreno+outfinal$distIHSgnecchi
    tmp1<-outfinal[ ! outfinal$ISNA & ! outfinal$ENSID %in% listOrtho$ENSID ,]
    #tmp1<-outfinal
    kept<-tmp1[0,]
    listDist<-sort(unique(tmp1$distIHS))
    dist=0
    while(nrow(kept)<numGenes){
      dist=dist+1
      tmp=tmp1[ tmp1$distIHS==listDist[dist],]
      kept<-rbind(kept,head(tmp[order(tmp$distRecomb),],numGenes-nrow(kept)))
    }
    kept<-kept[ sample(c(1:nrow(kept)),nrow(kept),replace = F),]
    ENSID=line$ENS
    for(set in c(1:numGenes)){
      matchedSets[[paste("Set",set,".P",sep="")]]<-rbind(matchedSets[[paste("Set",set,".P",sep="")]],cbind(ENSID,kept[set,paste("P_Zf_",listPops,sep="")]))
    }
  }
    
  print(paste("estimating enrichment for individual combinations with ",numGenes,"sets"))
  for(popops in listPops){
    testPval<-na.omit(outfinal[ outfinal$ENSID %in% listOrtho$ENSID,c("ENSID",paste("P_Zf_",popops,sep=""))])
    outtable[paste(popops,sumStat,sep="_"),"N"]<-nrow(testPval)
    if(length(popops)==1){
      outtable[paste(popops,sumStat,sep="_"),"Nsignif"]<-sum(testPval[,paste("P_Zf_",popops,sep="")]<signifLevel)
    }else{
      outtable[paste(popops,sumStat,sep="_"),"Nsignif"]<-sum(apply(testPval[,paste("P_Zf_",popops,sep="")],1,allSignif,th=signifLevel))
    }
    outtable[paste(popops,sumStat,sep="_"),"Group.proportion"]<-outtable[paste(popops,sumStat,sep="_"),"Nsignif"]/outtable[paste(popops,sumStat,sep="_"),"N"]
    countProp=0
    listPropSets=c()
    for(set in c(1:numGenes)){
      match<-matchedSets[[paste("Set",set,".P",sep="")]]
      match<-match[ match$ENSID %in% testPval$ENSID,paste("P_Zf_",popops,sep="")]
      prop=sum(match<signifLevel)/length(match)
      if(prop>outtable[paste(popops,sumStat,sep="_"),"Group.proportion"]){
        countProp=countProp+1
      }
      listPropSets=c(listPropSets,prop)
    }
    outtable[paste(popops,sumStat,sep="_"),"P.proportion"]<-countProp/numGenes
    proporSets[[paste(popops,sumStat,sep="_")]]<-listPropSets
  }
  print(paste("test the union of combinations of N pops,with ",numGenes," sets"))
  for(n in c(1:length(listPopsCombination))){
    orthoTable<-na.omit(outfinal[ outfinal$ENSID %in% listOrtho$ENSID,paste("P_Zf_",listPopsCombination,sep="")])
    listToTest<-apply(combn(listPopsCombination,n),2,paste,collapse="_")
    comb=0
    for(popToTest in listToTest){
      comb=comb+1
      orthoTable<-cbind(orthoTable,NA)
      names(orthoTable)[length(orthoTable)]<-paste("Comb",comb,sep="")
      popops=strsplit(popToTest,split="_")[[1]]
      if(length(popops)==1){
        orthoTable[,paste("Comb",comb,sep="")]<-orthoTable[,paste("P_Zf_",popops,sep="")]<signifLevel
      }else{
        orthoTable[,paste("Comb",comb,sep="")]<-apply(orthoTable[,paste("P_Zf_",popops,sep="")],1,allSignif, th=signifLevel)
      }
    }
    if(comb>1){
      orthoTable$AtLeastOne<-apply(orthoTable[,paste("Comb",seq(1,comb),sep="")],1,atleastOne)
    }else{
      orthoTable$AtLeastOne=orthoTable[,paste("Comb",comb,sep="")]
    }
        
    outtable[paste("InAtLeast",n,"Comparisons_",sumStat,sep=""),"N"]<-nrow(orthoTable)
    outtable[paste("InAtLeast",n,"Comparisons_",sumStat,sep=""),"Nsignif"]<-sum(orthoTable$AtLeastOne)
    outtable[paste("InAtLeast",n,"Comparisons_",sumStat,sep=""),"Group.proportion"]<-
    outtable[paste("InAtLeast",n,"Comparisons_",sumStat,sep=""),"Nsignif"]/outtable[paste("InAtLeast",n,"Comparisons_",sumStat,sep=""),"N"]
      
    ###now get the proprtion of gene set with higher proportion of signif genes
    countProp=0
    listPropSets=c()
    for(set in c(1:numGenes)){
      match<-matchedSets[[paste("Set",set,".P",sep="")]]
      match<-match[ match$ENSID %in% testPval$ENSID,paste("P_Zf_",listPops,sep="")]
      comb=0
      for(popToTest in listToTest){
        comb=comb+1
        match<-cbind(match,NA)
        names(match)[length(match)]<-paste("Comb",comb,sep="")
        popops=strsplit(popToTest,split="_")[[1]]
        if(length(popops)==1){
          match[,paste("Comb",comb,sep="")]<-match[,paste("P_Zf_",popops,sep="")]<signifLevel
        }else{
          match[,paste("Comb",comb,sep="")]<-apply(match[,paste("P_Zf_",popops,sep="")],1,allSignif, th=signifLevel)
        }
      }
      if(comb>1){
        match$AtLeastOne<-apply(match[,paste("Comb",seq(1,comb),sep="")],1,atleastOne)
      }else{
        match$AtLeastOne<-match[,paste("Comb",comb,sep="")]
      }
      if(sum(match$AtLeastOne/nrow(match)) > outtable[paste("InAtLeast",n,"Comparisons_",sumStat,sep=""),"Group.proportion"]){
        countProp=countProp+1
      }
      listPropSets=c(listPropSets,sum(match$AtLeastOne/nrow(match)))
    }
    outtable[paste("InAtLeast",n,"Comparisons_",sumStat,sep=""),"P.proportion"]<-countProp/numGenes
    proporSets[[paste("InAtLeast",n,"Comparisons_",sumStat,sep="")]]<-listPropSets
        
  }
}
    
write.table(outtable,paste("FINAL_Enrichment_PopCombinations_",buf,"KB_numSet",numGenes,".tsv",sep=""),col.names = T,row.names = T,sep="\t",quote=F)

max=max(c(outtable$Group.proportion,unlist(proporSets)))*1.05
pdf(paste("FINAL_Enrichment_PopCombinations_",buf,"KB_numSet",numGenes,".pdf",sep=""),width=15,height=5)
plot(0,0,"n",axes=F,ann=F,xlim=c(1,nrow(outtable)),ylim=c(0,max))
abline(h=0)
for(xPlot in c(1:nrow(outtable))){
  c("Wichi"="deeppink","Yanesha"="#9751A1")
  c("Yukpa"="#FD7F23","Bari"="orangered")
  col=ifelse(grepl("Comparisons",row.names(outtable)[xPlot]),"grey",
             ifelse(grepl("Bari",row.names(outtable)[xPlot]),"orangered",
                    ifelse(grepl("Yukpa",row.names(outtable)[xPlot]),"#FD7F23",
                           ifelse(grepl("Wichi",row.names(outtable)[xPlot]),"deeppink",
                                  ifelse(grepl("Yanesha",row.names(outtable)[xPlot]),"#9751A1",stop())))))
  violin_plot(proporSets[[row.names(outtable)[xPlot]]],at=xPlot,add=T,axes=F,col=col,main=NA,violin_width = 0.8)
  segments(x0=xPlot-0.3,x1=xPlot+0.3,
           y0=outtable[ row.names(outtable)[xPlot], "Group.proportion" ],
           y1=outtable[ row.names(outtable)[xPlot], "Group.proportion" ],col="black",lwd=4)
  #ax=str_replace(row.names(outtable)[xPlot],"_mean","\n(mean)")
  #ax=str_replace(ax,"_median","\n(median)")
  #ax=str_replace_all(ax,"_","/")
  #ax=str_replace(ax,"InAtLeast","At least\n")
  #ax=str_replace(ax,"Comparisons"," Comparisons")
  #axis(1,at=xPlot,labels = ax,cex.axis=0.8,las=1,tick = F)
      
  P=outtable[ row.names(outtable)[xPlot], "P.proportion" ]
  signif=ifelse(P<0.01,"**",ifelse(P<0.05,"*",ifelse(P<0.1,"-",NA)))
  text(x=xPlot,y=max,labels = signif,cex=2)
}
axis(2)

change<-function(string,split,pos){
  return(strsplit(string,split=split)[[1]][pos])
}
sumStat<-sapply(row.names(outtable),change,split="_",pos=2)
Comparison<-sapply(row.names(outtable),change,split="_",pos=1)
Comparison=str_replace(Comparison,"InAtLeast","At least ")
Comparison=str_replace(Comparison,"Comparisons"," Comparisons")
Comparison=str_replace(Comparison,"1 Comparisons","1 Comparison")
Comparison=str_replace(Comparison,"At least 3 Comparisons","All 3 Comparisons")
Comparison[ ! grepl("Comparison",Comparison)]<-paste(Comparison[ ! grepl("Comparison",Comparison)]," Reference",sep="")
axis(1,at=c(1:nrow(outtable)),labels = sumStat,cex.axis=1,las=1,tick = F,line=0)
axis(1,at=seq(1.5,nrow(outtable)-0.5,length.out=nrow(outtable)/2),labels = unique(Comparison),cex.axis=1,las=1,tick = F,line=-1)

compa<-which(grepl("Compar",row.names(outtable)))
axis(1,at=mean(compa),labels = "Not considering comparisons with Yanesha as reference population",cex.axis=1,las=1,tick = F,line=2)
legend("topleft",pch=NA,legend = c("-  : P-value < 0.1","*   : P-value < 0.05","** : P-value < 0.01"))
dev.off()
