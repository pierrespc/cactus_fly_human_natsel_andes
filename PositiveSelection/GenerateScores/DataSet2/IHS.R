#!/bin/Rscript

require(rehh)
params<-commandArgs(trailingOnly = T)
params<-c("/Users/pierrespc/Documents/PostDoc/DetoxFly/LAST/MaoMoreno/Reformate/Outputs/IHS/Andean",
          0.01,0.02,0.05,
          "/Users/pierrespc/Documents/PostDoc/DetoxFly/LAST/MaoMoreno/Reformate/Outputs/6pops.KinshipFilter.AdmFilter.CPFS.ancestral",
          "/Volumes/MARIOLOLO/Data/1KG_Haplotypes/Phase3/",
          "/Users/pierrespc/Documents/PostDoc/DetoxFly/LAST/MaoMoreno/IHS")

if(length(params)!=7){
  stop("<shapeITpref>  <maf> <geno> <mind> <ancFile> <folderGenetMap> <DirOut>")
}

####
flipAllele<-function(s){
  return(ifelse(s=="A","T",ifelse(s=="T","A",ifelse(s=="C","G",ifelse(s=="G","C",ifelse(s=="0","0",stop(s)))))))
}
###

###recode Hap
recodeHap<-function(vec_withBool){
  ##vec with bool is a vector of alleles encoded 0/1 with two first elements A1 and A2
  
  return(ifelse(vec_withBool[-c(1:2)]==0,vec_withBool[1],
                    ifelse(vec_withBool[-c(1:2)]==1,vec_withBool[2],-9)
                    )
            )
}

shapeITpref<-params[1]
maf<-as.numeric(params[2])
geno<-as.numeric(params[3])
mind<-as.numeric(params[4])
ancPref<-params[5]
folderMap=params[6]
DirOut<-params[7]


system(paste("mkdir ",DirOut,sep=""))
system(paste("mkdir ",DirOut,"/IntermediateFiles/",sep=""))
samples<-read.table(paste(shapeITpref,".sample",sep=""),stringsAsFactors = F,header=F,skip=2)
listPop<-unique(samples$V1)

haps<-read.table(paste(shapeITpref,".haps",sep=""),stringsAsFactors = F,header=F,skip=2)
row.names(haps)<-haps$V2


ancestral<-read.table(ancPref,stringsAsFactors=F,header=F)
names(ancestral)<-c("chr","pos","Anc","Der")

if(dim(haps)[2]!=dim(samples)[1]*2+5){
  stop("pb dimensions sample/haps")
}
map<-haps[,c(1:5)]
names(map)<-c("chr","id","pos","A1","A2")

if(sum(duplicated(map$id))>0){
  stop("duplicated snp IDs")
}
haps<-haps[,-c(1:5)]

outMAP<-c()
outHAP<-c()

mapINFO<-c()
for(chr in unique(map$chr)){
  mapCHR<-map[ map$chr==chr,]
  row.names(mapCHR)<-mapCHR$id
  print(paste("reading anc file for chr",chr,sep=""))
  #ancCHR<-read.table(paste(ancPref,chr,".ancestral",sep=""),stringsAsFactors = F,header=F)
  #names(ancCHR)<-c("chr","pos","Anc","Der")
  ancCHR<-ancestral[ ancestral$chr==chr,]
  merged<-merge(mapCHR,ancCHR,by=c("chr","pos"),all.x=T)
  row.names(merged)<-merged$id
  print(paste("actualize map allele annot with checks on  Anc/Der allele info for chr",chr,sep=""))
  #first remove no SNP variants
  noSNP<-row.names(merged)[ ! (merged$A1 %in% c("0","A","T","G","C") & merged$A2 %in% c("0","A","T","G","C"))]
  #first do not consider ambigous genotyps
  #ambiguous<-row.names(merged)[ merged$A1 == sapply(merged$A2,flipAllele,USE.NAMES = F)]
  #merged$Anc[ambiguous]<-NA
  #merged$Der[ambiguous]<-NA
  #snp with no Ancestral info f
  noAncInfo<-row.names(merged)[  is.na(merged$Anc) | is.na(merged$Der) ]
  merged<-merged[ !row.names(merged) %in% noAncInfo,]
  ##check were we have to flip ancestral info
  #flippedAnc<-sapply(merged$Anc,flipAllele,USE.NAMES =F)
  #flippedDer<-sapply(merged$Der,flipAllele,USE.NAMES =F)
  #toFlip<-(flippedAnc ==merged$A1 & flippedDer ==merged$A2) |  
  #                             (flippedAnc ==merged$A2 & flippedDer ==merged$A1)
  #merged$Anc[toFlip]<-flippedAnc[toFlip]
  #merged$Der[toFlip]<-flippedDer[toFlip]
  
  ###check where ancestral info doesn't coindice A1/A2
  pb<-(merged$A1!=merged$Anc & merged$A2!=merged$Anc) | (merged$A1!=merged$Der & merged$A2!=merged$Der)
  merged<-merged[! pb,]
  
  ##check all is good
  OK<-(merged$A1==merged$Anc & merged$A2==merged$Der) | (merged$A1==merged$Der & merged$A2==merged$Anc)
  if(sum(OK)!=dim(merged)[1]){
    stop("pb recoding ancestral alleles")
  }
  
  print(paste("reorder the map as read for chr",chr,sep=""))
  mapCHR<-mapCHR[ row.names(mapCHR)%in% row.names(merged),]
  merged<-merged[row.names(mapCHR),]
  
  print(paste("recode alleles in hap chr",chr,sep=""))
  hapCHR<-haps[ row.names(haps)%in% row.names(merged),]
  hapCHRrecoded<-t(apply(cbind(merged$A1,merged$A2,hapCHR),1,recodeHap))
  
  ##wirting rehh formated files
  for(pop in unique(listPop)){
    whichHAPcol<-which(samples$V1 ==pop)
    whichHAPcol<-rep(whichHAPcol,each=2)*2+c(-1,0)
    hapPOP<-data.frame(hapCHRrecoded[,whichHAPcol],stringsAsFactors = F)
    write.table(hapPOP,paste(DirOut,"/IntermediateFiles/rehh.",pop,".chr",chr,".haps",sep=""),row.names = F,col.names = F,sep=" ",quote=F)
  }
  
  print("make genetic map instead of physical map") 
  genetMap<-read.table(paste(folderMap,"/genetic_map_chr",chr,"_combined_b37.txt",sep=""),stringsAsFactors = F,header=T)
  merged$cM<-approx(x=genetMap$position,y=genetMap$Genetic_Map.cM.,xout=merged$pos,rule=2)$y
  ###avoid duplicated
  for(dup in unique(merged$cM[duplicated(merged$cM)])){
    merged$cM[merged$cM == dup]<-merged$cM[merged$cM == dup]+c(1:sum(merged$cM == dup))*1e-8
  }
  write.table(merged[,c("id","chr","cM","Anc","Der")],paste(DirOut,"/IntermediateFiles/rehh.chr",chr,".map",sep=""),row.names = F,col.names = F,sep=" ",quote=F)
  mapINFO<-rbind(mapINFO,merged[,c("id","chr","pos")])
}

mapINFO<-data.frame(mapINFO,stringsAsFactors = F)

row.names(mapINFO)<-mapINFO$id

for(pop in unique(listPop)){
  if(! pop %in% c("Andean")){
    next
  }
  scanLIST<-c()
  for(chr in unique(map$chr)){
    print(paste("run ehh for ",pop,"chr",chr,sep=""))
    haplo<-data2haplohh(hap_file=paste(DirOut,"/IntermediateFiles/rehh.",pop,".chr",chr,".haps",sep=""),
                          map_file=paste(DirOut,"/IntermediateFiles/rehh.chr",chr,".map",sep=""),
                          #recode.allele=T,
                          allele_coding = "map",
                          haplotype.in.columns=T,
                          min_maf=maf,
                          min_perc_geno.hap=(1-mind)*100,
                          min_perc_geno.mrk =(1-geno)*100
                        )
    ihh<-scan_hh(haplo)
    scanLIST<-rbind(scanLIST,ihh)
  }
  print(paste("normalize ihh by freq bin of 0.1 for ",pop))
  ihs<-ihh2ihs(scanLIST,freqbin=0.1,min_maf = 0)
  table<-na.omit(cbind(row.names(ihs$ihs),ihs$ihs[,c(1:4)]))
  names(table)<-c("SNP","CHR","POS","IHS","LOGPVALUE")
  
  print(paste("empirical pvalue for ",pop))
  table$iHS_p<-(rank(as.numeric(table$IHS),ties.method="min")+1)/(dim(table)[1]+1)
  
  table$abs_iHS<-abs(table$IHS)
  table$abs_iHS_p<-(rank(-as.numeric(table$abs_iHS),ties.method="min")+1)/(dim(table)[1]+1)
  if(sum(table$SNP %in% row.names(mapINFO))!=nrow(table)){
    stop("pb getting back physical position")
  }
  table$SNP<-as.character(table$SNP)
  if(sum(table$SNP == mapINFO[table$SNP,"id"])!=nrow(table)){
	stop("pb getting back physical position 2")
  }
  print(head(table))
  table$POS<-mapINFO[table$SNP,"pos"]
  print(head(table))
  write.table(table[,c("SNP","CHR","POS","IHS","iHS_p")],paste(DirOut,"/",pop,".IHS",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
  write.table(table[,c("SNP","CHR","POS","abs_iHS","abs_iHS_p")],paste(DirOut,"/",pop,".absIHS",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
  png(paste(DirOut,"/",pop,".IHS.Distribution.png",sep=""))
  distribplot(data=as.numeric(table$IHS),main="iHS distribution")
  dev.off()

}

