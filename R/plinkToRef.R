args<-commandArgs(trailingOnly = T)

plinkFile<-unlist(args[1])

if(length(args)==0){

    print("Arguments have to be supplied: ")
    print("1. plinkFile to turn in to ref panel, 2. if remove Dups (1: yes, 0: no (default)), 3. Maf filter")
    q()
}

rmDups<-0
maf<-0

if(length(args)>1){

    if(!as.numeric(args[2])%in%c(0,1)){
        print("2nd argument has to be  if remove Dups, and has to be 0 (no) or 1 (yes), 3rd argument is MAF cutoff")
        q()
    }
    else{
        rmDups<-as.numeric(args[2])
    } 

}

if(length(args)>2){
  maf<-as.numeric(args[3])
  if(maf>=0.5){
    cat("Has to be minor allele frequency - meaning < 0.5")
    q() 
  }
}

if(grepl("~",plinkFile)){
  print("Migt not work with '~' in path to plink files")
}

require(snpStats)

plinkV2<-function(plinkFile){
  pl<-snpStats::read.plink(plinkFile)
  pl2<-matrix(as(pl$genotypes,"numeric"),nrow=nrow(pl$genotypes),ncol=ncol(pl$genotypes))
  colnames(pl2)<-colnames(pl$genotypes)
  snp<-colnames(pl2)
  ##    geno<-as.integer(as.integer(pl)-1)
  ##dim(geno)<-dim(pl)
  bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
  fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
  rownames(pl2)<-fam$V2
  ind<-rownames(pl2)
  list(geno=pl2,bim=bim,fam=fam)
}

pl<-plinkV2(plinkFile)

## 0 major/major 1 major/minor 2 minor/minor
## fam: FamilyID IndividualID
l<-list()
for(pop in unique(pl$fam$V1)){
    indis<-pl$fam[ pl$fam$V1%in%pop,2]
    y<-pl$geno[indis,]
    if(class(y)=="numeric"){
        ## calculates freq from values without NAs
        l[[pop]]<-1-sum(y,na.rm=T)/(sum(!is.na(y))*2)
    } else{
        l[[pop]]<-1-colSums(y,na.rm=T)/(colSums(!is.na(y))*2)
    }
}
f2<-do.call(cbind,l)


## ref looks like this, example with French, Han, Karitiana, Papuan and Yoruba
## id chr pos name A0_freq A1 French Han Karitiana Papuan Yoruba
if(length(args)>2 & rmDups>0){
    ## bim has same number of sites as geno with same ordering
    alleles<-apply(pl$bim,1,function(x) paste(sort(x[5:6]),collapse="_"))
    idKeep<-!duplicated(paste0(pl$bim$V1,"_",pl$bim$V4,"_",alleles))
    keep2<-apply(pl$geno,2,function(x) sum(x,na.rm = T)/(2*sum(!is.na(x))) >= maf & sum(x,na.rm = T)/(2*sum(!is.na(x))) <= 1-maf)
    keep<-keep2 & idKeep
    f3<-cbind(id=as.vector(paste(trimws(pl$bim$V1[keep]),trimws(pl$bim$V4[keep]),sep="_")),chr=pl$bim$V1[keep],pos=pl$bim$V4[keep],name=pl$bim$V2[keep],A0_freq=pl$bim$V5[keep],A1=pl$bim$V6[keep],f2[keep,])
} else if(length(args)>1 & rmDups>0){
    ## removes sites based on frequncy of site across all pops
    alleles<-apply(pl$bim,1,function(x) paste(sort(x[5:6]),collapse="_"))
    idKeep<-!duplicated(paste0(pl$bim$V1,"_",pl$bim$V4,"_",alleles))
    keep<-idKeep
    f3<-cbind(id=as.vector(paste(trimws(pl$bim$V1[keep]),trimws(pl$bim$V4[keep]),sep="_")),chr=pl$bim$V1[keep],pos=pl$bim$V4[keep],name=pl$bim$V2[keep],A0_freq=pl$bim$V5[keep],A1=pl$bim$V6[keep],f2[keep,])
} else{
    
    f3<-cbind(id=as.vector(paste(trimws(pl$bim$V1),trimws(pl$bim$V4),sep="_")),chr=pl$bim$V1,pos=pl$bim$V4,name=pl$bim$V2,A0_freq=pl$bim$V5,A1=pl$bim$V6,f2)  
}


nInd<-rbind(names(table(pl$fam$V1)),as.vector(table(pl$fam$V1)))
if(grepl("/",plinkFile)){
    name<-tail(unlist(strsplit(plinkFile,"/")),1)
} else{
    name<-plinkFile
}
sites<-cbind(f3[,2:3],f3[,5:6])

write.table(f3,paste("refPanel_",name,".txt",sep=""),col=T,row=F,quote=F)
write.table(nInd,paste("nInd_",name,".txt",sep=""),col=F,row=F,quote=F)
write.table(sites,paste(name,".sites",sep=""),col=F,row=F,quote=F)

