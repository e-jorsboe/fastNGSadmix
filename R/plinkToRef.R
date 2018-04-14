args<-commandArgs(trailingOnly = T)

plinkFile<-unlist(args[1])

if(length(args)==0){

    print("Arguments have to be supplied: ")
    print("1. plinkFile to turn in to ref panel, 2. if remove Dups (1: yes (default), 0: no), 3. Maf filter")
    print("It is HIGHLY recommended to do missingness filter on data prior to creating a ref panel")
    q()
}

rmDups<-1
maf<-0

if(length(args)>1){

    if(!as.numeric(args[2])%in%c(0,1)){
        print("2nd argument has to be if remove Dups, and has to be 0 (no) or 1 (yes), 3rd argument is MAF cutoff")
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

pl<-snpStats::read.plink(plinkFile)


## 0 major/major 1 major/minor 2 minor/minor
## fam: FamilyID IndividualID
l<-list()
for(pop in unique(pl$fam[,1])){
    indis<-pl$fam[ pl$fam[,1]%in%pop,2]
    y<-as(pl$genotypes[indis,],"numeric")
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
    alleles<-apply(pl$map,1,function(x) paste(sort(x[5:6]),collapse="_"))
    idKeep<-!duplicated(paste0(pl$map[,1],"_",pl$map[,4],"_",alleles))
    geno<-as(pl$genotypes,"numeric")
    keep2<-apply(geno,2,function(x) sum(x,na.rm = T)/(2*sum(!is.na(x))) >= maf & sum(x,na.rm = T)/(2*sum(!is.na(x))) <= 1-maf)
    keep<-keep2 & idKeep
    f3<-cbind(id=as.vector(paste(trimws(pl$map[,1][keep]),trimws(pl$map[,4][keep]),sep="_")),chr=pl$map[,1][keep],pos=pl$map[,4][keep],name=pl$map[,2][keep],A0_freq=pl$map[,5][keep],A1=pl$map[,6][keep],f2[keep,])
} else if(length(args)>1 & rmDups>0){
    ## removes sites based on frequncy of site across all pops
    alleles<-apply(pl$map,1,function(x) paste(sort(x[5:6]),collapse="_"))
    idKeep<-!duplicated(paste0(pl$map[,1],"_",pl$map[,4],"_",alleles))
    keep<-idKeep
    f3<-cbind(id=as.vector(paste(trimws(pl$map[,1][keep]),trimws(pl$map[,4][keep]),sep="_")),chr=pl$map[,1][keep],pos=pl$map[,4][keep],name=pl$map[,2][keep],A0_freq=pl$map[,5][keep],A1=pl$map[,6][keep],f2[keep,])
} else{
    
    f3<-cbind(id=as.vector(paste(trimws(pl$map[,1]),trimws(pl$map[,4]),sep="_")),chr=pl$map[,1],pos=pl$map[,4],name=pl$map[,2],A0_freq=pl$map[,5],A1=pl$map[,6],f2)  
}


nInd<-rbind(names(table(pl$fam[,1])),as.vector(table(pl$fam[,1])))
if(grepl("/",plinkFile)){
    name<-tail(unlist(strsplit(plinkFile,"/")),1)
} else{
    name<-plinkFile
}
sites<-cbind(f3[,2:3],f3[,5:6])

write.table(f3,paste("refPanel_",name,".txt",sep=""),col=T,row=F,quote=F)
write.table(nInd,paste("nInd_",name,".txt",sep=""),col=F,row=F,quote=F)
write.table(sites,paste(name,".sites",sep=""),col=F,row=F,quote=F)

