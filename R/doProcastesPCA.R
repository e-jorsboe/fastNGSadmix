arguments<-commandArgs(trailingOnly=T)

## will take all _covar.txt files in that diretory
dir<-arguments[1]
## should ONLY be with individuals you want to plot, SHOULD be fam used for generating PCA
famFile<-arguments[2]
## Should have column of IDs  which is same as filename of *_covar.txt, way to tell which samples should be included
infoFile<-arguments[3]
## name of group of smaples that you want to plot together via procrustes
groupName<-arguments[4]

if(length(arguments)==0){

    print("Arguments have to be supplied: ")
    print("1. directory of *_covar.txt files, 2. .fam file of plink file used for ref individuals (geno=), 3. file of individuals to be included, must have same name as *_covar.txt files (except _covar.txt), 4. name of group of *_covar.txt files to be plotted")
    q()
}


fam<-read.table(famFile,as.is=T)
## translate pops to colours like in fastNGSadmixPCA
cols<-as.data.frame(cbind(pop=unique(fam$V1),col=as.integer(as.factor(unique(fam$V1)))))

files2<-list.files(dir,pattern="_covar.txt")
files2<-sapply(files2,function(x) unlist(strsplit(x,"_covar.txt"))[1])
info<-read.table(infoFile,as.is=T)

files<-list.files(dir,pattern="_covar.txt",full=T)

covar<-read.table(files[1],as.is=T,h=T)
 
## the ref PCA has order of individuals found in _covar.txt files
fam<-fam[ fam$V2%in%rownames(covar),]
## colours of PCA based on order of appearance in .fam file

fam<-fam[order(match(fam$V2,rownames(covar)[1:(length(covar)-1)])),]

filesName <- as.vector(sapply(files, function(x) unlist(strsplit(basename(x),"_covar.txt")[1])))

files2 <- cbind(files=files,filesName=filesName)

int<-intersect(info[,1],files2[,"filesName"])

info[,1]<-info[ info[,1]%in%int,]
files2<-files2[ files2[,"filesName"]%in%int,]

print("creating tmp procrustesPCs directory")
## returns 1 if dir already exists and 0 if not
createdError<-system("mkdir procustesPCs")

pc <- function(f){
  r<-read.table(f,head=T)
  e<-eigen(as.matrix(r))
  outname <-  unlist(strsplit(basename(f),"_covar.txt")[1])
  ##outname <- sapply(strsplit(basename(f),".sor"),function(x)x[1])
  write.table(e$vectors[,1:2],paste0("procustesPCs/",outname),col=F,row=F,qu=F)
}

## look at how many cores avaible

cores<-parallel:::detectCores()
a<-parallel::mclapply(files2[,"files"],pc,mc.cores=min(ceiling(cores/2),20))

##########################################
library(vegan)
pcs<-list.files("procustesPCs",full=T)

pcs2 <-  sapply(pcs,function(x) basename(x))

print("This many files were succesfully read and will be procustered:")
print(length(pcs))

m0<-read.table(pcs[1])

m0<-m0[1:(nrow(m0)-1),]

## sumperimposes PCA points of same individuals for first PCA pcs[1], and then predicts PCA of new individual
pro <- function(x){
  m2<-read.table(x)
  
  r<-procrustes(m0,m2[1:(nrow(m2)-1),])
  predict(r,m2[nrow(m2),])
}
res <- t(sapply(pcs,pro))

ylim<-range(c(res[,2],m0[,2]))

## same colours as PCA
ccol <- c("darkgreen","darkorange","goldenrod2","#A6761D","darkred","lightgreen","darkblue","lightblue")
grDevices::palette(ccol)
gar<-grDevices::dev.off()

## pca coordinates of ref individuals with pop label
r1<-cbind(fam$V1,fam$V2,m0)
## pca coordinates of vikings with id
r2<-cbind(rep(groupName,length(pcs)),sapply(basename(pcs),function(x) gsub("_covar.txt","",x)),res)
rownames(r2)<-NULL
colnames(r2)<-colnames(r1)
pcaPoints<-rbind(r1,r2)

bitmap(paste0("procrustesPCA",groupName,".png"),res=300)
write.table(pcaPoints,paste0("procrustesPCA",groupName,".txt"),col=F,row=F,qu=F)

## should we also be able to have different labels for the samples being procrustered on?

pcaColours<-sapply(fam$V1, function(x) cols[ cols$pop==x,"col"])

plot(m0,col=pcaColours,lwd=2,ylim=ylim,ylab="PC2",xlab="PC1",main=paste0("procrustes PCA with ",groupName))
points(res,pch=4)

legend("bottomright",cex=1,pch=c(rep(15,length(unique(fam$V1))),4),col=c(unique(pcaColours),"black"),legend=c(paste0(unique(fam$V1)),groupName))

dev.off()
if(!createdError){
    system("rm -r procustesPCs")
}
