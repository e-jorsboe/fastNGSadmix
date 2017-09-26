## du -cksh *

########### do not change ################3
l<-commandArgs(TRUE)
getArgs<-function(x,l){
    unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
}
Args<-function(l,args){
    if(length(l)%%2!=0){
        cat("Error -> not all options have an argument!\n")
        q("no")       
    }
    if(sum(grepl("^-",l))!=length(l)/2){
        cat("Error -> not correct number of options - or argument starting with -\n")
        q("no")       
    }
    ## this assumes -option argument structure!
    l<-sapply(seq(1,(length(l)-1),2),function(x) paste0(gsub("^-","",l[x]),"=",l[x+1]))
    
    if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
        cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument")
        q("no")
    }
    arguments<-list()
    for(a in names(args)){
        arguments[[a]]<-getArgs(a,l)
    }
    
  ## check for plinkFile or beagle file
  if(!any(c("plinkFile","likes")%in%names(arguments))){
    cat("Error -> plinkFile or likes argument has to be supplied!\n")
    q("no")  
  } else if(all(c("plinkFile","likes")%in%names(arguments))){
    cat("Error -> plinkFile and likes argument have both been supplied, only one please!\n")
    q("no")  
  } else if(!all(c("qopt")%in%names(arguments))){
    cat("Error -> estimated admixture proporotions have to be supplied from fastNGSadmix, as .qopt file!\n")
    q("no")  
  } else if(!all(c("ref")%in%names(arguments))){
    cat("Error -> genotypes of reference populations have to be supplied as plinkFile binary files!\n")
    q("no")    
  } else if(c("plinkFile")%in%names(arguments)){
    arguments$likes<-""
  } else if(c("likes")%in%names(arguments)){
    arguments$plinkFile<-""
  } 
  
  if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
    cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
    q("no")
  }
  
  for(a in names(args))
    if(is.null(arguments[[a]]))
      arguments[[a]]<-args[[match(a,names(args))]]
  return(arguments)
}

print.args<-function(args,des){
  if(missing(des)){
    des<-as.list(rep("",length(args)))
    names(des)<-names(args)
  }
  cat("->  needed arguments:\n")
  mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
  cat("->  optional arguments (defaults):\n")
  mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
  q("no")
}
###### ####### ###### ###### ###### #######
## choose your parameters and defaults
## NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments

## getting arguments for run

args<-list(likes=NULL,
           plinkFile=NULL,
           ## called genotype of reference panel 
           ref = NULL, 
           qopt = NULL,
           dryrun = FALSE,
           out = 'output',
           PCs="1,2",
           multiCores=1,
           saveCovar="0",
           ngsTools="0",
           onlyPrior="0",
           overlapRef="0",
           doPlots="1"
)
## if no argument aree given prints the need arguments and the optional ones with default
des<-list(likes="input GL in beagle format",
          plinkFile="input binary plink in bed format",
          ref = 'plink binary files filename with reference individuals for PCA',
          qopt= "estimated admixture proportions from fastNGSadmix, as .qopt file",
          dryrun = '',
          out= "output filename prefix",
          PCs= "which Principal components to be ploted default 1 and 2",
          multiCores= "using mclapply from the parallel package, denote how many cores to be used, if 1 normal lapply used - for number of cores use: parallel:::detectCores()",
          saveCovar="if covariance matrix of ref individuals should be stored for faster computation, depends on geno file used and pops analysed, 1 or 0 (default)",
          ngsTools="Use ngsTools' (Fumagalli et al., 2013) method where genotypes between individuals are assumed to be independent only give the data, set 1 or 0 (default)",
          onlyPrior="Run analyses with uniform genotype likelihoods, meaning that only the prior is used for inferring the covariance matrix, set 1 or 0 (default)",
          overlapRef="Bases covariance matrix for ref genos on only overlapping markers with input, set 1 or 0 (default)",
          doPlots="Option for if admixture plot and PCA plot should be generated or just covararinace matrix, eigenvectors and eigenvalues, set 1 (default) or 0:"
          
)

######################################
#######get arguments and add to workspace
### do not change
if(length(l)==0) print.args(args,des)
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
  cat("Arguments: output prefix\n")
  q("no")
}
###################################


if(!require(snpStats)){
    print("You must install the R package: 'snpStats'")
    if(as.numeric(multiCores) > 1 & !require(parallel)){
        print("You must install the R packages: 'parallel' and 'snpStats'")
    }
}

## for reading plink files using snpStats
plinkV2<-function(plinkFile){
  pl<-snpStats::read.plink(plinkFile)
  pl2<-matrix(methods::as(pl$genotypes,"numeric"),nrow=nrow(pl$genotypes),ncol=ncol(pl$genotypes))
  colnames(pl2)<-colnames(pl$genotypes)  
  ind<-rownames(pl2)
  snp<-colnames(pl2)
  bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
  fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
  ## fam has groupID and then individualID
  rownames(pl2)<-fam$V2
  list(geno=pl2,bim=bim,fam=fam,pl=pl)
}

pl<-plinkV2(ref)
admix<-read.table(qopt,h=T,as.is=T)
## refpops are analyzed pops for admixture estimation
refpops<-colnames(admix)

cols<-as.data.frame(cbind(pop=unique(pl$fam$V1),col=as.integer(as.factor(unique(pl$fam$V1)))))

if(dryrun){
  print(table(pl$fam$V1))
  cat('dryrun must be assign to default, FALSE, to execute FastNGSAdmixPCA\n')
  stop()
}

if(!all(k<-refpops%in%unique(pl$fam$V1))){   
    cat("These are not part of the genos:\n")
    print(refpops[!k])
    cat("use these pops from the genos instead\n")
    print(unique(pl$fam$V1))
    stop()
}


if(!all(c(saveCovar,ngsTools,onlyPrior,overlapRef,doPlots)%in%c("0","1"))){
    cat("saveCovar, ngsTools, onlyPrior, overlapRef or doPlots arguments must be 0 or 1\n")    
    stop()        
}
        

ccol <- c("darkgreen","darkorange","goldenrod2","#A6761D","darkred","lightgreen","darkblue","lightblue")
grDevices::palette(ccol)
gar<-grDevices::dev.off()


## used for calculating the covariance entries between input data and ref indis without normalizing
glfunc <- function(x,G_mat,my2,pre_norm,geno_test2) {
  freq <- my2/2
  ## calculates (g0-2f)*(g'-2f)*P(g0|X,h), (g1-2f)*(g'-2f)*P(g1|X,h), (g2-2f)*(g'-2f)*P(g2|X,h)
  ## that is enough as P(g'|X,h) != only when g' is actual genotype, no uncertainty
  abc <- (G_mat-my2)*(geno_test2[,x]-my2)*pre_norm
  ## takes sum for each site and divides by 2f(1-f)
  abcr <- (rowSums(abc))/((2*freq*(1-freq)))
  ## then takes sum for all sites and divides by number of sites
  return(sum(abcr)/length(my2))
}

## generates barplot of admixture proportions with conf intervals
generateBarplot<-function(admix,sorting,out){
    margins<-c(5.1, 4.1, 8.1, 2.1)
    admix<-admix[,match(sorting,colnames(admix))]
    if(nrow(admix)>10){
        
        m<-matrix(0,nrow=2,ncol=ncol(admix))
        m[1,]<-as.numeric(apply(admix,2,function(x) quantile(x[2:length(x)],probs=c(0.025))))
        m[2,]<-as.numeric(apply(admix,2,function(x) quantile(x[2:length(x)],probs=c(0.975))))
        
        bitmap(paste(out,"_","quantile_","admixBarplot.png",sep=""),res=300)
        par(mar=margins)
        b1<-barplot(as.numeric(admix[1,]),ylab="admixture proportion",col=cols[ cols$pop%in%colnames(admix),"col"],ylim=c(0,1.1))
        
        segments(b1,m[1,],b1,m[2,])
        segments(b1-0.2,m[1,],b1+0.2,m[1,])
        segments(b1-0.2,m[2,],b1+0.2,m[2,])
        par(xpd=T)
        legend("topright",inset=c(0.0,-0.2),legend=cols[ cols$pop%in%colnames(admix),"pop"],fill=cols[ cols$pop%in%colnames(admix),"col"],cex=1.5)
        garbage<-dev.off()
    } else{
        
        bitmap(paste(out,"_admixBarplot.png",sep=""),res=300)
        par(mar=margins)
        b1<-barplot(as.numeric(admix[1,]),ylab="admixture proportion",col=cols[ cols$pop%in%colnames(admix),"col"],ylim=c(0,1.1))
        par(xpd=T)
        legend("topright",inset=c(0.0,-0.2),legend=cols[ cols$pop%in%colnames(admix),"pop"],fill=cols[ cols$pop%in%colnames(admix),"col"],cex=1.5)
        
        garbage<-dev.off()
    }
}

estimateAdmixPCA<-function(likes=NULL,plinkFile=NULL,admix,refpops,out){

    ## if plink files reads and convert to beagle file
    if(plinkFile!=""){
        plInput<-plinkV2(paste(plinkFile,sep=""))
        GL.raw<-cbind(paste(plInput$bim$V1,plInput$bim$V4,sep="_"),plInput$bim$V6,plInput$bim$V5,0,0,0)
        GL.raw[ which(plInput$geno[1,]==2),4]<-1
        GL.raw[ which(plInput$geno[1,]==1),5]<-1
        GL.raw[ which(plInput$geno[1,]==0),6]<-1
        GL.raw<-GL.raw[ !is.na(plInput$geno[1,]),]
        GL.raw2<-as.data.frame(GL.raw,stringsAsFactors=F)
        colnames(GL.raw2)<-c("marker", "allele1", "allele2", "Ind0", "Ind0.1", "Ind0.2")
        
    } else{
        GL.raw2 <- read.table(paste(likes,sep=""),as.is=T,h=T,colC=c("character","integer","numeric")[c(1,1,1,3,3,3)])
    }
    if(any(duplicated(GL.raw2[,1]))){    
        print("Duplicate markers in beagle or plinkFile input file - fix this!")
        stop()
    }

    if(any(duplicated(paste(pl$bim$V1,pl$bim$V4,sep="_")))){    
        print("Duplicate markers in reference panel genotypes - fix this!")
        stop()
    }

    
    ## overlapping sites with ref genos
    rownames(GL.raw2) <- GL.raw2[,1]
   
    out1<-dirname(out)
    ref1<-basename(ref)
    ind <- pl$fam[ pl$fam$V1%in%refpops,"V1"]
    table(ind)

    covarFilename<-paste0(out1,"/",ref1,paste0(refpops,collapse=""),".Rdata")
    if(!file.exists(covarFilename)){

        geno_test<-pl$geno[ pl$fam[  pl$fam$V1%in%refpops,"V2"],]
        if(overlapRef=="1"){
            print(paste0("Only using overlap of markers for ref genos - cannot save covariance matrix!"))
            keep<-paste(pl$bim$V1,pl$bim$V4,sep="_")%in%GL.raw2[,1]
            geno_test<-geno_test[ ,keep]            
        }
        
        ## first PCA for ref pops based on all SNPs
        ##snp row::sample col
        geno_test <- t(geno_test)
        ## too many NA, so put NA to 2 (major major) instead of removing column
        geno_test[is.na(geno_test)] = 2 
        my <- rowMeans(geno_test,na.rm=T)
        freq<-my/2            
        keep<-freq>0 & freq < 1
        geno_test<- geno_test[keep,]
        freq<-freq[keep]
        my<-my[keep]        
        ##normalizing the genotype matrix, 
        ## sart because we square both de- and nominator in matrix multi
        M <- (geno_test-my)/sqrt(2*freq*(1-freq))      
        ##M[is.na(M)] <- 2
        ##get the (almost) covariance matrix
        print("Calculating covarinace matrix for reference individuals")
        print("")
        Xtmp<-crossprod(M,M)
        ##Xtmp<-(t(M)%*%M)
        ## normalizing the covariance matrix
        X<-(1/nrow(geno_test))*Xtmp
        if(saveCovar=="1" & overlapRef!="1"){
            print(paste0("Saving covariance matrix of ref individuals to: ",covarFilename))
            print(paste0("This can be disabled by setting saveCovar=0"))
            save(X,file=covarFilename)
        }
    } else{
        load(paste0(covarFilename))
    }
    
    GL.raw2<-GL.raw2[ GL.raw2[,1]%in%paste(pl$bim$V1,pl$bim$V4,sep="_"),]
    bim2<-pl$bim[ paste(pl$bim$V1,pl$bim$V4,sep="_")%in%GL.raw2[,1],]
    geno2<-pl$geno[ ,colnames(pl$geno)%in%bim2$V2]

    ## makes sure beagle or plinkFile input file ordered as reference genotypes
    GL.raw2<-GL.raw2[order(match(GL.raw2[,1],paste(bim2[,1],bim2[,4],sep="_"))),]
    
    ## if alleles coded as 0,1,2,3 instead of A,C,G,T
    if(any(c(0,1,2,3)%in%GL.raw2[,2]) | any(c(0,1,2,3)%in%GL.raw2[,3])){
        GL.raw2[,2]<-sapply(GL.raw2[,2], function(x) ifelse(x==0,"A",ifelse(x==1,"C",ifelse(x==2,"G",ifelse(x==3,"T",x)))))
        GL.raw2[,3]<-sapply(GL.raw2[,3], function(x) ifelse(x==0,"A",ifelse(x==1,"C",ifelse(x==2,"G",ifelse(x==3,"T",x)))))
        
    }
    
    print("The overlap between input and genos is:")
    print(ncol(geno2))
    print("")
    
    ## those were alleles agree should be flipped like for refPanel, so all genotypes point in same direction
    flip<-GL.raw2[,2]==bim2[,6] & GL.raw2[,3]==bim2[,5]
    geno_test2<-geno2[ pl$fam[  pl$fam$V1%in%refpops,"V2"],]
        
    ## hereby only constructing ref panel of individuals/pops in admix file
    popFreqs<-sapply(colnames(admix), function(x) colMeans(geno_test2[pl$fam[ pl$fam[,1]==x,2],],na.rm=T)/2)
    popFreqs2<-popFreqs[,match(colnames(popFreqs),colnames(admix))]
    ##snp row::sample col
    geno_test2 <- t(geno_test2) 
    geno_test2[is.na(geno_test2)] = 2  
    geno_test2[flip,]<-2-geno_test2[flip,]
    popFreqs2[!flip,]<-1-popFreqs2[!flip,]

    ## again because some freqs might be zero
    my2 <- rowMeans(geno_test2,na.rm=T)
    freq2<-my2/2

    keep2<-freq2>0 & freq2<1
    GL.raw2<-GL.raw2[keep2,]
    freq2<-freq2[keep2]
    my2<-my2[keep2]
    popFreqs2<-popFreqs2[keep2,]
    geno_test2<-geno_test2[keep2,]

    if(onlyPrior=="1"){
        ## uniform GLs
        GL.raw2[,4:6]<-1/3
    }
    
    ## calculating admixture adjusted freqs
    hj <- as.matrix(popFreqs2) %*% t(admix[1,])
    hj_inv <- 1-hj ## sites X 1
    gs <- cbind(hj**2,2*hj*hj_inv,hj_inv**2) ## Sites X 3
    ## likelihood P(X|G=g)P(G=g|Q,F)
    if(ngsTools=="1"){
        pre <- as.numeric(as.matrix(GL.raw2[,4:6]))*cbind(freq2**2,2*freq2*(1-freq2),(1-freq2)**2)
    } else{
        pre <- as.numeric(as.matrix(GL.raw2[,4:6]))*gs
    }
    ## normalizing likelihoods - the denominator sum(P(X_j|G_j=g)P(G_j=g|h_j))
    pre_norm <- pre/rowSums(pre)
    ## can have issues of pre being only 0's and thereby dividing by 0, yeilding NAs
    pre_norm[is.na(pre_norm)]<-0 
    size <- nrow(pre_norm)
    G_mat <- data.frame(x=rep(0,size),y=rep(1,size),z=rep(2,size))
    
    ## calculating covariances between input individual and ref individuals

    print("Calculating covarinace matrix for input individual")
    print("")
    if(as.numeric(multiCores)>1){
        GL_called <- unlist(parallel:::mclapply(colnames(geno_test2),glfunc,G_mat=G_mat,my=my2,pre_norm=pre_norm,geno_test=geno_test2,mc.cores=as.numeric(multiCores)))
    } else{
        GL_called <- unlist(lapply(colnames(geno_test2),glfunc,G_mat=G_mat,my=my2,pre_norm=pre_norm,geno_test=geno_test2))
    }
    
    ## calculating covariances, between individual itself
    abc_single <- (G_mat-my2)*(G_mat-my2)*pre_norm
    
    ## takes sum for each site and divides by 2f(1-f)
    abcr_single <- (rowSums(abc_single))/((2*freq2*(1-freq2)))
    ## takes big sum and divides by number of sites
    abcr_single <-sum(abcr_single)/length(my2)
    
    GL_called_diag <- c(as.numeric(GL_called),abcr_single)  
    ## normalizing input, putting on last row of data, X is normalized covariance matrix only from ref genos
    X_1 <- rbind(X,'SAMPLE'=as.numeric(GL_called))
    X_2 <- as.data.frame(cbind(X_1,SAMPLE=GL_called_diag))   
    ## final normalized covariance matrix
    X_norm <- X_2  
    return(list(covar = X_norm,indi=ind))
    
}




PCAplotV2 = function(cova,ind,admix,out,PCs,onlyPlot=F) {
    ## eigen decomposition of covariance matrix for PCA
    E<-eigen(cova)
    ## writes out eigenvectors and values
    write.table(cbind(rownames(cova),E$vectors),paste0(out,'_eigenvecs.txt'),row=T,quote=F,col=F)
    write.table(E$values,paste0(out,'_eigenvals.txt'),row=T,quote=F,col=F)
    if(onlyPlot){
        return()
    }
    ## extracts chosen PCs
    PC_12 <- round(as.numeric(E$values/sum(E$values))[PCs],3)*100
    a <- data.frame(E$vectors[,PCs])
    a$pop <- c(ind,'SAMPLE')
    colnames(a) <- c(paste('PC',1,sep=""),paste('PC',2,sep=""),'pop')
    pdf(paste0(out,'_PCAplot.pdf'))
    par(mar=c(5, 4, 4, 8) + 0.1)

    ## this maps pop into some colour number
    pcaColours<-sapply(a$pop[1:(nrow(a)-1)], function(x) cols[ cols$pop==x,"col"])
    
    plot(a$PC1[1:(nrow(a)-1)],a$PC2[1:(nrow(a)-1)],xlab=paste('PC',PCs[1],' (%)',PC_12[PCs[1]]),ylab=paste('PC',PCs[2],' (%)',PC_12[2]),col=pcaColours,pch=16,ylim=c(min(a$PC2),max(a$PC2)),xlim=c(min(a$PC1),max(a$PC1)))
    points(a$PC1[nrow(a)],a$PC2[nrow(a)],pch=4,cex=2,lwd=4)
    print("Input individual ('SAMPLE') is plotted at in PCA plot:")
    print(paste(a$PC1[nrow(a)],a$PC2[nrow(a)]))
    par(xpd=TRUE)
    legend("topright",inset=c(-0.3,0),legend=unique(a$pop[1:(nrow(a)-1)]),fill=unique(pcaColours))
    garbage<-dev.off()
}

pop_list<- estimateAdmixPCA(likes=likes,plinkFile=plinkFile,admix=admix,refpops = refpops,out = out)
write.table(pop_list$covar, file=paste0(out,'_covar.txt'),quote=F)
write.table(cbind(rownames(pop_list$covar),c(pop_list$ind,"SAMPLE")), file=paste0(out,'_indi.txt'),quote=F,col=F,row=F)
PCAplotV2(pop_list$covar,pop_list$indi,admix=admix,out,PCs=sort(as.numeric(unlist(strsplit(PCs,",")))),onlyPlot=(doPlots=="1"))
if(doPlots=="1"){
    generateBarplot(admix=admix,sorting = unique(pop_list$ind),out = out)    
}
