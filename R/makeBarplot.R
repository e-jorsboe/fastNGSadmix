
args<-commandArgs(trailingOnly = T)

if(length(args)==0){

    print("Arguments have to be supplied: ")
    print("1. .qopt file to plot")
    q()
}

qoptFile<-args[1]

if(grepl("~",qoptFile)){
  print("Migt not work with '~' in path to .qopt file")
}

ccol <- c("darkgreen","darkorange","goldenrod2","#A6761D","darkred","lightgreen","darkblue","lightblue")
grDevices::palette(ccol)
gar<-grDevices::dev.off()

admix<-read.table(qoptFile,h=T,as.is=T)
## refpops are analyzed pops for admixture estimation
refpops<-colnames(admix)
cols<-as.data.frame(cbind(pop=colnames(admix),col=as.integer(as.factor(colnames(admix)))))

generateBarplot<-function(admix,out){
    margins<-c(5.1, 4.1, 8.1, 2.1)   
    if(nrow(admix)>10){
        
        m<-matrix(0,nrow=2,ncol=ncol(admix))
        m[1,]<-as.numeric(apply(admix,2,function(x) quantile(x[2:length(x)],probs=c(0.025))))
        m[2,]<-as.numeric(apply(admix,2,function(x) quantile(x[2:length(x)],probs=c(0.975))))
        
        bitmap(paste(out,"_","quantile_","admixBarplot.png",sep=""),res=300)
        par(mar=margins)
        b1<-barplot(as.numeric(admix[1,]),ylab="admixture proportion",col=cols$col,ylim=c(0,1.1))
        
        segments(b1,m[1,],b1,m[2,])
        segments(b1-0.2,m[1,],b1+0.2,m[1,])
        segments(b1-0.2,m[2,],b1+0.2,m[2,])
        par(xpd=T)
        legend("topright",inset=c(0.0,-0.2),legend=cols$pop,fill=cols$col,cex=1.5)
        garbage<-dev.off()
    } else{
        
        bitmap(paste(out,"_admixBarplot.png",sep=""),res=300)
        par(mar=margins)
        b1<-barplot(as.numeric(admix[1,]),ylab="admixture proportion",col=cols$col,ylim=c(0,1.1))
        par(xpd=T)
        legend("topright",inset=c(0.0,-0.2),legend=cols$pop,fill=cols$col,cex=1.5)
        
        garbage<-dev.off()
    }
}


generateBarplot(admix=admix,out = unlist(strsplit(basename(qoptFile),".qopt"))[1])
