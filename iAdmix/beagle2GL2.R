GLfile<-commandArgs(trailingOnly = T)[1]

GL<-read.table(GLfile,as.is=T,h=T)

GL2<-cbind(log10(GL[,4]),log10(GL[,5]),log10(GL[,6]),GL[,1],GL[,2],GL[,3])

name<-unlist(strsplit(GLfile,"\\.txt$"))[1]

write.table(GL2,paste(name,".GL2",sep=""),col=F,row=F,quote=F)
