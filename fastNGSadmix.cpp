/*
  log:
  g++ fastNGSadmix.cpp -lz -lpthread  -O3 -o fastNGSadmix

  log: (with readplink function)
  g++ fastNGSadmix.cpp readplink.c -lz -lpthread  -O3 -o fastNGSadmix


  debug:
  g++ fastNGSadmix.cpp -lz -lpthread -ggdb -O3 -o fastNGSadmix

  debug: (with readplink function)
  g++ fastNGSadmix.cpp readplink.c -lz -lpthread -ggdb -O3 -o fastNGSadmix

*/
 
//optimazation parameteres for maf estimation
#define MAF_START 0.3
#define MAF_ITER 20

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <limits>
#include <zlib.h>
#include <vector>
#include <pthread.h>
#include <signal.h>
#include <vector>
#include <sys/stat.h>
// not really necessary any longer
//#include "readplink.h"

//This is taken from here:
//http://blog.albertarmea.com/post/47089939939/using-pthread-barrier-on-mac-os-x
#ifdef __APPLE__

#ifndef PTHREAD_BARRIER_H_
#define PTHREAD_BARRIER_H_

#include <pthread.h>
#include <errno.h>

/////////////////////////////////////////////////////////////////////
// globale ting: functioner + variable(ikke god stil)


#endif // PTHREAD_BARRIER_H_
#endif // __APPLE__



//global stuff below, this is very beautifull
char **keeps=NULL;
#define LENS 100000 //this is the max number of bytes perline, should make bigger
int SIG_COND =1;//if we catch signal then quit program nicely
double tol=1e-5; //stopping criteria


double errTolMin=1e-9;
double errTolStart=0.1;
double errTol=errTolStart;//frequencies and admixture coef cannot be less than this or more than 1-this
double misTol=0.05;
// for reading in plink files, so that GL has no 0 values
double seqError=1e-04;


double **allocDouble(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}

void dalloc(double **ret,size_t x){
  for(size_t i=0;i<x;i++){
    delete [] ret[i] ;
  }
  delete [] ret;
}

// for the emAcc, in order to avoid static pointers
typedef struct{
 double **F_em1;
 double *Q_em1;
 double **F_diff1;
 double *Q_diff1;
 double **F_em2;
 double *Q_em2;
 double **F_diff2;
 double *Q_diff2;
 double **F_diff3;
 double *Q_diff3;
 double **F_tmp;
 double *Q_tmp;
}accFQ;

void dallocAccFQ(accFQ &a, int nSites){
  dalloc(a.F_em1,nSites);
  dalloc(a.F_em2,nSites);
  dalloc(a.F_diff1,nSites);
  dalloc(a.F_diff2,nSites);
  dalloc(a.F_diff3,nSites);
  dalloc(a.F_tmp,nSites);
  delete [] a.Q_em1;
  delete [] a.Q_em2;
  delete [] a.Q_diff1;
  delete [] a.Q_diff2;
  delete [] a.Q_diff3;
  delete [] a.Q_tmp;
 
}


accFQ createAccFQ(int nPop, int nSites){
  accFQ a;
  a.F_em1 = allocDouble(nSites,nPop);
  a.Q_em1 = new double[nPop];
  a.F_diff1 = allocDouble(nSites,nPop);
  a.Q_diff1 = new double[nPop];
  a.F_em2 = allocDouble(nSites,nPop);
  a.Q_em2 = new double[nPop];
  a.F_diff2 = allocDouble(nSites,nPop);
  a.Q_diff2 = new double[nPop];
  a.F_diff3 = allocDouble(nSites,nPop);
  a.Q_diff3 = new double[nPop];
  a.F_tmp = allocDouble(nSites,nPop);
  a.Q_tmp = new double[nPop];
  return(a);

}

void minus1dEmil(double *fst,double *sec,size_t x,double *res){
  for(size_t i=0;i<x;i++){
      res[i] = fst[i]-sec[i];
    }
}

void minusEmil(double **fst,double **sec,size_t x,size_t y,double **res){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++){
      res[i][j] = fst[i][j]-sec[i][j];
    }
  }
}

double sumSquareEmil(double **mat,size_t x,size_t y){
  double tmp=0;
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++){
      tmp += mat[i][j]*mat[i][j];
    }
  }
  return tmp;
}

double sumSquare1dEmil(double *mat,size_t x){
  double tmp=0;
  for(size_t i=0;i<x;i++){
      tmp += mat[i]*mat[i];
  }
  return tmp;
}

 
double likeFixedMinor(double p,double *likes,int numInds,char *keepInd){
  // should these actually be normed genotype likelihoods? Or not normed?
  // only used for checking minLrt
  // returns minus the log of the likelihood
  double totalLike=0;
  for(int i=0;i<numInds;i++){
    if(keepInd[i])
      totalLike+=log(likes[i*3+0]*(1-p)*(1-p)+likes[i*3+1]*2.0*p*(1-p)+likes[i*3+2]*p*p);
  }
  return -totalLike;
}


double emFrequency(double *loglike,int numInds, int iter,double start,char *keep,int keepInd){
  if(keepInd == 0)
    return 0.0; 
  float W0;
  float W1;
  float W2;
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;

  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=(loglike[i*3+0])*(1-p)*(1-p);
      W1=(loglike[i*3+1])*2*p*(1-p);
      W2=(loglike[i*3+2])*p*p;
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      if(0&&std::isnan(sum)){
	exit(0);
      }
    }

    p=sum/keepInd;
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }
  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    fprintf(stderr,"keepList (nInd)\n");
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);
    
    for(int ii=0;0&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1)
	fprintf(stderr,"1\t");
      for(int gg=0;gg<3;gg++)
	fprintf(stderr,"%f\t",loglike[ii*3+gg]);
      fprintf(stderr,"\n");
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0){
	// continue jumps to next iteration of loop
        continue;
      }
      W0=(loglike[i*3+0])*(1-p)*(1-p);
      W1=(loglike[i*3+1])*2*p*(1-p);
      W2=(loglike[i*3+2])*p*p;
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
    }
    p=-999;
  } 
  return(p);
}


void map2domainQEmil(double* Q, int nPop){  
  double sum=0;
  for(int k=0;k<nPop;k++){
    if(Q[k]<errTol){
      Q[k] = errTol;
    }
    if(Q[k]>(1-errTol)){
      Q[k] = 1-errTol;
    }
    sum+=Q[k];
  }
  for(int k=0;k<nPop;k++){
    Q[k]=Q[k]/sum;
  }
}

void map2domainFEmil(double** F, int nSites, int nPop){
  for(int s=0;s<nSites;s++)
    for(int k=0;k<nPop;k++){
      if(F[s][k]<errTol){
	F[s][k] = errTol;
      }
      if(F[s][k]>1-errTol){
	F[s][k] = 1-errTol;
      }
    }
}

////////////////////////

int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

std::vector<char *> dumpedFiles;
FILE *openFile(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  if(0&&fexists(c)){//ANDERS DAEMON DRAGON HATES THIS
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  dumpedFiles.push_back(strdup(c));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}

gzFile openFileGz(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  if(0&&fexists(c)){//ANDERS DAEMON DRAGON HATES THIS
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  dumpedFiles.push_back(strdup(c));
  gzFile fp = gzopen(c,"w");
  delete [] c;
  return fp;
}

//some struct will all the data from the beagle file
typedef struct{
  double **genos;
  char *major;
  char *minor;
  // for snp names pointer to pointer
  char **ids;
  int nSites;
  int nInd;
  char **keeps; //matrix containing 0/1 indicating if data or missing
  int *keepInd; //keepInd[nSites] this is the number if informative samples
  float *mafs;
}bgl;



bgl allocBeagle(int nSites){
  bgl b;

  b.nSites = nSites;
  b.major = new char[nSites];
  b.minor = new char[nSites];
  b.ids = new char*[nSites];
  b.nInd = 1;
  b.genos= allocDouble(nSites,3);
  b.keeps = new char*[nSites]; // array nSites x nInd 0 if missing info
  b.keepInd = new int[nSites];
  b.mafs = new float[nSites];

  for(int s=0;SIG_COND&& (s<nSites);s++){
    b.ids[s] = strdup("0");
    b.major[s] = '0';
    b.minor[s] = '0';
    b.keeps[s] = strdup("0");
    b.keepInd[s] = 0;
    b.mafs[s] = 0.0;
  
    for(int i=0;i<3;i++){
      b.genos[s][i] = seqError;
    }
  }
  return(b);
} 
      
//utility function for cleaning up out datastruct
void dallocBeagle(bgl &b){
  for(int i=0;i<b.nSites;i++){
    delete [] b.genos[i];
    //    delete [] b.ids[i];
    free(b.ids[i]);
  }
  delete [] b.minor;
  delete [] b.major;
  delete [] b.genos;
  delete [] b.ids;
}

/*
  plug in random variables

// Q - an nInd x K matrix of population proportions (random)
// F - a nSites x K matrix of frequencies (random)

Q sums to one for every sample
F doesn't do crap
 */

/*
  Returns the bgl struct containing all data from a beagle file.

  It find the nsamples from counting the header
  It finds the number of sites by queing every line in a std::vector
  After the file has been read intotal it reloops over the lines in the vector and parses data
 */

bgl readBeagle(const char* fname) {
  // denoting delimeters in file
  const char *delims = "\t \n";
  gzFile fp = NULL;
  // checking if file can be opened
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  // creating bgl struct
  bgl ret;
  // creaing string of length LENS
  char buf[LENS];
  //find number of columns
  // reading string from compressed file
  gzgets(fp,buf,LENS);
  // a bit like the java Scanner, read line from file
  // cuts a string into tokens, where tokens seperated by delimeter
  strtok(buf,delims);
  int ncols=1;
  // reading first line in order to see nCol
  while(strtok(NULL,delims))
    ncols++;
  if(0!=(ncols-3) %3 ){
    fprintf(stderr,"ncols=%d\n",ncols);
    exit(0);
  }
  ret.nInd = (ncols-3)/3;//this is the number of samples
  //read every line into a vector
  //allocate vector, array like random access, 
  // insertion at end in constant time
  std::vector<char*> tmp;
  // while still some left of string from file
  while(gzgets(fp,buf,LENS)){
    // duplicates string and puts it into tmp vector
    tmp.push_back(strdup(buf));
  }
  //now we now the number of sites
  ret.nSites=tmp.size();
  ret.major= new char[ret.nSites];
  ret.minor= new char[ret.nSites];
  ret.ids = new char*[ret.nSites];
  ret.genos= new double*[ret.nSites];
  //then loop over the vector and parsing every line
  for(int s=0;SIG_COND&& (s<ret.nSites);s++){
    ret.ids[s] = strdup(strtok(tmp[s],delims));
    ret.major[s] = strtok(NULL,delims)[0];
    ret.minor[s] = strtok(NULL,delims)[0];
    ret.genos[s] = new double[3*ret.nInd];
    for(int i=0;i<ret.nInd*3;i++){
      ret.genos[s][i] = atof(strtok(NULL,delims));
      if(ret.genos[s][i]<0){
	fprintf(stderr,"Likelihoods must be positive\n");
	fprintf(stderr,"site %d ind %d geno %d has value %f\n",s,int(i*1.0/3),i%3,ret.genos[s][i]);
	exit(0);
      }
    }
    for(int i=0;i<ret.nInd;i++){
      double tmpS = 0.0;
      for(int g=0;g<3;g++)
	tmpS += ret.genos[s][i*3+g];
      if(!(tmpS>0)){
	fprintf(stderr,"The sum of likelihoods for a genotypes must be positive\n");
	fprintf(stderr,"individual %d site %d has sum %f\n",i,s,tmpS);
	exit(0);
      } 
    }
    free(tmp[s]);
  }
  // Here additional stuff is calculated from the likelihoods
  // this must be done again while filtering later in main
  ret.keeps=new char*[ret.nSites]; // array nSites x nInd 0 if missing info
  ret.keepInd = new int[ret.nSites];
  ret.mafs = new float[ret.nSites];
  for(int s=0;s<ret.nSites;s++){
    ret.keeps[s] = new char[ret.nInd];
    int nKeep =0;
    for(int i=0;i<ret.nInd;i++){
      double mmin=std::min(ret.genos[s][i*3],std::min(ret.genos[s][i*3+1],ret.genos[s][i*3+2]));
      double mmax=std::max(ret.genos[s][i*3],std::max(ret.genos[s][i*3+1],ret.genos[s][i*3+2]));
      if(fabs(mmax-mmin)<misTol)
	ret.keeps[s][i] =0;
      else{
	ret.keeps[s][i] =1;
	nKeep++;
      }
    }
    ret.keepInd[s] = nKeep;
    ret.mafs[s] = emFrequency(ret.genos[s],ret.nInd,MAF_ITER,MAF_START,ret.keeps[s],ret.keepInd[s]);
  }
  gzclose(fp); //clean up filepointer
  return ret;
}

void readDouble(double **d,int x,int y,const char*fname,int neg){
  fprintf(stderr,"opening : %s with x=%d y=%d\n",fname,x,y);
  const char*delims=" \n\t";
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000;
  char buf[lens];
  for(int i=0;i<x;i++){
    if(NULL==fgets(buf,lens,fp)){
      fprintf(stderr,"Increase buffer\n");
      exit(0);
    }
    if(neg)
      d[i][0] = -atof(strtok(buf,delims));
    else
      d[i][0] = atof(strtok(buf,delims));
    for(int j=1;j<y;j++){
      if(neg)
	d[i][j] = -atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
  fclose(fp);
}


/* function to read plink data into bgl format, not necessary anymore
bgl readPlinkToBeagle(char* plinkName) {
  // denoting delimeters in file
  // annoying that pointer to struct
  plink* pl = readplink(plinkName);


  // I do not think it gets that the beagle file has been create
  // and therefore the d.genos is null
  bgl b = allocBeagle(pl->y);
  for(int i=0;i<pl->y;i++){
    if(pl->d[0][i]==0){
      b.genos[i][2]=1.0-seqError;
    } else if(pl->d[0][i]==1){
      b.genos[i][1]=1.0-seqError;
    } else if(pl->d[0][i]==2){
      b.genos[i][0]=1.0-seqError;
    }
  }
  kill_plink(pl);
  return(b);
  

} */


void readDoubleGZ(double **d,int nSites,int nPop,const char*fname,int neg){
  fprintf(stderr,"opening : %s with x=%d y=%d\n",fname,nSites,nPop);
  const char*delims=" \n";
  gzFile fp = NULL;
  if((fp=gzopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000;
  char buf[lens];
  for(int i=0;i<nSites;i++){
    if(NULL==gzgets(fp,buf,lens)){
      fprintf(stderr,"Increase buffer\n");
      exit(0);
    }
    if(neg){
      d[i][0] = -atof(strtok(buf,delims));
    }
    else{
      d[i][0] = atof(strtok(buf,delims));
    }
    for(int j=1;j<nPop;j++){
      if(neg){
	d[i][j] = -atof(strtok(NULL,delims));
      }
      else{
	d[i][j] = atof(strtok(NULL,delims));
      }
    }

  }
  gzclose(fp);
}

void readDouble1d(double *d,int x,const char*fname){
  fprintf(stderr,"opening : %s with x=%d\n",fname,x);
  const char*delims=" \n\t";
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
  if(NULL==fgets(buf,lens,fp)){
    fprintf(stderr,"Increase buffer\n");
    exit(0);
  }
  d[0] = atof(strtok(buf,delims));
  for(int j=1;j<x;j++){
    d[j] = atof(strtok(NULL,delims));
  }
  fclose(fp);
}

void printDouble(double **ret,size_t x,size_t y,FILE *fp){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++)
      fprintf(fp,"%.20f ",ret[i][j]);
    fprintf(fp,"\n");
  }
}

void printDoubleGz(double **ret,size_t x,size_t y,gzFile fp){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++)
      gzprintf(fp,"%.20f ",ret[i][j]);
    gzprintf(fp,"\n");
  }
}

double likelihoodEmil(double* Q, double** F,int nSites, int nPop,double **genos){
  double prod_ind = 0.0;
  for(int j = 0; j < nSites; j++) {
    double *gg = genos[j];
    double freq = 0.0;
    for(int k = 0; k < nPop; k++) {
      freq += (F[j][k])*Q[k];
    }
    // has to be like this, as I sort freqs
    // prior to running this program
    double f = freq;
    double sum = gg[0] * f * f;
    sum += gg[1]*2*f*(1-f);
    sum += gg[2]*(1-f)*(1-f);
    prod_ind += log(sum);
  }
  return prod_ind;
}

void bootstrap(bgl dOrg, bgl &d, double** F_orgOrg, double** F_org, double** F, int nPop) {
  for(int j=0;j<dOrg.nSites;j++){
    // generate random Int from 0 to (nSites-1)
    int row = std::rand() % dOrg.nSites;
    for(int k = 0; k < nPop; k++) {
      F[j][k] = F_orgOrg[row][k];
      F_org[j][k] = F_orgOrg[row][k];
      if(k<=2){
	d.genos[j][k] = dOrg.genos[row][k];
	// pointers getting same adress...
	d.ids[j] = strdup(dOrg.ids[row]);
      }
    }
  }
}

void emUnadjusted(double* Q, double** F, int nSites, int nPop,double **genos,double *Q_1) {
  double sumAG[nPop];
  double sumBG[nPop];
  map2domainFEmil(F, nSites, nPop);
  // do this for each site
  for(int k=0;k<nPop;k++){ //time killar
    sumAG[k]=0;
    sumBG[k]=0;
  }
  for(int j=0;j<nSites;j++){   
    double fpart=0;
    double fpartInv=0;
    double expGG=0;
    for(int k=0;k<nPop;k++){ //time killar
      // admixture adjusted freq, for each pop
      // this has to be F
      fpart += F[j][k] * Q[k];
      fpartInv += (1-F[j][k]) * Q[k];
      // fpartAdj=1-fpartAdj;
      // pre GL (sites x 3) * (adjusted freq)
      double pp0=(fpartInv)*(fpartInv)*genos[j][2];
      double pp1=2*(fpartInv)*fpart*  genos[j][1];
      double pp2=fpart*fpart*        genos[j][0];
      // in order to do the sum
      double sum=pp0+pp1+pp2;
      // for calculating H
      expGG =(pp1+2*pp2)/sum;//range 0-2, this is the expected genotype	  
      // nPop is number of ancestral populations
    }
    for(int k=0;k<nPop;k++){ //time killar
      sumAG[k] += expGG/(fpart) * (Q[k] * F[j][k]); //proteckMe
      sumBG[k] += (2-expGG)/fpartInv * (Q[k] * (1-F[j][k])); //proteckMe
    }
  }
  for(int k=0;k<nPop;k++){ //time killar
    Q_1[k]=(sumAG[k] + sumBG[k])/(2.0*nSites);
  }
  map2domainQEmil(Q_1,nPop);
  // should this be checked before using F_1??
}

// genos - genotype likelihoods
void emEmil(double* Q, double** F, int nSites, double* nInd, int nPop,double **genos,double **F_1,double *Q_1, double** F_org) {
  double sumAG[nPop];
  double sumBG[nPop];
  double sumAGadj[nPop]; // for after freqs adjusted with input
  double sumBGadj[nPop];
  map2domainFEmil(F, nSites, nPop);
  // do this for each site
  for(int k=0;k<nPop;k++){
    sumAGadj[k]=0;
    sumBGadj[k]=0;
  }
  for(int j=0;j<nSites;j++){   
    for(int k=0;k<nPop;k++){ //time killar
      sumAG[k]=0;
      sumBG[k]=0;
    }
    double fpart=0;
    double fpartInv=0;
    double expGG=0;
    for(int k=0;k<nPop;k++){ //time killar
      // admixture adjusted freq, for each pop
      fpart += F[j][k] * Q[k];
      fpartInv += (1-F[j][k]) * Q[k];
      // fpart = 1-fpart;
      // pre GL (sites x 3) * (adjusted freq)
      double pp0=(fpartInv)*(fpartInv)*genos[j][2];
      double pp1=2*(fpartInv)*fpart*  genos[j][1];
      double pp2=fpart*fpart*        genos[j][0];
      // in order to do the sum
      double sum=pp0+pp1+pp2;
      // for calculating H
      expGG=(pp1+2*pp2)/sum;//range 0-2, this is the expected genotype	
      // nPop is number of ancestral populations
    }
    for(int k=0;k<nPop;k++){
      // similar to (H/(q*f))*q, for jth marker
      // H is added to many times :/
      sumAG[k] = (expGG) / (fpart) * (Q[k]*F[j][k]); //time killar
      sumBG[k] = (2-expGG) / fpartInv * (Q[k]*(1-F[j][k]));//time killar
      // adjust with ref panel, so we have input + ref expected number of alleles
      sumAG[k]+=nInd[k]*2*F_org[j][k];
      sumBG[k]+=2*nInd[k]-(2*nInd[k]*F_org[j][k]);
    }
    double fpartAdj=0;
    double fpartAdjInv=0;
    for(int k=0;k<nPop;k++){ //time killar
      // admixture adjusted freq, for each pop
      // this has to be F
      F_1[j][k]=sumAG[k]/(sumAG[k]+sumBG[k]);
      fpartAdj += F_1[j][k] * Q[k];
      fpartAdjInv += (1-F_1[j][k]) * Q[k];
      // fpartAdj=1-fpartAdj;
      // pre GL (sites x 3) * (adjusted freq)
      double pp0=(fpartAdjInv)*(fpartAdjInv)*genos[j][2];
      double pp1=2*(fpartAdjInv)*fpartAdj*  genos[j][1];
      double pp2=fpartAdj*fpartAdj*        genos[j][0];
      // in order to do the sum
      double sum=pp0+pp1+pp2;
      // for calculating H
      expGG =(pp1+2*pp2)/sum;//range 0-2, this is the expected genotype	
      // nPop is number of ancestral populations
    }
    for(int k=0;k<nPop;k++){ //time killar
      sumAGadj[k] += expGG/(fpartAdj) * (Q[k] * F_1[j][k]); //proteckMe
      sumBGadj[k] += (2-expGG)/fpartAdjInv * (Q[k] * (1-F_1[j][k])); //proteckMe
    }
  }
  for(int k=0;k<nPop;k++){ //time killar
    Q_1[k]=(sumAGadj[k] + sumBGadj[k])/(2.0*nSites);
  }
  map2domainQEmil(Q_1,nPop);
  // should this be checked before using F_1??
  map2domainFEmil(F_1,nSites,nPop);
}

int emAccelUnadjusted(const bgl &d, int nPop, double **F,double *Q,double *Q_new,double &lold, accFQ &a){
  //maybe these should be usersettable?
  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  double objfnInc=1;
  // first EM run
  emUnadjusted(Q, F, d.nSites, nPop,d.genos, a.Q_em1);
  // diff between 0 and 1 assigned to diff1
  minus1dEmil(a.Q_em1,Q,nPop,a.Q_diff1);
  // calculates "norm" of F_diff1 and Q_diff1 and sums them
  double sr2 = sumSquare1dEmil(a.Q_diff1,nPop);
  // checks if convergence
  if(sqrt(sr2)<tol){
    return 0;
  }
  // second EM run
  emUnadjusted(a.Q_em1, F, d.nSites, nPop,d.genos, a.Q_em2);
  // diff between 1 and 2 assigned to diff2  
  minus1dEmil(a.Q_em2,a.Q_em1,nPop,a.Q_diff2);
  // calculates "norm" of F_diff2 and Q_diff2 and sums them 
  double sq2 = sumSquare1dEmil(a.Q_diff2,nPop);
  // checks if convergence - a second time
  if(sqrt(sq2)<tol){
    return 0;
  }
  // diff between diff1 and diff2
  minus1dEmil(a.Q_diff2,a.Q_diff1,nPop,a.Q_diff3);
  // WHY NOT SQ2 - because it is v = diff2 - diff1 actually
  double sv2 = sumSquare1dEmil(a.Q_diff3,nPop);
  // does not have to make alpha negative because changed accordingly in equations
  double alpha = sqrt(sr2/sv2);
  // makes sure alpha does not go below 1 and above stepMax
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  // update with the linear combination
  // updates F with new values via alpha - where alpha controls rate of change
  // check that F has no entries of 0
  for(size_t i=0;i<nPop;i++){
    // because alpha pos in code, does not need minus -(-alpha) = alpha
    Q_new[i] = Q[i]+2*alpha*a.Q_diff1[i]+alpha*alpha*a.Q_diff3[i];
    a.Q_tmp[i] = 1.0;
  }
  double* tmpAdrQ=a.Q_tmp;
  map2domainQEmil(Q_new,nPop);
  // if alpha not too close (0.01 close) to 1 
  if (fabs(alpha - 1) > 0.01){
    // we estimate new Q and F, with our inferred Q and F via alpha
    emUnadjusted(Q_new, F, d.nSites, nPop,d.genos,a.Q_tmp);
    // assign value of Q_tmp and F_tmp to Q_new and F_new
    std::swap(Q_new,a.Q_tmp);
    // in order to avoid a.Q_tmp or a.F_tmp pointing to Q and F
    a.Q_tmp=tmpAdrQ;
  }
  double lnew =0;
  if (alpha == stepMax) {
    stepMax = mstep*stepMax;
  }
  lold=lnew;
  return 1;
}

//returnval =1 continue, rturnval =0, convergence has been achieved
// lold in this function is a reference to the argument passed on to it, so every change in this function is also applied to function outside
int emAccelEmil(const bgl &d, double* nInd, int nPop, double **F,double *Q,double **F_new,double *Q_new,double &lold, double** F_org, accFQ &a){
  //maybe these should be usersettable?
  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  double objfnInc=1;
  // first EM run
  emEmil(Q, F, d.nSites, nInd, nPop,d.genos, a.F_em1, a.Q_em1, F_org);
  // diff between 0 and 1 assigned to diff1
  minusEmil(a.F_em1,F,d.nSites,nPop,a.F_diff1);
  minus1dEmil(a.Q_em1,Q,nPop,a.Q_diff1);
  // calculates "norm" of F_diff1 and Q_diff1 and sums them
  //  double sr2 = sumSquareEmil(a.F_diff1,d.nSites,nPop)+sumSquare1dEmil(a.Q_diff1,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
  double sr2 = sumSquare1dEmil(a.Q_diff1,nPop);
  // checks if convergence
  if(sqrt(sr2)<tol){
    return 0;
    //break;
  }
  // second EM run
  emEmil(a.Q_em1, a.F_em1, d.nSites, nInd, nPop,d.genos, a.F_em2, a.Q_em2, F_org);
  // diff between 1 and 2 assigned to diff2  
  minusEmil(a.F_em2,a.F_em1,d.nSites,nPop,a.F_diff2);
  minus1dEmil(a.Q_em2,a.Q_em1,nPop,a.Q_diff2);
  // calculates "norm" of F_diff2 and Q_diff2 and sums them 
  //double sq2 = sumSquareEmil(a.F_diff2,d.nSites,nPop)+sumSquare1dEmil(a.Q_diff2,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
  double sq2 = sumSquare1dEmil(a.Q_diff2,nPop);
  // checks if convergence - a second time
  if(sqrt(sq2)<tol){
    return 0;
    //    break;
  }
  // diff between diff1 and diff2
  minusEmil(a.F_diff2,a.F_diff1,d.nSites,nPop,a.F_diff3);
  minus1dEmil(a.Q_diff2,a.Q_diff1,nPop,a.Q_diff3);
  // WHY NOT SQ2 - because it is v = diff2 - diff1 actually
  //double sv2 = sumSquareEmil(a.F_diff3,d.nSites,nPop)+sumSquare1dEmil(a.Q_diff3,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
  double sv2 = sumSquare1dEmil(a.Q_diff3,nPop);
  // does not have to make alpha negative because changed accordingly in equations
  double alpha = sqrt(sr2/sv2);
  // makes sure alpha does not go below 1 and above stepMax
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  // updates F with new values via alpha - where alpha controls rate of change
  for(size_t i=0;i<d.nSites;i++){
    for(size_t j=0;j<nPop;j++){
      // (*F_new)[i][j] wtf??
      F_new[i][j] = F[i][j]+2*alpha*a.F_diff1[i][j]+alpha*alpha*a.F_diff3[i][j];
      a.F_tmp[i][j] = 1.0;
    }
  }
  // check that F has no entries of 0
  map2domainFEmil(F_new,d.nSites,nPop);
  for(size_t i=0;i<nPop;i++){
    // because alpha pos in code, does not need minus -(-alpha) = alpha
    Q_new[i] = Q[i]+2*alpha*a.Q_diff1[i]+alpha*alpha*a.Q_diff3[i];
    a.Q_tmp[i] = 1.0;
  }
  double** tmpAdrF=a.F_tmp;
  double* tmpAdrQ=a.Q_tmp;
  map2domainQEmil(Q_new,nPop);
  // if alpha not too close (0.01 close) to 1 
  if (fabs(alpha - 1) > 0.01){
    // we estimate new Q and F, with our inferred Q and F via alpha
    emEmil(Q_new, F_new, d.nSites, nInd, nPop,d.genos,a.F_tmp,a.Q_tmp,F_org);
    // assign value of Q_tmp and F_tmp to Q_new and F_new
    std::swap(Q_new,a.Q_tmp);
    std::swap(F_new,a.F_tmp);
    // in order to avoid a.Q_tmp or a.F_tmp pointing to Q and F
    a.Q_tmp=tmpAdrQ;
    a.F_tmp=tmpAdrF;
  }
  double lnew =0;
  if (alpha == stepMax) {
    stepMax = mstep*stepMax;
  }
  lold=lnew;
  return 1;
}


void info(){
  
  fprintf(stderr,"Arguments:\n");
  fprintf(stderr,"\t-likes Beagle likelihood filename\n");
  fprintf(stderr,"\t-K Number of ancestral populations\n"); 
  fprintf(stderr,"\t-Nname Number of individuals in each reference populations\n"); 
  fprintf(stderr,"Optional:\n");
  fprintf(stderr,"\t-fname Ancestral population frequencies\n"); 
  fprintf(stderr,"\t-qname Admixture proportions of each ancestral population\n"); 
  fprintf(stderr,"\t-outfiles Prefix for output files\n"); 
  fprintf(stderr,"\t-printInfo print ID and mean maf for the SNPs that were analysed\n"); 
  fprintf(stderr,"\t-printFreq print admixture adjusted allele frequencies of reference panel + input individual (1: yes, 0: no (default))\n"); 

  fprintf(stderr,"Setup:\n"); 
  fprintf(stderr,"\t-doAdjust Adjusts the frequencies in the reference populations with the input (1: yes (default), 0: no)\n");
  fprintf(stderr,"\t-seed Seed for initial guess in EM\n"); 
  fprintf(stderr,"\t-P Number of threads\n"); 
  fprintf(stderr,"\t-method If 0 no acceleration of EM algorithm (1: yes (default), 0: no)\n"); 
  fprintf(stderr,"\t-misTol Tolerance for considering site as missing\n");

  fprintf(stderr,"Stop chriteria:\n"); 
  fprintf(stderr,"\t-Qconv Stopping criteria based on change in Q (works best when using doAdjust) (1: yes, 0: no (default))\n"); 
  fprintf(stderr,"\t-Qtol Tolerance value for stopping criteria based on change in Q (0.001 (default))\n"); 
  fprintf(stderr,"\t-tolLike50 Loglikelihood difference in 50 iterations\n"); 
  fprintf(stderr,"\t-tol Tolerance for convergence\n"); 
  fprintf(stderr,"\t-dymBound Use dymamic boundaries (1: yes (default) 0: no)\n"); 
  fprintf(stderr,"\t-maxiter Maximum number of EM iterations\n"); 
  fprintf(stderr,"\t-boot Number of bootstrapping iterations, default 50\n"); 

  fprintf(stderr,"Filtering\n"); 
  fprintf(stderr,"\t-minMaf Minimum minor allele frequency - does not really work!\n"); 


  exit(0);
}



int VERBOSE =1;
void handler(int s) {
  if(VERBOSE)
    fprintf(stderr,"Caught SIGNAL: %d will try to exit nicely (no more threads are created, we will wait for the current threads to finish)\n",s);
  VERBOSE=0;
  SIG_COND=0;
}

float calcThresEmil(double *d1,double *d2, int x){
  // finds the largest difference between 2 arrays
  // arrays has dimention x times y
  float diff=0;
  for(int i=0;i<x;i++){
    if(fabs(d1[i]-d2[i])>diff){
      diff=fabs(d1[i]-d2[i]);
    }
  }
  return diff;
}

void printLikes(bgl &d){
  for(int s=0;s<d.nSites;s++){
    double *g  = d.genos[s];
    for(int i=0;i<d.nInd;i++){
      for(int j=0;j<3;j++)
	fprintf(stdout,"%f\t",g[3*i+j]);
    }
    fprintf(stdout,"\n");
  }
  exit(0);
}

void printKeepSites(bgl &d,FILE *ffilter){
  fprintf(ffilter,"marker\tmajor\tminor\tmaf\tnonMis\n");
 for(int s=0;s<d.nSites;s++){
   fprintf(ffilter,"%s\t%d\t%d\t%f\t%d\n",d.ids[s],d.major[s],d.minor[s],d.mafs[s],d.keepInd[s]);

 }
}


void filterMinMaf(bgl &d,float minMaf, int* &badMaf){
  int posi = 0;
  for(int s=0;s<d.nSites;s++){
    if(d.mafs[s]>minMaf&&d.mafs[s]<1-minMaf){
      d.genos[posi] = d.genos[s];
      d.major[posi] = d.major[s];
      d.minor[posi] = d.minor[s];
      d.ids[posi] = d.ids[s];
      d.keeps[posi] = d.keeps[s];
      d.keepInd[posi] = d.keepInd[s];
      d.mafs[posi] = d.mafs[s];
      // store which one filtered away
      badMaf[s] = 1;
      posi++;   
    } 
    
  }
  d.nSites=posi;
}


  void filterMinLrt(bgl &d,float minLrt,int* &badLrt){
  int posi =0;
  for(int s=0;s<d.nSites;s++){
    float lik=likeFixedMinor(d.mafs[s],d.genos[s],d.nInd,d.keeps[s]);
    float lik0=likeFixedMinor(0.0,d.genos[s],d.nInd,d.keeps[s]);
    if(2.0*(lik0-lik)>minLrt){
      d.genos[posi] = d.genos[s];
      d.major[posi] = d.major[s];
      d.minor[posi] = d.minor[s];
      d.ids[posi] = d.ids[s];
      d.keeps[posi] = d.keeps[s];
      d.keepInd[posi] = d.keepInd[s];
      d.mafs[posi] = d.mafs[s];
      // store which one filtered away
      badLrt[s] = 1;
      posi++;   
    } 
  }
  d.nSites=posi;
}

////////////////////////// it begins 
 int main(int argc, char **argv){ 
  if(argc==1){// if no arguments, print info on program
    info();
    return 0;
  }
  //below for catching ctrl+c, and dumping files
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

  //initial values
  int dymBound = 0;
  int maxIter = 2000;
  int method = 1;
  int printInfo = 0;
  int printFreq = 0;
  int doAdjust = 1;
  float minMaf =0.00;
  const char* lname = NULL;
  // is this an issue that it is not const??
  char* plinkName = NULL;
  const char* fname = NULL;
  const char* qname = NULL;
  const char* Nname = NULL;
  const char* outfiles = NULL;
  int nPop = 3;
  int seed = time(NULL);
  // float tolLike50=0.1;
  float tolLike50=0.01; // changed by Emil
  int nBoot = 50; // for bootstrapping
  int Qconv = 0;
  double Qtol = 0.001;
  // reading arguments
  argv++;
  while(*argv){
    // GL in the shape of beagle file
    if(strcmp(*argv,"-likes")==0 || strcmp(*argv,"-l")==0) lname=*++argv; //name / char arrays
    else if(strcmp(*argv,"-plink")==0 || strcmp(*argv,"-p")==0) plinkName=*++argv; 
    // probably will need this, implicitly in the freqs file, by default 3
    else if(strcmp(*argv,"-K")==0) nPop=atoi(*++argv); 
    // would also need number of individuals in each ref category
    // to read start values from output from previous run 
    else if(strcmp(*argv,"-fname")==0 || strcmp(*argv,"-f")==0) fname=*++argv; 
    // starting guess not really sure we need this
    else if(strcmp(*argv,"-qname")==0 || strcmp(*argv,"-q")==0) qname=*++argv;
    // for reading in number of individauls in each population
    else if(strcmp(*argv,"-Nname")==0 || strcmp(*argv,"-N")==0) Nname=*++argv;
    // prefix for output files
    else if(strcmp(*argv,"-outfiles")==0 || strcmp(*argv,"-o")==0) outfiles=*++argv; 
    // settings: seed, threads and if method==0 not accelerated
    else if(strcmp(*argv,"-seed")==0||strcmp(*argv,"-s")==0) seed=atoi(*++argv); //int - atoi - char array to integer
    // flag for printing adjusted freqs
    else if(strcmp(*argv,"-printFreq")==0) printFreq=atoi(*++argv); 
    else if(strcmp(*argv,"-printInfo")==0) printInfo=atoi(*++argv); 
    // flag for doing Adjustment
    else if(strcmp(*argv,"-doAdjust")==0) doAdjust=atoi(*++argv);
    else if(strcmp(*argv,"-method")==0 || strcmp(*argv,"-m")==0) method=atoi(*++argv); 
    // different stop criteria - whether based on diff in Q values
    else if(strcmp(*argv,"-Qconv")==0) Qconv=atoi(*++argv); 
    // tolerance for Q stopping criteria
    else if(strcmp(*argv,"-Qtol")==0) Qtol=atof(*++argv); 
    else if(strcmp(*argv,"-tolLike50")==0||strcmp(*argv,"-lt50")==0) tolLike50=atof(*++argv); //float/double - atof - char array to double/float
    // do I need those when I have only one individual??
    else if(strcmp(*argv,"-bootstrap")==0||strcmp(*argv,"-boot")==0) nBoot=atoi(*++argv);
    else if(strcmp(*argv,"-tol")==0||strcmp(*argv,"-t")==0) tol=atof(*++argv);
    else if(strcmp(*argv,"-maxiter")==0 || strcmp(*argv,"-i")==0) maxIter=atoi(*++argv); 
    // different filterings 
    else if(strcmp(*argv,"-misTol")==0 || strcmp(*argv,"-mt")==0) misTol=atof(*++argv);
    // might not really work, not sure if should apply to refPanel?
    else if(strcmp(*argv,"-minMaf")==0||strcmp(*argv,"-maf")==0) minMaf=atof(*++argv);

    // different genotype callers - tolerance which we do use
    else if(strcmp(*argv,"-dymBound")==0) dymBound=atoi(*++argv);
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      info();
      return 0;
    }
    ++argv;
  }


  //check that non optional options have been used. 
  if(lname==NULL and plinkName==NULL){
    fprintf(stderr,"Please supply either beagle or plink input file: -likes or -plink");
    info();
  }
  if(outfiles==NULL){
    fprintf(stderr,"Will use beagle/plink fname as prefix for output\n");
    outfiles=lname;
  }

  //out put files
  FILE *flog=openFile(outfiles,".log");
  FILE *ffilter=openFile(outfiles,".filter");

  fprintf(stderr,"Input: lname=%s nPop=%d, fname=%s qname=%s outfiles=%s\n",lname,nPop,fname,qname,outfiles);
  fprintf(stderr,"Setup: seed=%d method=%d\n",seed,method);
  if(method==0){
    fprintf(stderr,"The unaccelerated EM has been chosen\n");
  } else{
    fprintf(stderr,"The acceleted EM has been chosen\n");
  }
  if(doAdjust==0){
    fprintf(stderr,"The unadjusted method has been chosen\n");
  } else{
    fprintf(stderr,"The adjusted method has been chosen\n");
  }
  fprintf(stderr,"Convergence: maxIter=%d tol=%f tolLike50=%f dymBound=%d\n",maxIter,tol,tolLike50,dymBound);
  fprintf(stderr,"The following number of bootstraps have been chosen: %i\n",nBoot);
  if(Qconv>0){
    fprintf(stderr,"Convergence via difference in Q values chosen, threshold of: %f\n",Qtol);    
  }
  fprintf(stderr,"Filters: misTol=%f minMaf=%f\n",misTol,minMaf);

  fprintf(flog,"Input: lname=%s nPop=%d, fname=%s qname=%s outfiles=%s\n",lname,nPop,fname,qname,outfiles);
  fprintf(flog,"Setup: seed=%d method=%d\n",seed,method);
  if(method==0){
    fprintf(flog,"The unacceleted EM has been chosen\n");
  } else{
    fprintf(flog,"The acceleted EM has been chosen\n");
  }
  if(doAdjust==0){
    fprintf(flog,"The unadjusted method has been chosen\n");
  } else{
    fprintf(flog,"The adjusted method has been chosen\n");
  }
  fprintf(flog,"Convergence: maxIter=%d tol=%f tolLike50=%f dymBound=%d\n",maxIter,tol,tolLike50,dymBound);
  fprintf(flog,"The following number of bootstraps have been chosen: %i\n",nBoot);
  if(Qconv>0){
    fprintf(flog,"Convergence via difference in Q values chosen, threshold of: %f\n",Qtol);    
  }
  fprintf(flog,"Convergence: maxIter=%d tol=%f tolLike50=%f dymBound=%d\n",maxIter,tol,tolLike50,dymBound);
  fprintf(flog,"Filters: misTol=%f minMaf=%f\n",misTol,minMaf);

  if(dymBound==0){
    errTolStart = errTolMin;
    errTol = errTolMin;
  }
    
  clock_t t=clock();//how long time does the run take
  time_t t2=time(NULL);

  //read BEAGLE likelihood file  
  // made into object to give it additional info
  bgl d;
  bgl dOrg;

  if(lname!=NULL){
    d=readBeagle(lname);
    dOrg=readBeagle(lname);
  } else {
    info();
  }
   
  fprintf(stderr,"Input file has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);
  fprintf(flog,"Input file has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);

  // START HERE
  // in order to have same sites in freq file
  // the genotypes should be read in and converted to GL
  // see README for details on how to do this
  // should be flipped in Rscript
  // it seems ADMIXTURE output - has diff direction of freqs than plink
  // do some function that reads plink and converts to GL
 


  // filter sites based on MAF - 

  if(printInfo){
    printKeepSites(d,ffilter);
  }
  fflush(stderr);

  // seed for bootstrapping, have only one seed!
  std::srand(seed);
  accFQ a = createAccFQ(nPop, d.nSites);
 
  //unknown parameters
  // we only have a vector of Qs and then a matrix of freqs
  // so allocates a an array of pointers, each pointing to a second array
  // this gives matrix like structure of F
  double **F = allocDouble(d.nSites,nPop);
  double **F_new = allocDouble(d.nSites,nPop);
  // F_org is for storing initial freqs from ref panel
  double **F_org = allocDouble(d.nSites,nPop);
  // for sampling from when doing bootstrap
  double **F_orgOrg = allocDouble(d.nSites,nPop);
  double **F_1stRun = allocDouble(d.nSites,nPop);  

  //get start values
  readDoubleGZ(F,d.nSites,nPop,fname,0);
  
  for(int i=0;i<d.nSites;i++){
    for(int k=0;k<nPop;k++){
      F_org[i][k]= F[i][k];
      F_orgOrg[i][k]= F[i][k];
    }
     
  }

  // because has to have value for normal data
  // and then nBoot bootstrapped values
  double **Q = allocDouble(nBoot+1,nPop);
  double **Q_new = allocDouble(nBoot+1,nPop);

  double *N = new double[nPop];

  // putting in an intial guess 
  double *sum = new double[nBoot+1];
  // nBoot + 1 rows
  for(int j=0;j<=nBoot;j++){
    sum[j]=0;
    for(int k=0;k<nPop;k++){
      Q[j][k]=1.0/nPop;
      sum[j]+=Q[j][k];
    }
  }
  // nBoot + 1 rows
  for(int j=0;j<=nBoot;j++){
    for(int k=0;k<nPop;k++) {
      // to make sure that proportions sum to 1
      Q[j][k] = Q[j][k]/sum[j]; 
      Q_new[j][k] = Q[j][k];
    }
  }
  delete [] sum;

  if(Nname==NULL){
    fprintf(stderr,"Please supply number of individauls file: -Nname");
    info();
  }
  // reading nInd
  readDouble1d(N,nPop,Nname);

   
  //update the global stuff NOW
  //update the internal stuff in the pars for the threading
  double lold = likelihoodEmil(Q[0], F, d.nSites, nPop,d.genos);
  
  fprintf(stderr,"iter[start] like is=%f\n",lold);

  //////////////////////////////////////// em ///////////////////////////////////  
  //below is the main looping trhought the iterations.
  // we have 4 possible ways, unAdjusted/Adjusted basic EM/Accelerated EM
  int nit;
  double likeLast = lold;
  double lastQthres = 0;


  for(int b=0;SIG_COND and b<=nBoot;b++) { 
    for(nit=1;SIG_COND and nit<maxIter;nit++) { 
      if(doAdjust==0){
	if(method==0){
	  emUnadjusted(Q[b], F, d.nSites, nPop,d.genos,Q_new[b]);
	} else{
	  if(emAccelUnadjusted(d, nPop, F, Q[b], Q_new[b],lold, a)==0){    
	    if(errTol>errTolMin){
	      errTol=errTol/5;
	      if(errTol<errTolMin){
		errTol=errTolMin;
	      }
	    }
	  }
	}
      } else{   
	if(method==0){
	  emEmil(Q[b], F, d.nSites, N, nPop,d.genos,F_new,Q_new[b],F_org);
	} else{
	  
	  if(emAccelEmil(d, N, nPop, F, Q[b], F_new, Q_new[b],lold,F_org, a)==0){
	    if(errTol>errTolMin){
	      errTol=errTol/5;
	      if(errTol<errTolMin){
		errTol=errTolMin;
	      }
	    }
	    
	    else{
	      if(b==0){
		fprintf(stderr,"EM accelerated has reached convergence with tol %f\n",tol);
	      }
	      break; //if we have achieved convergence
	    }
	  }
	}
	std::swap(F,F_new);
      }
      std::swap(Q,Q_new);
      if((nit%10)==0 ){ //stopping criteria
	double lik = likelihoodEmil(Q[b], F, d.nSites, nPop,d.genos);
	// thres is largest differense in admixture fractions
	if(b==0){
	  fprintf(stderr,"iter[%d] like is=%f thres=%f\n",nit,lik,calcThresEmil(Q[b],Q_new[b],nPop));
	  fprintf(stderr,"iter[%d] diff in likelihood is=%f\t",nit,std::abs(lik-likeLast));      
	  fprintf(stderr,"iter[%d] Q is=%f, %f, %f\t",nit,Q[b][0],Q[b][1],Q[b][2]);      
	}
	if(errTol>errTolMin){
	  errTol=errTol/10;
	  if(errTol<errTolMin)
	    errTol=errTolMin;
	} else if(calcThresEmil(Q[b],Q_new[b],nPop) < Qtol and Qconv>0) {
	  if(b==0){
	    fprintf(stderr,"Convergence achived because diffence in Q values less than %f\n",Qtol);
	  }
	  break;
	} else if(fabs(calcThresEmil(Q[b],Q_new[b],nPop)-lastQthres)<1e-07){
	  if(b==0){
	    fprintf(stderr,"Convergence achived because diffence in Q values less than %f\n",Qtol);
	  }
	  break;
	}
	
	else if(std::abs(lik-likeLast) < tolLike50 and Qconv==0) {
	  if(b==0){
	    fprintf(stderr,"Convergence achived becuase log likelihooditer difference for 50 iteraction is less than %f\n",tolLike50);
	  }
	  if(lik+likeLast<-1){
	    if(b==0){
	      fprintf(stderr,"Convergence achived because log likelihooditer difference was NEGATIVE\n");
	    }
	  }
	  break;
	}
	likeLast=lik; 
	lastQthres=calcThresEmil(Q[b],Q_new[b],nPop);
      } 
    }
    // first run
    if(b==0){
      // nBoot + 1 rows
      for(int j=1;j<=nBoot;j++){
	for(int k=0;k<nPop;k++){
	  //	  make first estimated Q starting guess for rest
	  Q[j][k] = Q[0][k];
	  Q_new[j][k] = Q[j][k];
	  
	}
      }
      // same for GL input!!
      std::swap(F,F_1stRun);
    }
    // F_orgOrg original F sampled from
    bootstrap(dOrg,d,F_orgOrg,F_org,F,nPop); 

    if(b>0){
      fprintf(stderr,"At this bootstrapping: %i out of: %i\n",b,nBoot);
    }
  }
  lold = likelihoodEmil(Q[0], F, d.nSites, nPop,d.genos);
  fprintf(stderr,"best like=%f after %d iterations\n",lold,nit);
  fprintf(flog,"estimated  Q = ");
  for(int i=0;i<nPop;i++){
    fprintf(flog,"%f ",Q[0][i]);
  }
  fprintf(flog,"\n");
  fprintf(flog,"best like=%f after %d iterations\n",lold,nit);
  
  /////////////////////////////////////////////////////////// done - make output and clean /////////////////////////////////  
  // Print F and Q in files

  // change to mean CIlower CIupper
  FILE *fp=openFile(outfiles,".qopt");
  // nBoot + 1, because first value is estimated Q
  printDouble(Q,nBoot+1,nPop,fp);
  fclose(fp);

  // only if certain flag
  if(printFreq>0){
    gzFile fpGz=openFileGz(outfiles,".fopt.gz");
    printDoubleGz(F_1stRun,d.nSites,nPop,fpGz);
    gzclose(fpGz);
  }
  
  //deallocate memory 
  dalloc(F,d.nSites); 
  dalloc(F_1stRun,d.nSites);  
  dalloc(F_new,d.nSites);  
  dalloc(Q,nBoot); 
  dalloc(Q_new,nBoot); 

  dalloc(F_org,d.nSites);
  dalloc(F_orgOrg,d.nSites);

  delete [] N;

  // problem with F and F_tmp being same pointer 
  // not there any more really

  dallocAccFQ(a,d.nSites);
  dallocBeagle(d);
  dallocBeagle(dOrg);

  for(int i=0;1&&i<dumpedFiles.size();i++){
    //    fprintf(stderr,"dumpedfiles are: %s\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  

  // print to log file
  fprintf(flog, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(flog, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fclose(flog); 
  
  fclose(ffilter);
  return 0;
    
 }

