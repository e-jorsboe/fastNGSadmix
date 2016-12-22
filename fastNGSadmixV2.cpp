/*
  log:
  g++ fastNGSadmixV2.cpp -lz -lpthread  -O3 -o fastNGSadmixV2

  log: (with readplink function)
  g++ fastNGSadmixV2.cpp readplinkV2.c -lz -lpthread  -O3 -o fastNGSadmixV2


  debug:
  g++ fastNGSadmixV2.cpp -lz -lpthread -ggdb -O3 -o fastNGSadmixV2

  debug: (with readplink function)
  g++ fastNGSadmixV2.cpp readplinkV2.c -lz -lpthread -ggdb -O3 -o fastNGSadmixV2

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
#include <map>


// not really necessary any longer
#include "readplinkV2.h"

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
double tol=1e-7; //stopping criteria

double maxAlpha = 30.0;
double errTolMin=1e-5;
double errTolStart=0.05;
double errTol=errTolStart;//frequencies and admixture coef cannot be less than this or more than 1-this
double misTol=0.05;



double **allocDouble(size_t x,size_t y){
  double **ret = new double*[x];
  for(size_t i=0;i<x;i++){
    ret[i] = new double[y];
  }
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

// function for keeping sure Q values do not become
// 0.0 as then division by 0 might occur, errTol is limit
void map2domainQEmil(double* &Q, int nPop){  
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

// function for keeping sure F values do not become
// 0.0 as then division by 0 might occur, errTol is limit
void map2domainFEmil(double** &F, int nSites, int nPop){
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

struct cmp_char {
  bool operator()(const char *first,const char* second) const {
    int tmp = std::strcmp(first, second);
    return tmp<0;
  }
};




//some struct will all the data from the beagle file
typedef struct{
  double **genos;
  char *major;
  char *minor;
  // for snp names pointer to pointer
  char **ids;
  int nSites;
  int nInd;
  //matrix containing 0/1 indicating if data or missing
  char **keeps;
  //keepInd[nSites] this is the number if informative samples
  int *keepInd; 
  float *mafs;
  // map of ids in beagle file for finding overlap with ref
  // id is like this: chr_pos 
  std::map <char*,int,cmp_char> idMap;
  
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
    b.keeps[s] = "0";
    b.keepInd[s] = 0;
    b.mafs[s] = 0.0;
  
  }
  return(b);
} 
      
//utility function for cleaning up out datastruct
void dallocBeagle(bgl &b){
  for(int i=0;i<b.nSites;i++){
    delete [] b.genos[i];
    
    free(b.ids[i]);
    
    //    delete [] b.ids[i];
  }
  

 
  // run through all the keys and then free them

  
  delete [] b.minor;
  delete [] b.major;
  delete [] b.genos;
  delete [] b.ids;
  delete [] b.keeps;

  delete [] b.keepInd;
  delete [] b.mafs;

  

}

// refPanel struct for reading in refPanel with header
typedef struct{
  char **id;
  int *chr;
  int *pos;
  char **name;
  char *A0;
  char *A1;
  double **freqs;
  int refSites;
  char **populations;
  // has map of column to keep in ref for calculations
  // is coded so key is old column number
  // value is new column number - 1 (cause has to be above 0)
  std::map <int,int> colsToKeep;
  std::map <char *,int,cmp_char> popsToKeep;
}refPanel;

refPanel allocRefPanel(int nSites, int refPops, int nPop){
  refPanel ref;

  ref.refSites = nSites;

  ref.id = new char*[nSites];
  ref.chr = new int[nSites];
  ref.pos = new int[nSites];

  ref.name = new char*[nSites];
  ref.populations = new char*[nPop];
  ref.A0 = new char[nSites];
  ref.A1 = new char[nSites];
  ref.freqs= allocDouble(nSites,refPops);

  for(int s=0;SIG_COND&& (s<nSites);s++){

    ref.id[s] = strdup("0");
    ref.chr[s] = 0;
    ref.pos[s] = 0;
    ref.name[s] = strdup("0");
    
    ref.A0[s] = '0';
    ref.A1[s] = '0';
  
    for(int i=0;i<refPops;i++){
      ref.freqs[s][i] = 0.0;
      ref.populations[i] = strdup("0");
    }
  }
  return(ref);
} 

void dallocRefPanel(refPanel &ref, int nPop){
  for(int i=0;i<ref.refSites;i++){
    delete [] ref.freqs[i];
    //    free(ref.id[i]);
    free(ref.name[i]);
    
    
  }
  for(int j=0;j<nPop;j++){
    free(ref.populations[j]);
  }
  delete [] ref.chr;
  delete [] ref.pos;
  delete [] ref.A0;
  delete [] ref.A1;
  delete [] ref.freqs;
  delete [] ref.name;
  delete [] ref.populations;
  delete [] ref.id;
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

// convert 0,1,2,3 to A,C,G,T if beagle coded thusly
char intToChar(char intLike){
  if(intLike=='A' or intLike=='C' or intLike=='G' or intLike=='T'){
    return(intLike);
  } else {
  if(intLike == '0'){
    return('A');
  } else if(intLike == '1'){
    return('C');
  } else if(intLike == '2'){
    return('G');
  } else if(intLike == '3'){
    return('T');
  } else{
    fprintf(stderr,"Beagle nucleotide not valid must be A,C,G,T or 0,1,2,3\n");
    exit(0);
  }
  }

}

bgl readBeagle(const char* fname, std::map <char*,int,cmp_char> overlap) {
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
  while(strtok(NULL,delims)){
    ncols++;
  }
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
    // puts id of all sites in map for fast lookup
    //    ret.idMap.insert(std::pair<char*,int>(strtok(buf,delims),1));

    ret.idMap[strtok(strdup(buf),delims)] = 1;


    
  }

  
  //now we now the number of sites
  ret.nSites = overlap.size();
  ret.major = new char[overlap.size()];
  ret.minor = new char[overlap.size()];
  ret.ids = new char*[overlap.size()];
  ret.genos= new double*[overlap.size()];
  //then loop over the vector and parsing every line
  int  refIndex = 0;
  for(int s=0;SIG_COND&& (s<tmp.size());s++){
    
    char * id = strtok(tmp[s],delims);

    if(overlap.count(id)==0){
      continue;
    } 
    
    ret.ids[refIndex] = strdup(id);
    // because allele might be coded 0,1,2,3
    ret.major[refIndex] = intToChar(strtok(NULL,delims)[0]);
    ret.minor[refIndex] = intToChar(strtok(NULL,delims)[0]);
    ret.genos[refIndex] = new double[3*ret.nInd];
    for(int i=0;i<ret.nInd*3;i++){
      ret.genos[refIndex][i] = atof(strtok(NULL,delims));
      if(ret.genos[refIndex][i]<0){
	fprintf(stderr,"Likelihoods must be positive\n");
	fprintf(stderr,"site %d ind %d geno %d has value %f\n",s,int(i*1.0/3),i%3,ret.genos[refIndex][i]);
	exit(0);
      }
    }
    for(int i=0;i<ret.nInd;i++){
      double tmpS = 0.0;
      for(int g=0;g<3;g++)
	tmpS += ret.genos[refIndex][i*3+g];
      if(!(tmpS>0)){
	fprintf(stderr,"The sum of likelihoods for a genotypes must be positive\n");
	fprintf(stderr,"individual %d site %d has sum %f\n",i,s,tmpS);
	exit(0);
      } 
    }
    refIndex++;
    
  }
  // Here additional stuff is calculated from the likelihoods
  // this must be done again while filtering later in main
  ret.keeps = new char*[overlap.size()]; // array nSites x nInd 0 if missing info
  ret.keepInd = new int[overlap.size()];
  ret.mafs = new float[overlap.size()];
  for(int s=0;s<overlap.size();s++){
    ret.keeps[s] = new char[ret.nInd];
    int nKeep =0;
    for(int i=0;i<ret.nInd;i++){
      double mmin=std::min(ret.genos[s][i*3],std::min(ret.genos[s][i*3+1],ret.genos[s][i*3+2]));
      double mmax=std::max(ret.genos[s][i*3],std::max(ret.genos[s][i*3+1],ret.genos[s][i*3+2]));
      if(fabs(mmax-mmin)<misTol)
	ret.keeps[s][i] = 0;
      else{
	ret.keeps[s][i] = 1;
	nKeep++;
      }
    }
    ret.keepInd[s] = nKeep;
    ret.mafs[s] = emFrequency(ret.genos[s],ret.nInd,MAF_ITER,MAF_START,ret.keeps[s],ret.keepInd[s]);
  }
  for(int s=0;s<tmp.size();s++){
    free(tmp[s]);
  }
  gzclose(fp); //clean up filepointer
  return ret;
}




bgl readPlinkToBeagle(const char* plinkName, std::map <char*,int,cmp_char> overlap) {
  // denoting delimeters in file
  // annoying that pointer to struct
  plink* pl = readplink(plinkName);
  double seqError = 0.00;

  // make readBim, where gets id and alleles 
  // I do not think it gets that the beagle file has been create
  // and therefore the d.genos is null
  bgl b = allocBeagle(overlap.size());
  b.nSites = overlap.size();
  b.nInd = 1; // can only have one individual is this program

  int beagleIndex = 0;
  for(int i=0;i<(pl->y);i++){
  
  // put second plink allele first in bgl
  // and then first plink allele
  // then put id as chr_pos
 
    if(overlap.count(pl->bim.id[i])<1){

      continue;
    }
    if(beagleIndex+1 > overlap.size()){
      fprintf(stderr,"Duplicated sites - multiple sites with same position in plink file!! \n");
      fprintf(stderr,"Use plink2 with --list-duplicate-vars suppress-first ids-only - and then --exclude \n");
      exit(0);
	
    }
  
    if(pl->d[0][i]==0){
      b.genos[beagleIndex][2]=1.0-seqError;
      b.genos[beagleIndex][1]=seqError/2.0;
      b.genos[beagleIndex][0]=seqError/2.0;
      
    } else if(pl->d[0][i]==1){
      b.genos[beagleIndex][2]=seqError/2.0;
      b.genos[beagleIndex][0]=seqError/2.0;
      b.genos[beagleIndex][1]=1.0-seqError;
    } else if(pl->d[0][i]==2){
      b.genos[beagleIndex][2]=seqError/2.0;
      b.genos[beagleIndex][1]=seqError/2.0;
      b.genos[beagleIndex][0]=1.0-seqError;
    }
    b.major[beagleIndex]=pl->bim.major[i];
    b.minor[beagleIndex]=pl->bim.minor[i];
    b.ids[beagleIndex]=strdup(pl->bim.id[i]);
    b.idMap[strdup(pl->bim.id[i])] = 1;

    beagleIndex++;
  }
  kill_plink(pl);
  return(b);
  

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


// function for reading in ref panel, has to have beagle file for finding overlapping sites
// has to have included pops to be able to extract those columns and K
// has to have format id chr pos name A0_freq A1 pop1 pop2 ... (freq has to be of A0 allele)
// Beagle files assumed to be subset of ref (all sites in beagle also in ref)
std::vector<char *> tmp;
refPanel readRefPanel(const char* fname, bgl &b, std::map <char*,int,cmp_char> &includedPops, int nPop, std::map <char*,int,cmp_char> overlap) {
  // denoting delimeters in file
  const char *delims = "\t \n";
  gzFile fp = NULL;
  // checking if file can be opened
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  // creating ref struct
  refPanel ref;
  // creaing string of length LENS
  char buf[LENS];
  //find number of columns
  // reading string from compressed file 
  
  
  // while still some left of string from file
  while(gzgets(fp,buf,LENS)){
    // duplicates string and puts it into tmp vector
    tmp.push_back(strdup(buf));    
  }
  // copy first line into buffer, as empty now on some systems
  // Copies the values of num bytes from the location pointed to by source directly to the memory block pointed to by destination.
  memcpy(buf,tmp[0],strlen(tmp[0]));//tsk
  
  // number of sites in original ref
  int totalSites = tmp.size();
  fprintf(stderr,"Ref has this many sites %i\n",totalSites);
  // a bit like the java Scanner, read line from file
  // cuts a string into tokens, where tokens seperated by delimeter
  strtok(buf,delims);
  int ncols=1;
  // reading first line in order to see nCol
  while(strtok(NULL,delims)){
    ncols++;
  }
  if(ncols<7){
    // has to have at least 7 columns 
    fprintf(stderr,"Too few cols, ncols=%d\n",ncols);
    exit(0);
  }

  // for which columns to keep
  ref.populations = new char*[nPop];
  
  // keep track of which original column it is
  int orgCol = 0;
  // keeps track of which new column (index = newCol-1) 
  // has to be above 0 for lookup in map
  int newCol = 1;
  char* columnID = strtok(tmp[0],delims);
  
  while(columnID!=NULL){
    orgCol++;
    // first 6 columns not freqs and has to at most K new columns included in ref
    if(orgCol>6 and newCol<=nPop){

      // if in supplied populations or if no populations supplied include
      if(includedPops.count(columnID)>0 or includedPops.empty()){
	// prints out which populations chosen
	fprintf(stderr,"Chosen pop %s\n",columnID);

	ref.populations[newCol-1] = strdup(columnID);
 	ref.popsToKeep[columnID] = newCol;
	// so that can translate from org (0 first freq column) to new column
	ref.colsToKeep[orgCol-7] = newCol;
	newCol++;
      }
    }

    columnID = strtok(NULL,delims);

   }
  
  // if K smaller than pops in ref K first chosen
  if((orgCol-6>nPop) and includedPops.empty()){
    fprintf(stderr,"You have chosen a K smaller than refPanel only %i first ref columns included\n",nPop);
  }

  // has to be number of sites in beagle
  ref.id = new char*[b.nSites];
  ref.chr = new int[b.nSites];
  ref.pos = new int[b.nSites];
  ref.name = new char*[b.nSites];
  ref.A0 = new char[b.nSites];
  ref.A1 = new char[b.nSites];
  ref.freqs = new double*[b.nSites];
  ref.refSites = b.nSites;


  // for keeping track of which index in new ref with
  // same number sites as in Beagle input
  int refIndex = 0;
  //then loop over the vector and parsing every line
  fprintf(stderr,"START\n");
  //fprintf(stderr,"Is %s Here %i\n",b.ids[8],ids[b.ids[8]]);
  for(int s=0;SIG_COND&& (s<totalSites);s++){

    // looking at id value chr_pos for detecting overlap
    char* id = strtok(tmp[s],delims);
    
    // check if site is in beagle file
    // otherwise continues to next site in ref
    
    if(overlap.count(id) > 0){
      ref.id[refIndex] = id;    
    }	else {
      continue;
    }
    
    ref.chr[refIndex] = atoi(strtok(NULL,delims));
    ref.pos[refIndex] = atoi(strtok(NULL,delims));
    ref.name[refIndex] = strdup(strtok(NULL,delims));
    // ref has A,C,G,T alleles
    ref.A0[refIndex] = strtok(NULL,delims)[0];
    ref.A1[refIndex] = strtok(NULL,delims)[0];
    ref.freqs[refIndex] = new double[ref.colsToKeep.size()];    

    //    reading in ref freqs
    for(int i=0;i<(ncols-6);i++){
    
      // check if org column to keep and thereby pop to keep in ref
      if(ref.colsToKeep.count(i)>0){
	// if bgl rs1 A B GL(AA) GL(AB) GL(BB) Then ref 1 1 rs1 B A 1-f(B)
	// minor is last allele in beagle file
	if(ref.A0[refIndex]==b.minor[refIndex]){
	  // new col has to be - 1 for right index, at least for looking up properly
	  ref.freqs[refIndex][ref.colsToKeep[i]-1] = 1 - atof(strtok(NULL,delims));

	} else{
	  ref.freqs[refIndex][ref.colsToKeep[i]-1] = atof(strtok(NULL,delims));
	}
	
	if(ref.freqs[refIndex][ref.colsToKeep[i]-1]<0){
	  fprintf(stderr,"Frequencies must be positive\n");
	  fprintf(stderr,"site %d, pop %d, has value %f\n",s,i,ref.freqs[refIndex][ref.colsToKeep[i]-1]);
	  exit(0);
	}

      } else{
	// if not column to include move to next column
	strtok(NULL,delims);
      }     
    }

    // count up index of new freq with N input beagle sites
    // only here if site was included in new ref
    refIndex++;


            
  }
  fprintf(stderr,"SLOW?\n");

  

   if(refIndex!=b.nSites){
     fprintf(stderr,"refIndex %i, b.nSites %i\n",refIndex,b.nSites);
     fprintf(stderr,"Sites in beagle file there are not in reference panel, should be subset of ref panel!\n");
      exit(0);
    }
  // Here additional stuff is calculated from the likelihoods
  // this must be done again while filtering later in main


  /* ret.keeps=new char*[ret.nSites]; // array nSites x nInd 0 if missing info
  ret.keepInd = new int[ret.nSites];
  ret.mafs = new float[ret.nSites];
  for(int s=0;s<ret.nSites;s++){
    ret.keeps[s] = new char[ret.nInd];
    int nKeep =0;
    for(int i=0;i<ret.nInd;i++){
      double mmin=std::min(ret.genos[s][i*3],std::min(ret.genos[s][i*3+1],ret.genos[s][i*3+2]));
      double mmax=std::max(ret.genos[s][i*3],std::max(ret.genos[s][i*3+1],ret.genos[s][i*3+2]));
      if(fabs(mmax-mmin)<misTol)
	ret.keeps[s][i] = 0;
      else{
	ret.keeps[s][i] =1;
	nKeep++;
      }
    }
    ret.keepInd[s] = nKeep;
    ret.mafs[s] = emFrequency(ret.genos[s],ret.nInd,MAF_ITER,MAF_START,ret.keeps[s],ret.keepInd[s]);
    } */
  fprintf(stderr,"This many in ref %i\n",refIndex);

  /* FOR CLEARING UP tmp VECTOR
  for(std::vector<char*>::iterator it = tmp.begin(); it != tmp.end(); ++it){
     free(*it);
  } 
  tmp.clear(); */

  
 
  gzclose(fp); //clean up filepointer
  return ref;
}



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


  std::vector<char*> tmp;
  // while still some left of string from file
  while(gzgets(fp,buf,LENS)){
    // duplicates string and puts it into tmp vector
    tmp.push_back(strdup(buf));
  }

  // for reading in a bigger ref panel
  // and getting intersecting sites
  int freqSites = tmp.size();
  fprintf(stderr,"This many ref sites: %i\n",freqSites);

  
  for(int i=0;i<freqSites;i++){
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

// for reading number of individuals in each ref
void readDouble1d(double *d,int nPop,const char*fname, std::map<char*,int,cmp_char> &popsToKeep){
  fprintf(stderr,"opening : %s with nPop=%d\n",fname,nPop);
  const char*delims=" \n\t";
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
   
  std::vector<char*> tmp;
  std::map<int,int> toKeep;
  // while still some left of string from file
  while(fgets(buf,lens,fp)){
    // duplicates string and puts it into tmp vector

    tmp.push_back(strdup(buf));
    // puts id of all sites in map for fast lookup
  }
  // going to all values ind nInd file
  char * word = strtok(tmp[0],delims);
  // keeps track of which value at in nInd file
    
  int orgCol = 0;
  int newCol = 1;

  
 
  
  while(word!=NULL){

    if(popsToKeep.count(word)>0){
      // because map index has to start at 1 for count method to work
      toKeep[orgCol] = newCol;
      newCol++;
    }

    word = strtok(NULL,delims);
    orgCol++;
  }
  free(tmp[0]);
  
  word = strtok(tmp[1],delims);
  int index = 0;
  while(word!=NULL){
    if(toKeep.count(index)>0){
      // because map index has to start at 1 for count method to work
      d[toKeep[index]-1] = atof(word);
    }

    word = strtok(NULL,delims);
    index++;
  }
  free(tmp[1]);

  
  fclose(fp);
}

void printDouble(double **ret,size_t x,size_t y, int highestLike, char** populations ,FILE *fp){
  for(size_t i=0;i<x;i++){
    if(i==0){
      for(size_t j=0;j<y;j++){
	fprintf(fp,"%s ",populations[j]);
      }
      fprintf(fp,"\n");
    }
    if(i<10 and i == highestLike){
      for(size_t j=0;j<y;j++){
	fprintf(fp,"%.4f ",ret[i][j]);
      }
      fprintf(fp,"\n");
    } else if(i>=10){
      for(size_t j=0;j<y;j++){
	fprintf(fp,"%.4f ",ret[i][j]);
      }
      fprintf(fp,"\n");
    }
  }
  
}

void printDoubleGz(double **ret,size_t x,size_t y, char** populations ,gzFile fp){

 for(size_t i=0;i<x;i++){
    if(i==0){
      for(size_t j=0;j<y;j++){
	gzprintf(fp,"%s ",populations[j]);
      }
      gzprintf(fp,"\n");
    }
    for(size_t j=0;j<y;j++){
      gzprintf(fp,"%.4f ",ret[i][j]);
    }
    gzprintf(fp,"\n");
  }
}


double likelihoodEmil(double* Q, double** F,int nSites, int nPop,double **genos){
  //  map2domainFEmil(F, nSites, nPop);
  //  map2domainQEmil(Q,nPop);

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
	d.ids[j] = strdup(dOrg.ids[row]);
      }
    }
  }
}

void emUnadjusted(double* Q, double** F, int nSites, int nPop,double **genos,double *Q_1) {
  double sumAG[nPop];
  double sumBG[nPop];
  // makes sure neither F nor Q has 0 values
  map2domainFEmil(F, nSites, nPop);
  map2domainQEmil(Q,nPop);
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

  // makes sure neither F nor Q has 0 values
  map2domainFEmil(F, nSites, nPop);
  map2domainFEmil(F_org, nSites, nPop);
  map2domainQEmil(Q,nPop);
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
  map2domainFEmil(F,nSites,nPop);
  map2domainFEmil(F_1,nSites,nPop);
}

int emAccelUnadjusted(const bgl &d, int nPop, double **F,double *Q,double *Q_new,double &lold, accFQ &a){
  //maybe these should be usersettable?
  double stepMin = 1;
  double stepMax0 = 1;
  static double stepMax = stepMax0;
  double mstep = 2;
  double objfnInc = 1;
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
  // changed so that alpha does not go above 10
  // got inaccurate results sometimes if to big alpha jumps
  alpha = std::min(alpha,maxAlpha);  
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


int emAccelUnadjustedV2(const bgl &d, double* nInd, int nPop, double** F,double* Q,double* &Q_new,double &lold, double** F_org, int nit, int boot){
 
//maybe these should be usersettable?
  double stepMin = 1;
  double stepMax0 = 1;
  static double stepMax = stepMax0;
  double mstep = 4;
  double objfnInc = 1;

  //we make these huge structures static such that we just allocate them the first time
  // static means they are NOT dynamically deallocated when function is over
  // way to declare variables if you want them to exist after function

  static double *Q_em1 = NULL;
  static double *Q_diff1 = NULL;
  static double *Q_em2 = NULL;
  static double *Q_diff2 = NULL;
  static double *Q_diff3 = NULL;
  static double *Q_tmp = NULL;
  static double *Q_tmpDiff = NULL;
  
  if(Q_em1==NULL){
    Q_em1 = new double[nPop];
    Q_diff1 = new double[nPop];
    Q_em2 = new double[nPop];
    Q_diff2 = new double[nPop];
    Q_diff3 = new double[nPop];
    Q_tmp = new double[nPop];
    Q_tmpDiff = new double[nPop];


  }
  //should cleanup and exit
  if(Q==NULL){
    delete [] Q_em1;
    delete [] Q_em2;
    delete [] Q_diff1;
    delete [] Q_diff2;
    delete [] Q_diff3;
    delete [] Q_tmp;
    delete [] Q_tmpDiff;
    return 0;
  }


  
  // first EM run
  emUnadjusted(Q, F, d.nSites, nPop, d.genos, Q_em1);
  // diff between 0 and 1 assigned to diff1
 
  minus1dEmil(Q_em1,Q,nPop,Q_diff1);
  // calculates "norm" of F_diff1 and Q_diff1 and sums them
  //  double sr2 = sumSquareEmil(F_diff1,d.nSites,nPop)+sumSquare1dEmil(Q_diff1,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster

  // Have F in here somehow
  double sr2 = sumSquare1dEmil(Q_diff1,nPop);
  // checks if convergence
  // tol is 1e-05 - line 62
  if(sqrt(sr2)<tol){
    fprintf(stderr,"like is %f\n",likelihoodEmil(Q_new, F, d.nSites, nPop,d.genos));    
    return 0;  
  }
  // second EM run
  emUnadjusted(Q_em1, F, d.nSites, nPop, d.genos, Q_em2);
  // diff between 1 and 2 assigned to diff2  
  minus1dEmil(Q_em2,Q_em1,nPop,Q_diff2);
  // calculates "norm" of F_diff2 and Q_diff2 and sums them 
  //double sq2 = sumSquareEmil(F_diff2,d.nSites,nPop)+sumSquare1dEmil(Q_diff2,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
  double sq2 = sumSquare1dEmil(Q_diff2,nPop);
  // checks if convergence - a second time
  // tol is 1e-05 - line 62
  if(sqrt(sq2)<tol){
    fprintf(stderr,"like is %f\n",likelihoodEmil(Q_new, F, d.nSites, nPop,d.genos));
    return 0;
    //    break;
  }
  // diff between diff1 and diff2
  minus1dEmil(Q_diff2,Q_diff1,nPop,Q_diff3);
  // WHY NOT SQ2 - because it is v = diff2 - diff1 actually
  //double sv2 = sumSquareEmil(F_diff3,d.nSites,nPop)+sumSquare1dEmil(Q_diff3,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
  double sv2 = sumSquare1dEmil(Q_diff3,nPop);
  // does not have to make alpha negative because changed accordingly in equations
  double alpha = sqrt(sr2/sv2);
  
  // makes sure alpha does not go below 1 and above stepMax
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  // changed so that alpha does not go above 10
  // got inaccurate results sometimes if to big alpha jumps
 
  //alpha = std::min(alpha,maxAlpha); 
  // updates F with new values via alpha - where alpha controls rate of change
  

  // check that F has no entries of 0
  for(size_t i=0;i<nPop;i++){
    // because alpha pos in code, does not need minus -(-alpha) = alpha
    Q_new[i] = Q[i]+2*alpha*Q_diff1[i]+alpha*alpha*Q_diff3[i];
    Q_tmp[i] = 1.0;
  }
  
  //double** tmpAdrF = F_tmp;
  //double* tmpAdrQ = Q_tmp;
  map2domainQEmil(Q_new,nPop);
  // if alpha not too close (0.01 close) to 1 
  if (fabs(alpha - 1) > 0.01){
  
    // we estimate new Q and F, with our inferred Q and F via alpha
    emUnadjusted(Q_new, F, d.nSites, nPop,d.genos,Q_tmp);
    // assign value of Q_tmp and F_tmp to Q_new and F_new
    // in order to avoid Q_tmp or F_tmp pointing to Q and F
    //Q_tmp=tmpAdrQ;
    //    F_tmp=tmpAdrF;
   
    // I might have to create another one instead of _diff3, where _tmp - _new
    
    minus1dEmil(Q_tmp,Q_new,nPop,Q_tmpDiff);
    // WHY NOT SQ2 - because it is v = diff2 - diff1 actually
    //double sv2 = sumSquareEmil(F_diff3,d.nSites,nPop)+sumSquare1dEmil(Q_diff3,nPop);
    // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
    double res = sumSquare1dEmil(Q_tmpDiff,nPop);
    double parnorm = (1/std::sqrt(nPop))*sumSquare1dEmil(Q_tmpDiff,nPop);
    double kres = 1 + parnorm + sq2;
    if(res <= kres){
      std::swap(Q_new,Q_tmp);
    } else{
      std::swap(Q_new,Q_em2);
    }

    if(res > kres){
      if (alpha == stepMax){
	stepMax = std::max(stepMax0, stepMax/2);
      }
      alpha = 1;
    }
      
  }

  if (alpha == stepMax){ 
    stepMax = mstep * stepMax;
  }
  if (stepMin < 0 & alpha == stepMin) {
    stepMin = mstep * stepMin;
  }
  
  if(nit % 10 == 0){
    if(boot == 0){
      double lnew = likelihoodEmil(Q_new, F, d.nSites, nPop,d.genos);
      fprintf(stderr,"iter[%d] like=%f sr2=%f sq2=%f alpha=%f ",nit,lnew,sr2,sq2,alpha);
      for(int i=0;i<nPop;i++){	      
	fprintf(stderr,"Q=%f, ",Q_new[i]);
      }
      fprintf(stderr,"\n");
    }
  }
  
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
  
  double lold2 = likelihoodEmil(Q, F, d.nSites, nPop,d.genos);
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
  fprintf(stderr,"alpha: %f\n",alpha);
  // makes sure alpha does not go below 1 and above stepMax
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  // changed so that alpha does not go above 10
  // got inaccurate results sometimes if to big alpha jumps
  

  //alpha = std::min(alpha,maxAlpha); 



  // updates F with new values via alpha - where alpha controls rate of change
  fprintf(stderr,"alpha: %f %f\n",alpha,stepMax);
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

  fprintf(stderr,"Pointers: %p %p %p %p %p %p\n",a.Q_tmp,a.F_tmp,Q_new,F_new,Q,F);
  
  double** tmpAdrF = a.F_tmp;
  double* tmpAdrQ = a.Q_tmp;
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
    fprintf(stderr,"Pointers: %p %p %p %p %p %p\n",a.Q_tmp,a.F_tmp,Q_new,F_new,Q,F);
    //if (alpha == stepMax){
      //  stepMax = std::max(stepMax0, stepMax/2);
      //}

  }
  double lnew = likelihoodEmil(Q_new, F_new, d.nSites, nPop,d.genos);
  fprintf(stderr,"likes %f %f\n", lold2,lnew);
  
  if (lnew < lold2) {
    fprintf(stderr,"bad guess in %s\n", __FUNCTION__);
    //    *Q_new = Q_em2;
    //    *F_new =F_em2;
    if (alpha == stepMax){
      stepMax = std::max(stepMax0, stepMax/mstep);
    }
    alpha = 1;
  }

  if (std::abs(lnew - lold)<0.0001) {
    fprintf(stderr,"bad guess in %s\n", __FUNCTION__);
    //    *Q_new = Q_em2;
    //    *F_new =F_em2;
    if (alpha == stepMax){
      stepMax = std::max(stepMax0, stepMax/mstep);
    }
    alpha = 1;
  } 
  fprintf(stderr,"alpha %f\n", alpha);
  if (alpha == stepMax) {
    stepMax = mstep*stepMax;
  }
  lold=lnew;
  return 1;
}


// based on squarem1, from SQUAREM R package, by RAVI VARADHAN and CHRISTOPHE ROLAND Scandinavian Journal of Statistics, Vol. 35: 335–353, 2008
int emAccelEmilV2(const bgl &d, double* nInd, int nPop, double** F,double* Q,double** &F_new,double* &Q_new,double &lold, double** F_org, int nit, int boot){
  //maybe these should be usersettable?
  double stepMin = 1;
  double stepMax0 = 1;
  static double stepMax = stepMax0;
  double mstep = 4;
  double objfnInc = 1;

  //we make these huge structures static such that we just allocate them the first time
  // static means they are NOT dynamically deallocated when function is over
  // way to declare variables if you want them to exist after function
  static double **F_em1 = NULL;
  static double *Q_em1 = NULL;
  static double **F_diff1 = NULL;
  static double *Q_diff1 = NULL;
  static double **F_em2 = NULL;
  static double *Q_em2 = NULL;
  static double **F_diff2 = NULL;
  static double *Q_diff2 = NULL;
  static double **F_diff3 = NULL;
  static double *Q_diff3 = NULL;
  static double **F_tmp = NULL;
  static double *Q_tmp = NULL;
  static double **F_tmpDiff = NULL;
  static double *Q_tmpDiff = NULL;
  
  if(F_em1==NULL){
    fprintf(stderr,"I AM HERE!!!!!!!!!!\n");
    F_em1 = allocDouble(d.nSites,nPop);
    Q_em1 = new double[nPop];
    F_diff1 = allocDouble(d.nSites,nPop);
    Q_diff1 = new double[nPop];
    F_em2 = allocDouble(d.nSites,nPop);
    Q_em2 = new double[nPop];
    F_diff2 = allocDouble(d.nSites,nPop);
    Q_diff2 = new double[nPop];
    F_diff3 = allocDouble(d.nSites,nPop);
    Q_diff3 = new double[nPop];
    F_tmp = allocDouble(d.nSites,nPop);
    Q_tmp = new double[nPop];
    F_tmpDiff =  allocDouble(d.nSites,nPop);
    Q_tmpDiff = new double[nPop];


  }
  //should cleanup and exit
  if(F==NULL){
    dalloc(F_em1,d.nSites);
    dalloc(F_em2,d.nSites);
    dalloc(F_diff1,d.nSites);
    dalloc(F_diff2,d.nSites);
    dalloc(F_diff3,d.nSites);
    dalloc(F_tmp,d.nSites);
    dalloc(F_tmpDiff,d.nSites);
    delete [] Q_em1;
    delete [] Q_em2;
    delete [] Q_diff1;
    delete [] Q_diff2;
    delete [] Q_diff3;
    delete [] Q_tmp;
    delete [] Q_tmpDiff;
    return 0;
}
    
  // first EM run
  emEmil(Q, F, d.nSites, nInd, nPop,d.genos, F_em1, Q_em1, F_org);
  // diff between 0 and 1 assigned to diff1
  minusEmil(F_em1,F,d.nSites,nPop,F_diff1);
  minus1dEmil(Q_em1,Q,nPop,Q_diff1);
  // calculates "norm" of F_diff1 and Q_diff1 and sums them
  //  double sr2 = sumSquareEmil(F_diff1,d.nSites,nPop)+sumSquare1dEmil(Q_diff1,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster

  // Have F in here somehow
  double sr2 = sumSquare1dEmil(Q_diff1,nPop) + sumSquareEmil(F_diff1,d.nSites,nPop);
  // checks if convergence
  // tol is 1e-05 - line 62
  if(sqrt(sr2)<tol){
    fprintf(stderr,"like is %f\n",likelihoodEmil(Q_new, F_new, d.nSites, nPop,d.genos));    
    return 0;
    
    //break;
  }
  // second EM run
  emEmil(Q_em1, F_em1, d.nSites, nInd, nPop,d.genos, F_em2, Q_em2, F_org);
  // diff between 1 and 2 assigned to diff2  
  minusEmil(F_em2,F_em1,d.nSites,nPop,F_diff2);
  minus1dEmil(Q_em2,Q_em1,nPop,Q_diff2);
  // calculates "norm" of F_diff2 and Q_diff2 and sums them 
  //double sq2 = sumSquareEmil(F_diff2,d.nSites,nPop)+sumSquare1dEmil(Q_diff2,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
  double sq2 = sumSquare1dEmil(Q_diff2,nPop) + sumSquareEmil(F_diff2,d.nSites,nPop);
  // checks if convergence - a second time
  // tol is 1e-05 - line 62
  if(sqrt(sq2)<tol){
    fprintf(stderr,"like is %f\n",likelihoodEmil(Q_new, F_new, d.nSites, nPop,d.genos));
    return 0;
    //    break;
  }
  // diff between diff1 and diff2
  minusEmil(F_diff2,F_diff1,d.nSites,nPop,F_diff3);
  minus1dEmil(Q_diff2,Q_diff1,nPop,Q_diff3);
  // WHY NOT SQ2 - because it is v = diff2 - diff1 actually
  //double sv2 = sumSquareEmil(F_diff3,d.nSites,nPop)+sumSquare1dEmil(Q_diff3,nPop);
  // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
  double sv2 = sumSquare1dEmil(Q_diff3,nPop) + sumSquareEmil(F_diff3,d.nSites,nPop);
  // does not have to make alpha negative because changed accordingly in equations
  double alpha = sqrt(sr2/sv2);
  
  // makes sure alpha does not go below 1 and above stepMax
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  // changed so that alpha does not go above 10
  // got inaccurate results sometimes if to big alpha jumps
 
  //alpha = std::min(alpha,maxAlpha); 
  // updates F with new values via alpha - where alpha controls rate of change
  
  for(size_t i=0;i<d.nSites;i++){
    for(size_t j=0;j<nPop;j++){
      // based on ngsAdmix approach
      F_new[i][j] = F[i][j]+2*alpha*F_diff1[i][j]+alpha*alpha*F_diff3[i][j];
      F_tmp[i][j] = 1.0;
    }
  }
  // check that F has no entries of 0
  map2domainFEmil(F_new,d.nSites,nPop);
  for(size_t i=0;i<nPop;i++){
    // because alpha pos in code, does not need minus -(-alpha) = alpha
    Q_new[i] = Q[i]+2*alpha*Q_diff1[i]+alpha*alpha*Q_diff3[i];
    Q_tmp[i] = 1.0;
  }
  
  //double** tmpAdrF = F_tmp;
  //double* tmpAdrQ = Q_tmp;
  map2domainQEmil(Q_new,nPop);
  // if alpha not too close (0.01 close) to 1 
  if (fabs(alpha - 1) > 0.01){
  
    // we estimate new Q and F, with our inferred Q and F via alpha
    emEmil(Q_new, F_new, d.nSites, nInd, nPop,d.genos,F_tmp,Q_tmp,F_org);
    // assign value of Q_tmp and F_tmp to Q_new and F_new
    // in order to avoid Q_tmp or F_tmp pointing to Q and F
    //Q_tmp=tmpAdrQ;
    //    F_tmp=tmpAdrF;
   
    // I might have to create another one instead of _diff3, where _tmp - _new
    minusEmil(F_tmp,F_new,d.nSites,nPop,F_tmpDiff);
    minus1dEmil(Q_tmp,Q_new,nPop,Q_tmpDiff);
    // WHY NOT SQ2 - because it is v = diff2 - diff1 actually
    //double sv2 = sumSquareEmil(F_diff3,d.nSites,nPop)+sumSquare1dEmil(Q_diff3,nPop);
    // SQUAREM also only bases theta on Q, I do not change F very much - this is faster
    double res = sumSquare1dEmil(Q_tmpDiff,nPop) + sumSquareEmil(F_tmpDiff,d.nSites,nPop);


    double parnorm = (1/std::sqrt(nPop))*sumSquare1dEmil(Q_tmpDiff,nPop) + (1/std::sqrt(d.nSites*nPop))*sumSquareEmil(F_tmpDiff,d.nSites,nPop);
    double kres = 1 + parnorm + sq2;


    if(res <= kres){
      std::swap(Q_new,Q_tmp);
      std::swap(F_new,F_tmp);
    } else{
      std::swap(Q_new,Q_em2);
      std::swap(F_new,F_em2);
    }

    if(res > kres){
      if (alpha == stepMax){
	stepMax = std::max(stepMax0, stepMax/2);
      }
      alpha = 1;
    }
      
  }

  if (alpha == stepMax){ 
    stepMax = mstep * stepMax;
  }
  if (stepMin < 0 & alpha == stepMin) {
    stepMin = mstep * stepMin;
  }
  

  
  if(nit % 10 == 0){
    if(boot == 0){
      double lnew = likelihoodEmil(Q_new, F_new, d.nSites, nPop,d.genos);
      fprintf(stderr,"iter[%d] like=%f sr2=%f sq2=%f alpha=%f ",nit,lnew,sr2,sq2,alpha);
      for(int i=0;i<nPop;i++){	      
	fprintf(stderr,"Q=%f, ",Q_new[i]);
      }
      fprintf(stderr,"\n");
    }
  }
  
  return 1;
}


std::map <char*,int,cmp_char> findOverlap(const char* lname, const char* plinkName, const char* fname){
  std::map <char*,int,cmp_char> inputSites;
  static std::map <char*,int,cmp_char> overlap;
  const char *delims = "\t \n";
  std::vector<char*> tmp3;
  std::vector<char*> tmp4;
  if(plinkName==NULL){

    gzFile fp1 = NULL;
    // checking if file can be opened
    if(Z_NULL==(fp1=gzopen(lname,"r"))){
      fprintf(stderr,"Error opening file: %s\n",lname);
      exit(0);
    }
    // creaing string of length LENS
    char buf1[LENS];
    
    // while still some left of string from file
    while(gzgets(fp1,buf1,LENS)){
      // duplicates string and puts it into tmp vector
      tmp3.push_back(strdup(buf1));
    }
    //then loop over the vector and parsing every line
    for(int s=1;SIG_COND&& (s<tmp3.size());s++){
      inputSites[strtok(tmp3[s],delims)]=1;    
    }
    
    gzclose(fp1);
  } else{

    plink* pl_tmp = readplink(plinkName);
    
    for(int s=0;s<(pl_tmp->y);s++){
      if(pl_tmp->d[0][s]==3){
	continue;
      }
      inputSites[strdup(pl_tmp->bim.id[s])]=1;
      
    }
    
    kill_plink(pl_tmp);


    
  }
  gzFile fp2 = NULL;
  char buf2[LENS];
  
  // checking if file can be opened
  if(Z_NULL==(fp2=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }

  //find number of columns
  // reading string from compressed file 
  
   
  // while still some left of string from file
  while(gzgets(fp2,buf2,LENS)){
    // duplicates string and puts it into tmp vector
    tmp4.push_back(strdup(buf2));    
  }

 
  // starts at 1 to avoid header
  for(int s=1;SIG_COND&& (s<tmp4.size());s++){

    // looking at id value chr_pos for detecting overlap
    char* id = strtok(tmp4[s],delims);
    
    // check if site is in beagle file
    // otherwise continues to next site in ref
    
    if(inputSites.count(id) > 0){
      overlap[strdup(id)] = 1;    
    }	else {
      continue;
    }
  }
  gzclose(fp2);
  if(lname!=NULL){
    for(int s=0;SIG_COND&& (s<tmp4.size());s++){
      free(tmp4[s]);
    }
    for(int s=0;SIG_COND&& (s<tmp3.size());s++){
      free(tmp3[s]);
    }
  } else if(plinkName!=NULL){
    for(int s=0;SIG_COND&& (s<tmp4.size());s++){
      free(tmp4[s]);
    }
  }

  
  return(overlap);

}

void info(){
  
  fprintf(stderr,"Arguments:\n");
  fprintf(stderr,"\t-likes Beagle likelihood filename\n");
  fprintf(stderr,"\t-plinks Plink file in the binary bed format\n");
  fprintf(stderr,"\t-K Number of ancestral populations - default value is 3\n"); 
  fprintf(stderr,"\t-Nname Number of individuals in each reference populations\n");
  fprintf(stderr,"\t-fname Ancestral population frequencies\n");
  fprintf(stderr,"Optional:\n");
 
  fprintf(stderr,"\t-outfiles Prefix for output files\n"); 
  fprintf(stderr,"\t-printInfo print ID and mean maf for the SNPs that were analysed\n"); 
  fprintf(stderr,"\t-printFreq print admixture adjusted allele frequencies of reference panel + input individual (1: yes, 0: no (default))\n"); 

  fprintf(stderr,"Setup:\n");
  fprintf(stderr,"\t-whichPops Which populations from the reference panel to include in analysis, must be comma seperated (pop1,pop2,..)\n");
  fprintf(stderr,"\t-doAdjust Adjusts the frequencies in the reference populations with the input (1: yes (default), 0: no)\n");
  fprintf(stderr,"\t-seed Seed for initial guess in EM\n"); 
  fprintf(stderr,"\t-P Number of threads\n"); 
  fprintf(stderr,"\t-method If 0 no acceleration of EM algorithm (1: yes (default), 0: no)\n"); 
  fprintf(stderr,"\t-misTol Tolerance for considering site as missing\n");

  fprintf(stderr,"Stop chriteria:\n"); 
  fprintf(stderr,"\t-Qconv Stopping criteria based on change in Q (works best when using doAdjust) (1: yes, 0: no (default))\n"); 
  fprintf(stderr,"\t-Qtol Tolerance value for stopping criteria based on change in Q (0.001 (default))\n"); 
  fprintf(stderr,"\t-tolLike10 Loglikelihood difference in 10 iterations\n"); 
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

double calcThresEmil(double *d1,double *d2, int x){
  // finds the largest difference between 2 arrays
  // arrays has dimention x times y
  double diff=fabs(d1[0]-d2[0]);
  for(int i=1;i<x;i++){
    if(fabs(d1[i]-d2[i])<diff){
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

  int maxIter = 2000;
  int method = 1;
  int printInfo = 0;
  int printFreq = 0;
  int doAdjust = 1;
  float minMaf = 0.00;
  const char* lname = NULL;
  // is this an issue that it is not const??

  const char* fname = NULL;
  const char* pops = NULL;
  const char* Nname = NULL;
  const char* outfiles = NULL;
  const char* plinkName = NULL;
  int nPop = 3;
  int seed = time(NULL);
  // float tolLike50=0.1;
  float tolLike10=0.00001; // changed by Emil
  int nBoot = 0; // for bootstrapping
  int Qconv = 0;
  double Qtol = 0.001;
  // reading arguments
  argv++;
  while(*argv){
    // GL in the shape of beagle file
    if(strcmp(*argv,"-likes")==0 || strcmp(*argv,"-l")==0) lname=*++argv; //name / char arrays
    // probably will need this, implicitly in the freqs file, by default 3
    else if(strcmp(*argv,"-plink")==0 || strcmp(*argv,"-p")==0) plinkName=*++argv;
    else if(strcmp(*argv,"-K")==0) nPop=atoi(*++argv); 
    // would also need number of individuals in each ref category
    // to read start values from output from previous run 
    else if(strcmp(*argv,"-fname")==0 || strcmp(*argv,"-f")==0) fname=*++argv; 
    // for reading in number of individauls in each population
    else if(strcmp(*argv,"-Nname")==0 || strcmp(*argv,"-N")==0) Nname=*++argv;
    // prefix for output files
    else if(strcmp(*argv,"-outfiles")==0 || strcmp(*argv,"-o")==0) outfiles=*++argv; 
    else if(strcmp(*argv,"-whichPops")==0 || strcmp(*argv,"-pops")==0) pops=*++argv; 
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
    else if(strcmp(*argv,"-tolLike10")==0||strcmp(*argv,"-lt10")==0) tolLike10=atof(*++argv); //float/double - atof - char array to double/float
    // do I need those when I have only one individual??
    else if(strcmp(*argv,"-bootstrap")==0||strcmp(*argv,"-boot")==0) nBoot=atoi(*++argv);
    else if(strcmp(*argv,"-tol")==0||strcmp(*argv,"-t")==0) tol=atof(*++argv);
    else if(strcmp(*argv,"-maxiter")==0 || strcmp(*argv,"-i")==0) maxIter=atoi(*++argv); 
    // different filterings 
    else if(strcmp(*argv,"-misTol")==0 || strcmp(*argv,"-mt")==0) misTol=atof(*++argv);
    // might not really work, not sure if should apply to refPanel?
    else if(strcmp(*argv,"-minMaf")==0||strcmp(*argv,"-maf")==0) minMaf=atof(*++argv);

    // different genotype callers - tolerance which we do use
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      info();
      return 0;
    }
    ++argv;
  }


  //check that non optional options have been used. 
  if(lname==NULL and plinkName==NULL){
    fprintf(stderr,"Please supply a beagle or plink input file: -likes or -plink");
    info();
  } else if(fname==NULL){
    fprintf(stderr,"Please supply a reference panel: -fname");
    info();
  } else if(lname!=NULL and plinkName!=NULL){
    fprintf(stderr,"Please supply ONLY a beagle or plink input file, not BOTH: -likes or -plink");
    info();
  }
  if(outfiles==NULL){
    fprintf(stderr,"Will use beagle name as prefix for output\n");
    outfiles=lname;
  }

  //out put files
  FILE *flog=openFile(outfiles,".log");
  FILE *ffilter=openFile(outfiles,".filter");

  fprintf(stderr,"Input: lname=%s nPop=%d, fname=%s outfiles=%s\n",lname,nPop,fname,outfiles);
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
  fprintf(stderr,"Convergence: maxIter=%d tol=%f tolLike10=%f\n",maxIter,tol,tolLike10);
  fprintf(stderr,"The following number of bootstraps have been chosen: %i\n",nBoot);
  if(Qconv>0){
    fprintf(stderr,"Convergence via difference in Q values chosen, threshold of: %f\n",Qtol);    
  }
  fprintf(stderr,"Filters: misTol=%f minMaf=%f\n",misTol,minMaf);

  fprintf(flog,"Input: lname=%s nPop=%d, fname=%s outfiles=%s\n",lname,nPop,fname,outfiles);
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
  fprintf(flog,"Convergence: maxIter=%d tol=%f tolLike10=%f\n",maxIter,tol,tolLike10);
  fprintf(flog,"The following number of bootstraps have been chosen: %i\n",nBoot);
  if(Qconv>0){
    fprintf(flog,"Convergence via difference in Q values chosen, threshold of: %f\n",Qtol);    
  }
  fprintf(flog,"Convergence: maxIter=%d tol=%f tolLike10=%f\n",maxIter,tol,tolLike10);
  fprintf(flog,"Filters: misTol=%f minMaf=%f\n",misTol,minMaf);


  errTolStart = errTolMin;
  errTol = errTolMin;
    
  clock_t t = clock();//how long time does the run take
  time_t t2 = time(NULL);

 
  
  //read BEAGLE likelihood file  
  // made into object to give it additional info
  bgl d;
  bgl dOrg;
  std::map <char*,int,cmp_char> overlap;
  if(lname!=NULL){

    overlap = findOverlap(lname, NULL, fname);
    
    d=readBeagle(lname,overlap);
    dOrg=readBeagle(lname,overlap);
  } else if(plinkName!=NULL){
    overlap = findOverlap(NULL, plinkName, fname);
    d=readPlinkToBeagle(plinkName, overlap);
    dOrg=readPlinkToBeagle(plinkName, overlap);
    fprintf(stderr,"plink id %s\n",d.ids[0]);
    fprintf(stderr,"plink prob %f %f %f\n",d.genos[0][0],d.genos[0][1],d.genos[0][2]);
    
  } else{
    fprintf(stderr,"Must specify plink file\n");
    info();
  }
 fprintf(stderr,"Overlap: of %i sites between input and ref\n",overlap.size());
  
 
  
  //  fprintf(stderr,"plink last prob %f %f %f\n",d.genos[157946][0],d.genos[157946][1],d.genos[157946][2]);
   
  fprintf(stderr,"Input file has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);
  fprintf(flog,"Input file has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);

  // to get only some populations in refPanel
  std::map <char*,int,cmp_char> includedPops;
  
  // if pops are specified, reads which pops and construcs map with those
  if(pops!=NULL){

   
    char* pops1;
    strcpy(pops1,pops);
    char* temp = strtok(pops1,",");
    while(temp!=NULL){

      includedPops[temp] = 1;
      temp = strtok(NULL,",");
    }

    // checks pops and K do not differ in size
    if(includedPops.size()!=nPop){
      fprintf(stderr,"K and population selected differ: K=%i whichPops=%lu, must be equal!\n",nPop,includedPops.size());
      fprintf(flog,"K and population selected differ: K=%i whichPops=%lu, must be equal!\n",nPop,includedPops.size());
      info();

    }
  }
  // to find out which version of C++
  //fprintf(stderr,"V: [%ld] ", __cplusplus);
  
  // filter sites based on MAF - 

  if(printInfo){
    printKeepSites(d,ffilter);
  }
  fflush(stderr);

  // seed for bootstrapping, have only one seed!
  std::srand(seed);

  // reads in ref Panel
  refPanel ref;
  ref = readRefPanel(fname,d,includedPops,nPop,overlap);
  fprintf(stderr,"Done!\n");
    
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

  // gets freq values from ref struct, where read into
  for(int i=0;i<d.nSites;i++){
    for(int k=0;k<nPop;k++){
      F[i][k] = ref.freqs[i][k];
      F_org[i][k] = ref.freqs[i][k];
      F_orgOrg[i][k] = ref.freqs[i][k];
    }
     
  }

  // because has to have 10 values for normal data, in order to have converge
  // and then nBoot bootstrapped values
  double **Q = allocDouble(nBoot+10,nPop);
  double **Q_new = allocDouble(nBoot+10,nPop);
   
  double *N = new double[nPop];

  // putting in an intial guess 
  double *sum = new double[nBoot+10];
  // nBoot + 1 rows
  for(int j=0;j<nBoot+10;j++){
    sum[j]=0;
    for(int k=0;k<nPop;k++){
      Q[j][k]=rand()*1.0/RAND_MAX;
      sum[j]+=Q[j][k];
    }
  }
  // nBoot + 1 rows
  for(int j=0;j<nBoot+10;j++){
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


  // reading nInd, where colsToKeep from ref specified
  // to read in the same columns as in ref
  readDouble1d(N,nPop,Nname,ref.popsToKeep);

  for(int i=0;i<nPop;i++){
    fprintf(stderr,"N = %f\n",N[i]);
    fprintf(flog,"Chosen pop %s\n",ref.populations[i]);
  }

  double* bestLike = new double[10];
  int highestLike = 0;
   
  //update the global stuff NOW
  //update the internal stuff in the pars for the threading
  double lold = likelihoodEmil(Q[0], F_org, d.nSites, nPop,d.genos);
  
  fprintf(stderr,"iter[start] like is=%f\n",lold);

  //////////////////////////////////////// em ///////////////////////////////////  
  //below is the main looping through the iterations.
  // we have 4 possible ways, unAdjusted/Adjusted basic EM/Accelerated EM
  int nit = 0;
  double likeLast = lold;
  double lastQthres = 0;

  // first 10 runs for converge with new Q starting point
  // then nBoot new runs for bootstrapping with best Q starting point
  for(int b=0;SIG_COND and b<(nBoot+10);b++) {
    nit=0;

      for(nit=1;SIG_COND and nit<maxIter;nit++) {	
	if(doAdjust==0){
	  if(method==0){
	    emUnadjusted(Q[b], F, d.nSites, nPop,d.genos,Q_new[b]);
	  } else{
	  

	    if(0==emAccelUnadjustedV2(d, N, nPop, F, Q[b], Q_new[b],lold,F_org, nit, b)){
	      bestLike[b] = likelihoodEmil(Q[b], F_org, d.nSites, nPop,d.genos);
	      fprintf(stderr,"like after %f\n",bestLike[b]);

	      break;
	    }
	  }
	} else{   
	  if(method==0){
	    emEmil(Q[b], F, d.nSites, N, nPop,d.genos,F_new,Q_new[b],F_org);
	  } else{


	    //lold = likelihoodEmil(Q[b], F_org, d.nSites, nPop,d.genos);

	    if(0==emAccelEmilV2(d, N, nPop, F, Q[b], F_new, Q_new[b],lold,F_org, nit, b)){
	      bestLike[b] = likelihoodEmil(Q[b], F_new, d.nSites, nPop,d.genos);
	      fprintf(stderr,"like after %f\n", bestLike[b]);

	      break;
	    }

	  }
	  std::swap(F,F_new);
	}
      	std::swap(Q,Q_new);

	if((nit%10)==0 and method == 0){ //stopping criteria
      	  // thres is largest differense in admixture fractions
	  //fprintf(stderr,"delta Qs  %f last: %f\n",calcThresEmil(Q[b],Q_new[b],nPop),lastQthres);
	  //	  double lik = likelihoodEmil(Q[b], F_org, d.nSites, nPop,d.genos);
	  double lik = likelihoodEmil(Q[b], F, d.nSites, nPop,d.genos);
	  if(b==0){
	    fprintf(stderr,"iter[%d] last like is=%f thres=%f\t",nit,likeLast,calcThresEmil(Q[b],Q_new[b],nPop));
	    fprintf(stderr,"iter[%d] like is=%f thres=%f\t",nit,lik,calcThresEmil(Q[b],Q_new[b],nPop));
	    fprintf(stderr,"iter[%d] diff in likelihood is=%f\t",nit,std::abs(lik-likeLast));
	    fprintf(stderr,"iter[%d] ",nit);
	    for(int i=0;i<nPop;i++){	      
	      fprintf(stderr,"Q is=%f, ",Q[b][i]);

	    }
	    fprintf(stderr,"\n");      
	  }
	  if(calcThresEmil(Q[b],Q_new[b],nPop) < Qtol and Qconv>0) {
	    if(b==0){
	      fprintf(stderr,"Convergence achived because diffence in Q values less than %f\n",Qtol);
	    }
	    if(b<10){
	      bestLike[b] = lik;
	    }
  	    break;
	  } else /* if(fabs(calcThresEmil(Q[b],Q_new[b],nPop)-lastQthres)<1e-07){
	    fprintf(stderr,"delta Q %.20f\n",fabs(calcThresEmil(Q[b],Q_new[b],nPop)-lastQthres));
	    if(b==0){
	      fprintf(stderr,"Convergence achived because diffence in Q values less than 1e-07\n");
	  }
	    if(b<10){
	      bestLike[b] = lik;
	    }
	    break;
	    }
	    else */ if(std::abs(lik-likeLast) < tolLike10 and Qconv==0) {
	    if(b==0){
	      fprintf(stderr,"Convergence achived becuase log likelihooditer difference is less than %f\n",tol);
	    }
	    if(lik-likeLast<0){
	      if(b==0){
		fprintf(stderr,"Convergence achived because log likelihood difference was NEGATIVE\n");
	      }
	    }
	    if(b<10){
	      bestLike[b] = lik;
	    }
	    break;
	  }
	  
	  
	  likeLast = lik;

	}
	// only calc when run 9, 19 , ... (as just before conv check)
	/* if((nit%10)==9){
	  likeLast = likelihoodEmil(Q[b], F, d.nSites, nPop,d.genos); 
	  lastQthres = calcThresEmil(Q[b],Q_new[b],nPop);
	  } */
	
	
	if((nit+1)==maxIter){
	  bestLike[b] = likeLast;
	}
      }

     
      fprintf(stderr,"BOOM!\n");
      if(b<10){
	fprintf(stderr,"This many iterations %i for run %i\n",nit,b);
	fprintf(flog,"This many iterations %i for run %i\n",nit,b);
	
      }
      for(int i=0;i<d.nSites;i++){
	for(int k=0;k<nPop;k++){
	  // original ref panel freqs
	  F[i][k] = F_orgOrg[i][k];
	  F_new[i][k] = F_orgOrg[i][k];
	}
      }    
      // first run
      // nBoot + 1 rows
      // to find the Q with lowest likelihood, after 10 first runs
  if(b==9){
    double highLike = 0;
    for(int l=0;l<10;l++){
      if(l==0){
	highLike = bestLike[0];
      }else{
	if(highLike < bestLike[l]){
	  highLike = bestLike[l];
	  // highestLike is index of value with highest like
	  highestLike = l;
	}
      }
    }
    for(int j=10;j<(nBoot+10);j++){
      for(int k=0;k<nPop;k++){
	// make first estimated Q starting guess for rest
	Q[j][k] = Q[highestLike][k];
	Q_new[j][k] = Q[highestLike][k];
      }
    }
    // same for GL input!!
    std::swap(F,F_1stRun);
  }
  // F_orgOrg original F sampled from

  if(b>9){
    bootstrap(dOrg,d,F_orgOrg,F_org,F,nPop); 
    fprintf(stderr,"At this bootstrapping: %i out of: %i\n",b-9,nBoot);
  }
  }

  for(int i=0;i<10;i++){
    fprintf(stderr,"best like %f after %i!\n",bestLike[i],i);
    fprintf(flog,"best like %f after %i!\n",bestLike[i],i);
    for(int j=0;j<nPop;j++){
      fprintf(stderr,"Q %f ",Q[i][j]);
      fprintf(flog,"Q %f ",Q[i][j]);
    }
    fprintf(stderr," after %i!\n",i);
    fprintf(flog," after %i!\n",i);
  } 

  fprintf(flog,"estimated  Q = ");
  for(int i=0;i<nPop;i++){
    fprintf(flog,"%f ",Q[highestLike][i]);
  }
  fprintf(flog,"\n");

  fprintf(stderr,"best like %f after %i!\n",bestLike[highestLike],highestLike);
  fprintf(flog,"best like %f after %i!\n",bestLike[highestLike],highestLike);

  
  /////////////////////////////////////////////////////////// done - make output and clean /////////////////////////////////  
  // Print F and Q in files
  
  // change to mean CIlower CIupper
  FILE *fp=openFile(outfiles,".qopt");
  // nBoot + 1, because first value is estimated Q
  printDouble(Q,nBoot+10,nPop,highestLike,ref.populations,fp);
  fclose(fp);
  
  // only if certain flag
  if(printFreq>0){
    gzFile fpGz=openFileGz(outfiles,".fopt.gz");
    printDoubleGz(F_1stRun,d.nSites,nPop,ref.populations,fpGz);
    gzclose(fpGz);
  }
  double nop;
  if(doAdjust>0 and method>0){
    emAccelEmilV2(d, NULL, 0, NULL, NULL, F_new, Q_new[0],nop,NULL, 0,1);
  } else if(method > 0){
    emAccelUnadjustedV2(d, NULL, 0, NULL, NULL, Q_new[0],nop,NULL, 0, 1);   
  }

  //deallocate memory 
  dalloc(F,d.nSites); 
  dalloc(F_1stRun,d.nSites);  
  dalloc(F_new,d.nSites);  
  dalloc(Q,nBoot+10); 
  dalloc(Q_new,nBoot+10); 
  
  dalloc(F_org,d.nSites);
  dalloc(F_orgOrg,d.nSites);
  
  delete [] N;
  delete [] bestLike;
 

  // problem with F and F_tmp being same pointer 
  // not there any more really

  d.idMap.clear();

  dallocBeagle(d);
  dallocBeagle(dOrg);
  dallocRefPanel(ref,nPop);


  for(std::vector<char*>::iterator it = tmp.begin(); it != tmp.end(); ++it){
    free(*it);
  } 
  tmp.clear();

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

