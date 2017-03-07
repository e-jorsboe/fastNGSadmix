/*
  log:
  g++ fastNGSadmix.cpp -lz -lpthread  -O3 -o fastNGSadmix

  log: (with readplink function)
  g++ fastNGSadmix.cpp readplinkV2.c -lz -lpthread  -O3 -o fastNGSadmix


  debug:
  g++ fastNGSadmix.cpp -lz -lpthread -ggdb -O3 -o fastNGSadmix

  debug: (with readplink function)
  g++ fastNGSadmix.cpp readplinkV2.c -lz -lpthread -ggdb -O3 -o fastNGSadmix

*/


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
#include <iostream>
// stringApocalypse
#include <string>
#include "readplinkV2.h"


//this is the max number of bytes perline
#define LENS 100000
//if we catch signal then quit program nicely
int SIG_COND =1;

double errTolMin=1e-5;
double errTolStart=0.05;
//frequencies and admixture coef cannot be less than this or more than 1-this
double errTol=errTolStart;


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


void minus1d(double *fst,double *sec,size_t x,double *res){
  for(size_t i=0;i<x;i++){
      res[i] = fst[i]-sec[i];
    }
}

void minus(double **fst,double **sec,size_t x,size_t y,double **res){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++){
      res[i][j] = fst[i][j]-sec[i][j];
    }
  }
}

double sumSquare(double **mat,size_t x,size_t y){
  double tmp=0;
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++){
      tmp += mat[i][j]*mat[i][j];
    }
  }
  return tmp;
}

double sumSquare1d(double *mat,size_t x){
  double tmp=0;
  for(size_t i=0;i<x;i++){
      tmp += mat[i]*mat[i];
  }
  return tmp;
}

double calcThres(double *d1,double *d2, int x){
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

 

// function for keeping sure Q values do not become
// 0.0 as then division by 0 might occur, errTol is limit
void map2domainQ(double* &Q, int nPop){  
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
void map2domainF(double** &F, int nSites, int nPop){
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

///@param str Filename given as a string.
int fexists(const char* str){
  struct stat buffer ;
  /// @return Function returns 1 if file exists.
  return (stat(str, &buffer )==0 ); 
}

std::vector<char *> dumpedFiles;
FILE *openFile(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  if(0&&fexists(c)){
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  dumpedFiles.push_back(strdup(c));

  FILE *fp = fopen(c,"w");
  if (not(fp)){
    fprintf(stderr,"File: %s cannot be created specify valid path\n",c);
    fflush(stderr);
    exit(0);
  } 
  
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
  if(0&&fexists(c)){
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  dumpedFiles.push_back(strdup(c));
  gzFile fp = gzopen(c,"w");
  delete [] c;
  return fp;
}






//some struct with all the data from the beagle file
typedef struct{
  double **genos;
  char *major;
  char *minor;
  // for snp ids, chr_pos
  char **ids;
  int nSites;
  int nInd;
  // map of ids in beagle file for finding overlap with ref
  // id is like this: chr_pos 
  std::map <std::string,int> idMap;
  
}bgl;


bgl allocBeagle(int nSites){
  bgl b;
  b.nSites = nSites;
  b.major = new char[nSites];
  b.minor = new char[nSites];
  b.ids = new char*[nSites];
  b.nInd = 1;
  b.genos= allocDouble(nSites,3);
  
  for(int s=0;SIG_COND&& (s<nSites);s++){
    b.major[s] = '0';
    b.minor[s] = '0';
  
  }
  return(b);
} 
      
//utility function for cleaning up out datastruct
void dallocBeagle(bgl &b){
  for(int i=0;i<b.nSites;i++){
    delete [] b.genos[i];  
    free(b.ids[i]);
  }
  
  // run through all the keys and then free them
  delete [] b.minor;
  delete [] b.major;
  delete [] b.genos;
  delete [] b.ids;
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
  // is coded so key is old column number, from inputted ref
  // value is new column (in ref for analysis) number + 1 (cause has to be above 0)
  std::map <int,int> colsToKeep;
  std::map <std::string,int> popsToKeep;
}refPanel;
 

void dallocRefPanel(refPanel &ref, int nPop){
  for(int i=0;i<ref.refSites;i++){
    delete [] ref.freqs[i];
    free(ref.id[i]);
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


/*
  Returns the bgl struct containing all data from a beagle file.
  
  It find the nsamples from counting the header
  It finds the number of sites by queing every line in a std::vector
  After the file has been read intotal it reloops over the lines in the vector and parses data
*/


bgl readBeagle(const char* fname, std::map <std::string,int> overlap) {
  const char *delims = "\t \n";
  gzFile fp = NULL;
  // checking if file can be opened
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  bgl ret;
  char buf[LENS];
  gzgets(fp,buf,LENS);
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

  if(ncols > 6){
    fprintf(stderr,"Only one individual in beagle file, looks like there are=%d\n",(ncols-3) / 3);
    exit(0);
  }
  ret.nInd = (ncols-3)/3;
  std::vector<char*> tmp;
  while(gzgets(fp,buf,LENS)){
    tmp.push_back(strdup(buf));
    // puts id of all sites in map for fast lookup
    char* lol = strtok(buf,delims);
    std::string myString(lol, strlen(lol));
    ret.idMap[myString] = 1;
  }
  
  ret.nSites = overlap.size();
  ret.major = new char[overlap.size()];
  ret.minor = new char[overlap.size()];
  ret.ids = new char*[overlap.size()];
  ret.genos= new double*[overlap.size()];
  //then loop over the vector and parsing every line
  int  refIndex = 0;
  for(int s=0;SIG_COND&& (s<tmp.size());s++){
    
    char * id = strtok(tmp[s],delims);
    std::string idString(id, strlen(id));

    if(overlap.count(idString)==0){
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

  
  for(int s=0;s<tmp.size();s++){
    free(tmp[s]);
  }
  //clean up filepointer
  gzclose(fp); 
  return ret;
}



// read in plink file and converts to a bgl struct (beagle file)
bgl readPlinkToBeagle(const char* plinkName, std::map <std::string,int> overlap) {

  plink* pl = readplink(plinkName);
  if(pl->fam.individuals > 1){
    fprintf(stderr,"More than one individual in input plink file - should only be one! \n");
    exit(0);
  }
  bgl b = allocBeagle(overlap.size());
  b.nSites = overlap.size();
  // can only have one individual is this program
  b.nInd = 1; 
  int beagleIndex = 0;
  for(int i=0;i<(pl->y);i++){
    // put second plink allele first in bgl
    // and then first plink allele
    // then put id as chr_pos
    std::string bimString(pl->bim.id[i],strlen(pl->bim.id[i]));
    if(overlap.count(bimString)<1){

      continue;
    }
    if(beagleIndex+1 > overlap.size()){
      fprintf(stderr,"Duplicated sites - multiple sites with same position in plink file!! \n");
      fprintf(stderr,"Use plink2 with --list-duplicate-vars suppress-first ids-only - and then --exclude \n");
      exit(0);
    }
    
    if(pl->d[0][i]==0){
      b.genos[beagleIndex][2]=1.0     ;
      b.genos[beagleIndex][1]=0.0;
      b.genos[beagleIndex][0]=0.0;
      
    } else if(pl->d[0][i]==1){
      b.genos[beagleIndex][2]=0.0;
      b.genos[beagleIndex][0]=0.0;
      b.genos[beagleIndex][1]=1.0;
    } else if(pl->d[0][i]==2){
      b.genos[beagleIndex][2]=0.0;
      b.genos[beagleIndex][1]=0.0;
      b.genos[beagleIndex][0]=1.0;
    }
    b.major[beagleIndex]=pl->bim.major[i];
    b.minor[beagleIndex]=pl->bim.minor[i];
    b.ids[beagleIndex]=strdup(pl->bim.id[i]);
    // stores id of all sites from plink file
    b.idMap[bimString] = 1;

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


// function for reading in ref panel, has to have format id chr pos name A0_freq A1 pop1 pop2 ... (freq has to be of A0 allele)
std::vector<char *> tmpRef;
refPanel readRefPanel(const char* fname, bgl &b, std::map <std::string,int> &includedPops, int nPop, std::map <std::string,int> &overlap) {
  const char *delims = "\t \n";
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  char buf[LENS];
  refPanel ref;
  //find number of columns
  while(gzgets(fp,buf,LENS)){
    tmpRef.push_back(strdup(buf));    
  }
  // copy first line into buffer, as empty now on some systems
  // copies the values of num bytes from the location pointed to by source (2nd) directly to the memory block pointed to by destination (1st).
  memcpy(buf,tmpRef[0],strlen(tmpRef[0]));  
  int totalSites = tmpRef.size();
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
  if(nPop>0){
    ref.populations = new char*[nPop];
  } else{
    ref.populations = new char*[ncols-6];
  }

  // keep track of which original column it is
  int orgCol = 0;
  // keeps track of which new column (index = newCol-1) has to be above 0 for lookup in map
  int newCol = 1;
  char* columnID = strtok(tmpRef[0],delims);
  while(columnID!=NULL){
    orgCol++;
    // first 6 columns not freqs and has to at most K new columns included in ref
    if(orgCol>6){
      // if in supplied populations or if no populations supplied include
      if(includedPops.count(columnID)>0 or includedPops.empty()){
	// prints out which populations chosen
	fprintf(stderr,"Chosen pop %s\n",columnID);
	ref.populations[newCol-1] = strdup(columnID);
 	ref.popsToKeep[columnID] = newCol;
	// so that can translate from org column (where 7th column is first freq column) to new column
	ref.colsToKeep[orgCol-7] = newCol;
	newCol++;
      } 
    }
    columnID = strtok(NULL,delims);
   }
  

  ref.id = new char*[b.nSites];
  ref.chr = new int[b.nSites];
  ref.pos = new int[b.nSites];
  ref.name = new char*[b.nSites];
  ref.A0 = new char[b.nSites];
  ref.A1 = new char[b.nSites];
  ref.freqs = new double*[b.nSites];
  ref.refSites = b.nSites;
  // for keeping track of which index in new ref with
  int refIndex = 0;
  for(int s=0;SIG_COND&& (s<totalSites);s++){

    // looking at id value chr_pos for detecting overlap
    char* id = strtok(tmpRef[s],delims);
    std::string stringID(id,strlen(id));
    // check if site is in overlap with beagle file
    // otherwise continues to next site in ref
    if(overlap.count(stringID) > 0){
      ref.id[refIndex] = strdup(id);    
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
	// if bgl 1_1 A B GL(AA) GL(AB) GL(BB) Then ref 1 1 rs1 B A 1-f(B)
	// minor is last allele in beagle file
	if(ref.A0[refIndex]==b.minor[refIndex]){
	  // new col has to be - 1 for right index
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
    // only here if site was included in new ref
    refIndex++;       
  }

  if(refIndex!=overlap.size()){
    fprintf(stderr,"refIndex %i, overlap sites %lu\n",refIndex,overlap.size());
    fprintf(stderr,"Seems like there are duplicate ids in input or reference panel!\n");
    exit(0);
  }
  //clean up filepointer
  gzclose(fp); 
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
  while(gzgets(fp,buf,LENS)){
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

// for reading number of individuals in each ref - nInd file
void readDouble1d(double *d,int nPop,const char*fname, std::map<std::string,int> popsToKeep){
  fprintf(stderr,"Opening nInd file: %s with K=%d\n",fname,nPop);
  const char*delims=" \n\t";
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
   
  std::vector<char*> tmp;
  // keeps track of org index of pop to keep and new index of pop to keep
  std::map<int,int> toKeep;
  while(fgets(buf,lens,fp)){
    tmp.push_back(strdup(buf));
  }

  // first goes through the first line with names of pops, finds which should be included
  char * word = strtok(tmp[0],delims);
  // keeps track of which value at in nInd file, newCol has to be index+1, because has to be > 0 for lookup
  int orgCol = 0;

  while(word!=NULL){
    std::string stringWord(word,strlen(word));
    if(popsToKeep.count(stringWord)>0){
    
      // creates map of <nInd index, ref index> so can map from nInd order to ref order
      // so if first element in nInd is second in ref
      // the d array with nInd values has second element equal to nInd first value
      toKeep[orgCol] = popsToKeep[word];

    }
    word = strtok(NULL,delims); 
    orgCol++;
  }
  if(toKeep.size()!=popsToKeep.size()){
    fprintf(stderr,"nInd and ref panel do not have same size!\n");
    exit(0);
  }
  free(tmp[0]);
  // secondly goes through the second line with sizes of pops
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

  if(index!=orgCol){
    fprintf(stderr,"nInd has different number of elements between row 1 and 2\n");
    exit(0);
  }
  

 
  free(tmp[1]);
  fclose(fp);
}

void printDouble(double **ret,size_t x,size_t y, int highestLike, int nConv, char** populations ,FILE *fp){
  for(size_t i=0;i<x;i++){
    if(i==0){
      for(size_t j=0;j<y;j++){
	fprintf(fp,"%s ",populations[j]);
      }
      fprintf(fp,"\n");
    }
    if(i<nConv and i == highestLike){
      for(size_t j=0;j<y;j++){
	fprintf(fp,"%.4f ",ret[i][j]);
      }
      fprintf(fp,"\n");
    } else if(i>=nConv){
      for(size_t j=0;j<y;j++){
	fprintf(fp,"%.4f ",ret[i][j]);
      }
      fprintf(fp,"\n");
    }
  }
  
}


void printDoubleGz(double **ret,size_t x, size_t y,  char** id, char** populations ,gzFile fp){

 for(size_t i=0;i<x;i++){
    if(i==0){
      gzprintf(fp,"id ");
      for(size_t j=0;j<y;j++){
	gzprintf(fp,"%s ",populations[j]);
      }
      gzprintf(fp,"\n");
    }
    gzprintf(fp,"%s ",id[i]);
    for(size_t j=0;j<y;j++){
      gzprintf(fp,"%.4f ",1-ret[i][j]);
    }
    gzprintf(fp,"\n");
  }
}

// calculate log(likelihood) from likelihood function
double likelihood(double* Q, double** F,int nSites, int nPop,double **genos){
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

// does bootstrapping sampling nSites random sites with replacement
void bootstrap(bgl dOrg, bgl &d, double** F_orgOrg, double** F_org, double** F, int nPop) {
  for(int j=0;j<dOrg.nSites;j++){
    // generate random int from 0 to (nSites-1)
    int row = std::rand() % dOrg.nSites;
    for(int k = 0; k < nPop; k++) {
      F[j][k] = F_orgOrg[row][k];
      F_org[j][k] = F_orgOrg[row][k];
      if(k<=2){
	d.genos[j][k] = dOrg.genos[row][k];
      }
    }
  }
}

// em algorithm not adjusting F at every step
void emUnadjusted(double* Q, double** F, int nSites, int nPop,double **genos,double *Q_1) {
  double sumAG[nPop];
  double sumBG[nPop];
  // makes sure neither F nor Q has 0 values
  map2domainF(F, nSites, nPop);
  map2domainQ(Q,nPop);
  for(int k=0;k<nPop;k++){ 
    sumAG[k]=0;
    sumBG[k]=0;
  }
  for(int j=0;j<nSites;j++){   
    double fpart=0;
    double fpartInv=0;
    double expGG=0;
    for(int k=0;k<nPop;k++){ 
      // admixture adjusted freq, for each pop
      fpart += F[j][k] * Q[k];
      fpartInv += (1-F[j][k]) * Q[k];
      // pre GL (sites x 3) * (adjusted freq)
      double pp0=(fpartInv)*(fpartInv)*genos[j][2];
      double pp1=2*(fpartInv)*fpart*  genos[j][1];
      double pp2=fpart*fpart*        genos[j][0];
      double sum=pp0+pp1+pp2;
      // for calculating H range 0-2, this is the expected genotype	  
      expGG =(pp1+2*pp2)/sum;
    }
    for(int k=0;k<nPop;k++){ 
      sumAG[k] += expGG/(fpart) * (Q[k] * F[j][k]); 
      sumBG[k] += (2-expGG)/fpartInv * (Q[k] * (1-F[j][k])); 
    }
  }
  for(int k=0;k<nPop;k++){ 
    Q_1[k]=(sumAG[k] + sumBG[k])/(2.0*nSites);
  }
  map2domainQ(Q_1,nPop);
}

// em algorithm adjusting F at every step
void em(double* Q, double** F, int nSites, double* nInd, int nPop,double **genos,double **F_1,double *Q_1, double** F_org) {
  double sumAG[nPop];
  double sumBG[nPop];
  // makes sure neither F nor Q has 0 values
  map2domainF(F, nSites, nPop);
  map2domainF(F_org, nSites, nPop);
  map2domainQ(Q,nPop);
  double sumA[nPop];
  double sumB[nPop];
  for(int k=0;k<nPop;k++){ 
      sumA[k]=0;
      sumB[k]=0;
  }
  for(int j=0;j<nSites;j++){   
    for(int k=0;k<nPop;k++){ 
      sumAG[k]=0;
      sumBG[k]=0;
    }
    double fpart=0;
    double fpartInv=0;
    double expGG=0;
    for(int k=0;k<nPop;k++){ 
      // admixture adjusted freq, for each pop
      fpart += F[j][k] * Q[k];
      fpartInv += (1-F[j][k]) * Q[k];
      // pre GL (sites x 3) * (adjusted freq)
      double pp0=(fpartInv)*(fpartInv)*genos[j][2];
      double pp1=2*(fpartInv)*fpart*  genos[j][1];
      double pp2=fpart*fpart*        genos[j][0];
      double sum=pp0+pp1+pp2;
      // for calculating H range 0-2, this is the expected genotype	
      expGG=(pp1+2*pp2)/sum;
    }
    for(int k=0;k<nPop;k++){
      // similar to (H/(q*f))*q, for jth marker
      sumAG[k] = (expGG) / (fpart) * (Q[k]*F[j][k]); 
      sumBG[k] = (2-expGG) / fpartInv * (Q[k]*(1-F[j][k]));
      sumA[k] += sumAG[k];
      sumB[k] += sumBG[k];
      sumAG[k] += nInd[k]*2*F_org[j][k];
      sumBG[k] += 2*nInd[k]-(2*nInd[k]*F_org[j][k]);      
    }
    for(int k=0;k<nPop;k++){
      // adjust with ref panel, so we have input + ref expected number of alleles
      F_1[j][k]=sumAG[k]/(sumAG[k]+sumBG[k]);
    }
    
  }
  for(int k=0;k<nPop;k++){ 
    Q_1[k]=(sumA[k]+sumB[k])/(2.0*nSites);
  }
  map2domainQ(Q_1,nPop);
  map2domainF(F_1,nSites,nPop);
}


int emAccelUnadjustedV2(const bgl &d, double* nInd, int nPop, double** F,double* Q,double* &Q_new,double &lold, double** F_org, int nit, int boot, int Qconv, double Qtol, double tol){
 
  double stepMin = 1;
  double stepMax0 = 1;
  static double stepMax = stepMax0;
  double mstep = 4;
  double objfnInc = 1;
  //we make these huge structures static such that we just allocate them the first time
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
  minus1d(Q_em1,Q,nPop,Q_diff1);
  double sr2 = sumSquare1d(Q_diff1,nPop);
  // checks if convergence
  if(sqrt(sr2)<tol or (calcThres(Q,Q_em1,nPop) < Qtol and Qconv>0)){
    //fprintf(stderr,"like is %f\n",likelihood(Q_new, F, d.nSites, nPop,d.genos));    
    return 0;  
  }
  // second EM run
  emUnadjusted(Q_em1, F, d.nSites, nPop, d.genos, Q_em2);
  minus1d(Q_em2,Q_em1,nPop,Q_diff2);
  double sq2 = sumSquare1d(Q_diff2,nPop);
  // checks if convergence - a second time
  if(sqrt(sq2)<tol  or (calcThres(Q_em1,Q_em2,nPop) < Qtol and Qconv>0)){
    //fprintf(stderr,"like is %f\n",likelihood(Q_new, F, d.nSites, nPop,d.genos));
    return 0;
  }
  minus1d(Q_diff2,Q_diff1,nPop,Q_diff3);
  double sv2 = sumSquare1d(Q_diff3,nPop);
  double alpha = sqrt(sr2/sv2);
  // makes sure alpha does not go below 1 and above stepMax
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<nPop;i++){
    Q_new[i] = Q[i]+2*alpha*Q_diff1[i]+alpha*alpha*Q_diff3[i];
    Q_tmp[i] = 1.0;
  }
  map2domainQ(Q_new,nPop);
  // if alpha not too close (0.01 close) to 1 
  if (fabs(alpha - 1) > 0.01){
    // we estimate new Q and F, with our inferred Q and F via alpha
    emUnadjusted(Q_new, F, d.nSites, nPop,d.genos,Q_tmp);    
    minus1d(Q_tmp,Q_new,nPop,Q_tmpDiff);
    // adopted from squarem2 from SQUAREM package
    double res = sumSquare1d(Q_tmpDiff,nPop);
    double parnorm = (1/std::sqrt(nPop))*sumSquare1d(Q_tmpDiff,nPop);
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
      double lnew = likelihood(Q_new, F, d.nSites, nPop,d.genos);
      fprintf(stderr,"iter[%d] like=%f alpha=%f ",nit,lnew,alpha);
      for(int i=0;i<nPop;i++){	      
	fprintf(stderr,"Q=%f, ",Q_new[i]);
      }
      fprintf(stderr,"\n");
    }
  }
  return 1;
}

// based on squarem1, from SQUAREM R package, by RAVI VARADHAN and CHRISTOPHE ROLAND Scandinavian Journal of Statistics, Vol. 35: 335–353, 2008
int emAccelV2(const bgl &d, double* nInd, int nPop, double** F,double* Q,double** &F_new,double* &Q_new,double &lold, double** F_org, int nit, int boot, int Qconv, double Qtol, double tol){

  double stepMin = 1;
  double stepMax0 = 1;
  static double stepMax = stepMax0;
  double mstep = 4;
  double objfnInc = 1;

  //we make these huge structures static such that we just allocate them the first time
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
  em(Q, F, d.nSites, nInd, nPop,d.genos, F_em1, Q_em1, F_org);
  minus(F_em1,F,d.nSites,nPop,F_diff1);
  minus1d(Q_em1,Q,nPop,Q_diff1);
  double sr2 = sumSquare1d(Q_diff1,nPop) + sumSquare(F_diff1,d.nSites,nPop);
  // checks if convergence
  if(sqrt(sr2)<tol or (calcThres(Q,Q_em1,nPop) < Qtol and Qconv>0)){
    //fprintf(stderr,"like is %f\n",likelihood(Q_new, F_new, d.nSites, nPop,d.genos));    
    return 0;
  }
  // second EM run
  em(Q_em1, F_em1, d.nSites, nInd, nPop,d.genos, F_em2, Q_em2, F_org);
  minus(F_em2,F_em1,d.nSites,nPop,F_diff2);
  minus1d(Q_em2,Q_em1,nPop,Q_diff2);
  double sq2 = sumSquare1d(Q_diff2,nPop) + sumSquare(F_diff2,d.nSites,nPop);
  // checks if convergence - a second time
  if(sqrt(sq2)<tol or (calcThres(Q_em1,Q_em2,nPop) < Qtol and Qconv>0)){
    //fprintf(stderr,"like is %f\n",likelihood(Q_new, F_new, d.nSites, nPop,d.genos));
    return 0;
  }
  minus(F_diff2,F_diff1,d.nSites,nPop,F_diff3);
  minus1d(Q_diff2,Q_diff1,nPop,Q_diff3);
  double sv2 = sumSquare1d(Q_diff3,nPop) + sumSquare(F_diff3,d.nSites,nPop);
  double alpha = sqrt(sr2/sv2);  
  // makes sure alpha does not go below 1 and above stepMax
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<d.nSites;i++){
    for(size_t j=0;j<nPop;j++){
      // based on ngsAdmix approach
      F_new[i][j] = F[i][j]+2*alpha*F_diff1[i][j]+alpha*alpha*F_diff3[i][j];
      F_tmp[i][j] = 1.0;
    }
  }
  map2domainF(F_new,d.nSites,nPop);
  for(size_t i=0;i<nPop;i++){
    Q_new[i] = Q[i]+2*alpha*Q_diff1[i]+alpha*alpha*Q_diff3[i];
    Q_tmp[i] = 1.0;
  }
  map2domainQ(Q_new,nPop);
  // if alpha not too close (0.01 close) to 1 
  if (fabs(alpha - 1) > 0.01){
    // we estimate new Q and F, with our inferred Q and F via alpha
    em(Q_new, F_new, d.nSites, nInd, nPop,d.genos,F_tmp,Q_tmp,F_org);
    minus(F_tmp,F_new,d.nSites,nPop,F_tmpDiff);
    minus1d(Q_tmp,Q_new,nPop,Q_tmpDiff);
    double res = sumSquare1d(Q_tmpDiff,nPop) + sumSquare(F_tmpDiff,d.nSites,nPop);
    double parnorm = (1/std::sqrt(nPop))*sumSquare1d(Q_tmpDiff,nPop) + (1/std::sqrt(d.nSites*nPop))*sumSquare(F_tmpDiff,d.nSites,nPop);
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
      double lnew = likelihood(Q_new, F_new, d.nSites, nPop,d.genos);
      fprintf(stderr,"iter[%d] like=%f alpha=%f ",nit,lnew,alpha);
      for(int i=0;i<nPop;i++){	      
	fprintf(stderr,"Q=%f, ",Q_new[i]);
      }
      fprintf(stderr,"\n");
    }
  }
  
  return 1;
}

// finds overlap between input and ref panel based on id chr_pos, assumes no duplicate sites in terms of position
std::map <std::string,int> findOverlap(const char* lname, const char* plinkName, const char* fname, FILE* flog){
  std::map <std::string,int> inputSites;
  static std::map <std::string,int> overlap;
  const char *delims = "\t \n";
  std::vector<char*> tmpBeagle;
  std::vector<char*> tmpRef;
  // if beagle file input
  if(plinkName==NULL){
    gzFile fp1 = NULL;
    if(Z_NULL==(fp1=gzopen(lname,"r"))){
      fprintf(stderr,"Error opening file: %s\n",lname);
      exit(0);
    }
    char buf1[LENS];
    while(gzgets(fp1,buf1,LENS)){
      tmpBeagle.push_back(strdup(buf1));
    }
    // reads all beagle sites into map, checks for duplicates!
    for(int s=1;SIG_COND&& (s<tmpBeagle.size());s++){
      char* bglID = strtok(tmpBeagle[s],delims);
      std::string bglString(bglID,strlen(bglID));
      if(inputSites.count(bglString)>0){
	fprintf(stderr,"Duplicate sites in beagle file: %s - Go fix!\n",bglID);
	exit(0);
      } else{
	inputSites[bglString]=1;
      }
    }    
    gzclose(fp1);
    // if plink file input
  } else{
    plink* pl_tmp = readplink(plinkName);
    // reads all plink sites into map, checks for duplicates!
    for(int s=0;s<(pl_tmp->y);s++){
      std::string plinkString(pl_tmp->bim.id[s],strlen(pl_tmp->bim.id[s])); 
      if(inputSites.count(plinkString)>0){
	fprintf(stderr,"Duplicate sites in plink file: %s - Go fix!\n",pl_tmp->bim.id[s]);
	exit(0);
      }
      else if(pl_tmp->d[0][s]==3){
	continue;
      }
      inputSites[plinkString]=1;  
    }
    kill_plink(pl_tmp);
  }

  fprintf(stderr,"Input has this many sites %zu\n",inputSites.size());
  fprintf(flog,"Input has this many sites %zu\n",inputSites.size());
  // reads ref panel
  gzFile fp2 = NULL;
  char buf2[LENS];
  if(Z_NULL==(fp2=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  while(gzgets(fp2,buf2,LENS)){
    tmpRef.push_back(strdup(buf2));    
  }
  fprintf(stderr,"Ref has this many sites %zu\n",tmpRef.size());
  fprintf(flog,"Ref has this many sites %zu\n",tmpRef.size());
  // starts at 1 to avoid header
  for(int s=1;SIG_COND&& (s<tmpRef.size());s++){
    // looking at id value chr_pos for detecting overlap
    char* id = strtok(tmpRef[s],delims);
    std::string plinkString(id,strlen(id)); 
    // check if site is in beagle or plink file
    // otherwise continues to next site in ref
    if(overlap.count(plinkString)>0){
      fprintf(stderr,"Duplicate site in ref panel: %s - Go fix!\n",id);     
      exit(0);
    } else if(inputSites.count(plinkString) > 0){
      overlap[plinkString] = 1;    
    }	else {
      continue;
    }
  }
  if(overlap.size()==0){
    fprintf(stderr,"No overlapping sites where found!!\n");
    exit(0);
  }
  gzclose(fp2);
  // cleaning
  if(lname!=NULL){
    for(int s=0;SIG_COND&& (s<tmpRef.size());s++){
      free(tmpRef[s]);
    }
    for(int s=0;SIG_COND&& (s<tmpBeagle.size());s++){
      free(tmpBeagle[s]);
    }
  } else if(plinkName!=NULL){
    for(int s=0;SIG_COND&& (s<tmpRef.size());s++){
      free(tmpRef[s]);
    }
  } 
  return(overlap);
}

void info(){
  
  fprintf(stderr,"Arguments:\n");
  fprintf(stderr,"\t-likes Beagle likelihood filename\n");
  fprintf(stderr,"\t-plink Plink file in the binary bed format\n");
  fprintf(stderr,"\t-Nname Number of individuals in each reference populations\n");
  fprintf(stderr,"\t-fname Ancestral population frequencies\n");

  fprintf(stderr,"Optional:\n");
  fprintf(stderr,"\t-out Prefix for output files\n"); 
  fprintf(stderr,"\t-printFreq print admixture adjusted allele frequencies of reference panel + input individual (1: yes, 0: no (default))\n"); 

  fprintf(stderr,"Setup:\n");
  fprintf(stderr,"\t-whichPops Which populations from the ref panel to include in analysis, denotes number of populations (nPop) for admixture estimation\n \t if not specified all populations in ref are analyzed, must be comma seperated (pop1,pop2,..)\n");
  fprintf(stderr,"\t-doAdjust Adjusts the frequencies in the reference populations with the input (1: yes (default), 0: no)\n");
  fprintf(stderr,"\t-seed Seed for initial guess in EM and for bootstrap\n"); 
  fprintf(stderr,"\t-method If 0 no acceleration of EM algorithm (1: yes (default), 0: no)\n"); 

  fprintf(stderr,"Stop chriteria:\n"); 
  fprintf(stderr,"\t-Qconv Stopping criteria based on change in Q (works best when using doAdjust) (1: yes, 0: no (default))\n"); 
  fprintf(stderr,"\t-Qtol Tolerance value for stopping criteria based on change in Q (0.001 (default))\n"); 
  fprintf(stderr,"\t-tol Tolerance for convergence\n"); 
  fprintf(stderr,"\t-maxiter Maximum number of EM iterations\n"); 
  fprintf(stderr,"\t-boot Number of bootstrapping iterations, default 0, can at most be 10000, .qopt FIRST row BEST estimated Q, rest bootstraps!!\n"); 
  fprintf(stderr,"\t-conv Number of convergence iterations, each with random starting point, to check if has converged, default 1, can at most be 10\n"); 

  exit(0);
}


// for handling Ctrl + C
int VERBOSE =1;
void handler(int s) {
  if(VERBOSE)
    fprintf(stderr,"Caught SIGNAL: %d will try to exit nicely (no more threads are created, we will wait for the current threads to finish)\n",s);
  VERBOSE=0;
  SIG_COND=0;
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
  int printFreq = 0;
  int doAdjust = 1;
  float minMaf = 0.00;
  const char* lname = NULL;
  const char* fname = NULL;
  char* pops = NULL;
  const char* Nname = NULL;
  const char* outfiles = NULL;
  const char* plinkName = NULL;
  int nPop = 0;
  int seed = time(NULL);
  double tol=0.00001; 
  int nBoot = 0; 
  int nConv = 1;
  int Qconv = 0;
  double Qtol = 0.0000001;
  
  // reading arguments
  argv++;
  while(*argv){
    // GL in the shape of beagle file
    if(strcmp(*argv,"-likes")==0 || strcmp(*argv,"-l")==0) lname=*++argv; //name / char arrays
  
    else if(strcmp(*argv,"-plink")==0 || strcmp(*argv,"-p")==0) plinkName=*++argv;
    // ref panel
    else if(strcmp(*argv,"-fname")==0 || strcmp(*argv,"-f")==0) fname=*++argv; 
    // nInd file
    else if(strcmp(*argv,"-Nname")==0 || strcmp(*argv,"-N")==0) Nname=*++argv;
    // prefix for output files
    else if(strcmp(*argv,"-outfiles")==0 || strcmp(*argv,"-out")==0) outfiles=*++argv;
    // which populations in ref panel to be ananlyzed, must agree with K
    else if(strcmp(*argv,"-whichPops")==0 || strcmp(*argv,"-pops")==0) pops=*++argv; 
    else if(strcmp(*argv,"-seed")==0||strcmp(*argv,"-s")==0) seed=atoi(*++argv); //int - atoi - char array to integer
    // flag for printing adjusted freqs
    else if(strcmp(*argv,"-printFreq")==0) printFreq=atoi(*++argv); 
    // flag for doing Adjustment
    else if(strcmp(*argv,"-doAdjust")==0) doAdjust=atoi(*++argv);
    // flag for doing accelerated EM
    else if(strcmp(*argv,"-method")==0 || strcmp(*argv,"-m")==0) method=atoi(*++argv); 
    // different stop criteria - whether based on diff in Q values
    else if(strcmp(*argv,"-Qconv")==0) Qconv=atoi(*++argv); 
    // tolerance for Q stopping criteria
    else if(strcmp(*argv,"-Qtol")==0) Qtol=atof(*++argv);
    // tolerance for likelihood based stopping criteria
    else if(strcmp(*argv,"-tol")==0) tol=atof(*++argv);
    // number of boot straps
    else if(strcmp(*argv,"-bootstrap")==0||strcmp(*argv,"-boot")==0) nBoot=atoi(*++argv);
    // number of convergence runs with different starting points
    else if(strcmp(*argv,"-convergenceRuns")==0||strcmp(*argv,"-conv")==0) nConv=atoi(*++argv);
    // number of max total iterations
    else if(strcmp(*argv,"-maxiter")==0 || strcmp(*argv,"-i")==0) maxIter=atoi(*++argv); 
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      info();
      return 0;
    }
    ++argv;
  }

  //check that non optional options have been used. 
  if(lname==NULL and plinkName==NULL){
    fprintf(stderr,"Please supply a beagle or plink input file: -likes or -plink\n");
    info();
  } else if(fname==NULL){
    fprintf(stderr,"Please supply a reference panel: -fname\n");
    info();
  } else if(lname!=NULL and plinkName!=NULL){
    fprintf(stderr,"Please supply ONLY a beagle or plink input file, not BOTH: -likes or -plink\n");
    info();
  } else if(Nname==NULL){
    fprintf(stderr,"Please supply number of individauls file: -Nname\n");
    info();
  } 

  if(outfiles==NULL and lname!=NULL){
    fprintf(stderr,"Will use beagle name as prefix for output\n");
    outfiles=lname;
  }
  
  if(outfiles==NULL and plinkName!=NULL){
    fprintf(stderr,"Will use plink name as prefix for output\n");
    outfiles=plinkName;
  }
  // max 10000 bootstraps 10 conv runs and 
  nBoot = std::min(std::max(nBoot,0),10000);
  nConv = std::min(std::max(nConv,1),10);
  
  //out put files
  FILE *flog=openFile(outfiles,".log");

  fprintf(stderr,"Input: likes=%s plink=%s Nname=%s fname=%s outfiles=%s\n",lname,plinkName,Nname,fname,outfiles);
  fprintf(stderr,"Setup: seed=%d method=%d\n",seed,method);
  if(method==0){
    fprintf(stderr,"The unaccelerated EM has been chosen\n");
  } else{
    fprintf(stderr,"The accelerated EM has been chosen\n");
    tol=1e-7; //stopping criteria
  }
  if(doAdjust==0){
    fprintf(stderr,"The unadjusted method has been chosen\n");
  } else{
    fprintf(stderr,"The adjusted method has been chosen\n");
  }
  fprintf(stderr,"Convergence: maxIter=%d tol=%.8f\n",maxIter,tol);
  fprintf(stderr,"The following number of bootstraps have been chosen: %i\n",nBoot);
  if(Qconv>0){
    fprintf(stderr,"Convergence via difference in Q values chosen, threshold of: %f\n",Qtol);    
  }

  fprintf(flog,"Input: likes=%s plink=%s Nname=%s fname=%s outfiles=%s\n",lname,plinkName,Nname,fname,outfiles);
  fprintf(flog,"Setup: seed=%d method=%d\n",seed,method);
  if(method==0){
    fprintf(flog,"The unaccelerated EM has been chosen\n");
  } else{
    fprintf(flog,"The accelerated EM has been chosen\n");
  }
  if(doAdjust==0){
    fprintf(flog,"The unadjusted method has been chosen\n");
  } else{
    fprintf(flog,"The adjusted method has been chosen\n");
  }
  fprintf(flog,"Convergence: maxIter=%d tol=%.8f\n",maxIter,tol);
  fprintf(flog,"The following number of bootstraps have been chosen: %i\n",nBoot);
  if(Qconv>0){
    fprintf(flog,"Convergence via difference in Q values chosen, threshold of: %f\n",Qtol);    
  }


  // to get the populations from ref to be analyzed
  std::map <std::string,int> includedPops;
  
  // if pops are specified, reads which pops and construcs map with those
  if(pops!=NULL){
    char* temp = strtok(pops,",");
    while(temp!=NULL){
      nPop++;
      std::string popString(temp,strlen(temp));
      // check that population does not appear twice here
      if(includedPops.count(popString)>0){
	fprintf(stderr,"Same population selected twice with -whichPops, only each population once!\n");
	fprintf(flog,"Same population selected twice with -whichPops, only each population once!\n");
	info();
      }
      includedPops[popString] = 1;
      temp = strtok(NULL,",");
    }   
  }

  // to find out which version of C++
  //fprintf(stderr,"V: [%ld] ", __cplusplus);
  
  bgl d;
  bgl dOrg;
  std::map <std::string,int> overlap;
  if(lname!=NULL){
    // finds overlapping sites, then reads beagle
    overlap = findOverlap(lname, NULL, fname,flog);
    d=readBeagle(lname,overlap);
    dOrg=readBeagle(lname,overlap);
  } else if(plinkName!=NULL){
    // finds overlapping sites, then reads plink file
    overlap = findOverlap(NULL, plinkName, fname,flog);
    d=readPlinkToBeagle(plinkName, overlap);
    dOrg=readPlinkToBeagle(plinkName, overlap);  
  }
  
  fprintf(stderr,"Overlap: of %zu sites between input and ref\n",overlap.size());
  fprintf(flog,"Overlap: of %zu sites between input and ref\n",overlap.size());  
  
  
  // reads in ref Panel
  refPanel ref;

  // if nPop == 0 then read in all of them refs
  ref = readRefPanel(fname,d,includedPops,nPop,overlap);
  if(pops==NULL){
    nPop = ref.popsToKeep.size();
    // so checks that whichPops has same pops as ref, if whichPops given
  } else  if(includedPops.size()!=ref.popsToKeep.size()){
    fprintf(stderr,"Some populations given are not in the ref panel\n");
    info();
  }

  if(nPop<2){
    fprintf(stderr,"nPop has to be at least 2, nPop=%i\n",nPop);
    fprintf(flog,"nPop has to be at least 2, nPop=%i\n",nPop);
    info();
  }
  fprintf(stderr,"\n");
  fprintf(flog,"\n");
  fprintf(stderr,"nPop=%i\n",nPop);
  fprintf(flog,"nPop=%i\n",nPop);
  fprintf(stderr,"\n");
  fprintf(flog,"\n");
  // min tolerance for F and Q - so they are not 0
  errTolStart = errTolMin;
  errTol = errTolMin;
    
  clock_t t = clock();
  time_t t2 = time(NULL);

  // seed for bootstrapping and random starting points
  std::srand(seed);
      
  double **F = allocDouble(d.nSites,nPop);
  double **F_new = allocDouble(d.nSites,nPop);
  // F_org is for storing initial freqs from ref panel
  double **F_org = allocDouble(d.nSites,nPop);
  // for sampling from when doing bootstrap
  double **F_orgOrg = allocDouble(d.nSites,nPop);
  // for when having to write adjusted freqs 
  double **F_1stRun = allocDouble(d.nSites,nPop);  
  // gets freq values from ref struct, where read into
  for(int i=0;i<d.nSites;i++){
    for(int k=0;k<nPop;k++){
      F[i][k] = ref.freqs[i][k];
      F_org[i][k] = ref.freqs[i][k];
      F_orgOrg[i][k] = ref.freqs[i][k];
    }
     
  }

  // because has to have conv values for converge runs, and then nBoot bootstrapped values
  double **Q = allocDouble(nBoot+nConv,nPop);
  double **Q_new = allocDouble(nBoot+nConv,nPop);
  double *sum = new double[nBoot+nConv];
  // nBoot + conv rows
  for(int j=0;j<nBoot+nConv;j++){
    sum[j]=0;
    for(int k=0;k<nPop;k++){
      Q[j][k]=rand()*1.0/RAND_MAX;
      sum[j]+=Q[j][k];
    }
  }
  // nBoot + conv rows
  for(int j=0;j<nBoot+nConv;j++){
    for(int k=0;k<nPop;k++) {
      // to make sure that proportions sum to 1
      Q[j][k] = Q[j][k]/sum[j]; 
      Q_new[j][k] = Q[j][k];
    }
  }
  delete [] sum;
  
  double *N = new double[nPop];  
  // reading nInd, where colsToKeep from ref to read in the same columns as in ref
  readDouble1d(N,nPop,Nname,ref.popsToKeep);
  fprintf(flog,"Opening nInd file: %s with nPop=%d\n",fname,nPop); 
  for(int i=0;i<nPop;i++){
    fprintf(flog,"Chosen pop %s\n",ref.populations[i]);
    fprintf(flog,"N = %f\n",N[i]);
    fprintf(stderr,"N = %f\n",N[i]);
    // being printed in function to stderr

  }
  fprintf(stderr,"\n");
  fprintf(flog,"\n");

  double* bestLike = new double[nConv];
  int highestLike = 0;
  // initial likelihood   
  double lold = likelihood(Q[0], F_org, d.nSites, nPop,d.genos);
  fprintf(stderr,"iter[start] like is=%f\n",lold);
  int nit = 0;
  double likeLast = lold;
  double lastQthres = 0;

  //////////////////////////////////////// em ///////////////////////////////////
  
  //below is the main looping through the iterations.
  // we have 4 possible ways, unAdjusted/Adjusted basic EM/Accelerated EM
  // first conv runs for converge with new Q starting point
  // then nBoot new runs for bootstrapping with best Q starting point, with random sites
  for(int b=0;SIG_COND and b<(nBoot+nConv);b++) {
    // resets nit for each conv or bootstrap run
    nit=0;  
    if(b>(nConv-1)){
      // do bootstrapping when conv runs done
      bootstrap(dOrg,d,F_orgOrg,F_org,F,nPop); 
      fprintf(stderr,"At this bootstrapping: %i out of: %i\n",b-(nConv-1),nBoot);
      fprintf(flog,"At this bootstrapping: %i out of: %i\n",b-(nConv-1),nBoot);
    }
    for(nit=1;SIG_COND and nit<maxIter;nit++) {	
      if(doAdjust==0){
	if(method==0){
	  // unadjusted, unaccelerated EM
	  emUnadjusted(Q[b], F, d.nSites, nPop,d.genos,Q_new[b]);
	} else{
	  // unadjusted, accelerated EM
	  if(0==emAccelUnadjustedV2(d, N, nPop, F, Q[b], Q_new[b],lold,F_org, nit, b, Qconv, Qtol, tol)){
	    if(b<nConv){
	      // stores all likelihoods so max can be found
	      bestLike[b] = likelihood(Q[b], F_org, d.nSites, nPop,d.genos);
	      fprintf(stderr,"like after %f\n",bestLike[b]);
	    }
	    break;
	  }
	}
      } else{   
	if(method==0){
	  // adjusted, unaccelerated EM
	  em(Q[b], F, d.nSites, N, nPop,d.genos,F_new,Q_new[b],F_org);
	} else{
	  // adjusted, accelerated EM
	  if(0==emAccelV2(d, N, nPop, F, Q[b], F_new, Q_new[b],lold,F_org, nit, b, Qconv, Qtol, tol)){
	    if(b<nConv){
	      double tmpLike =  likelihood(Q_new[b], F_new, d.nSites, nPop,d.genos);
	      // stores F with max likelihood, so can be written later
	      if(b==0){
		for(int i=0;i<d.nSites;i++){
		  for(int j=0;j<nPop;j++){
			F_1stRun[i][j] = F_new[i][j];
		  }
		}
	      } else{
		if(tmpLike > bestLike[highestLike]){
		  highestLike = b;
		  for(int i=0;i<d.nSites;i++){
		    for(int j=0;j<nPop;j++){
		      F_1stRun[i][j] = F_new[i][j];
		    }
		  }
		}
	      }
	      // stores all likelihoods so max can be found
	      bestLike[b] = tmpLike;
	      fprintf(stderr,"like after %f\n", bestLike[b]);
	    }
	    break;
	  }
	  
	}
	// swaps adresses of F and F_new, as F_new has the results from the iteration just paseed
	std::swap(F,F_new);
      }
      // swaps adresses of Q and Q_new, as F_new has the results from the iteration just paseed
      std::swap(Q,Q_new);

      //stopping criteria, for EM unaccelerated
      if((nit%10)==0 and method == 0){ 
	double lik = likelihood(Q[b], F, d.nSites, nPop,d.genos);
	if(b==0){
	    fprintf(stderr,"iter[%d] last like is=%f thres=%f\t",nit,likeLast,calcThres(Q[b],Q_new[b],nPop));
	    fprintf(stderr,"iter[%d] like is=%f thres=%f\t",nit,lik,calcThres(Q[b],Q_new[b],nPop));
	    fprintf(stderr,"iter[%d] diff in likelihood is=%f\t",nit,std::abs(lik-likeLast));
	    fprintf(stderr,"iter[%d] ",nit);
	    for(int i=0;i<nPop;i++){	      
	      fprintf(stderr,"Q is=%f, ",Q[b][i]);
	      
	    }
	    fprintf(stderr,"\n");      
	}
	// if convergence based on Q
	if(calcThres(Q[b],Q_new[b],nPop) < Qtol and Qconv>0) {
	  if(b==0){
	    fprintf(stderr,"Convergence achived because diffence in Q values less than %f\n",Qtol);
	  }
	  if(b<nConv){
	    bestLike[b] = lik;
	  }
	  break;
	  // if convergence based on likelihood
	} else if(std::abs(lik-likeLast) < tol and Qconv==0) {
	  if(b==0){
	      fprintf(stderr,"Convergence achived becuase log likelihooditer difference is less than %.8f\n",tol);
	  }
	  if(lik-likeLast<0){
	    if(b==0){
	      fprintf(stderr,"Convergence achived because log likelihood difference was NEGATIVE\n");
	    }
	  }
	  // storing likelihoods from conv runs
	  if(b<nConv){
	    bestLike[b] = lik;
	  }
	  break;
	  }
	likeLast = lik;
      }
      // if no more iterations store likelihood no matter what
      if((nit+1)==maxIter){
	bestLike[b] = likeLast;
      }
    }     
    fprintf(stderr,"CONVERGENCE!\n");
    if(b<nConv){
      fprintf(stderr,"This many iterations %i for run %i\n",nit,b);
      fprintf(flog,"This many iterations %i for run %i\n",nit,b);
      fprintf(stderr,"\n");
      fprintf(flog,"\n");	
    }
    for(int i=0;i<d.nSites;i++){
      for(int k=0;k<nPop;k++){
	// original ref panel freqs
	F[i][k] = F_orgOrg[i][k];
	F_new[i][k] = F_orgOrg[i][k];
      }
    }    
    // to find the Q with lowest likelihood, after 10 first runs
    if(b==(nConv-1)){
      for(int j=nConv;j<(nBoot+nConv);j++){
	for(int k=0;k<nPop;k++){	
	  // make best estimated Q starting guess for bootstrap
	  Q[j][k] = Q[highestLike][k];
	  Q_new[j][k] = Q[highestLike][k];
	}
      }
    }
  }
  
  for(int i=0;i<nConv;i++){
    fprintf(stderr,"best like %f after %i!\n",bestLike[i],i);
    fprintf(flog,"best like %f after %i!\n",bestLike[i],i);
    for(int j=0;j<nPop;j++){
      fprintf(stderr,"Q %f ",Q[i][j]);
      fprintf(flog,"Q %f ",Q[i][j]);
    }
    fprintf(stderr," after %i!\n",i);
    fprintf(flog," after %i!\n",i);
  } 

  fprintf(stderr,"\n");
  fprintf(flog,"\n");
  
  fprintf(stderr,"Estimated  Q = ");
  fprintf(flog,"Estimated  Q = ");
  for(int i=0;i<nPop;i++){
    fprintf(flog,"%f ",Q[highestLike][i]);
    fprintf(stderr,"%f ",Q[highestLike][i]);
  }
  fprintf(flog,"\n");
  
  fprintf(stderr,"best like %f after %i runs!\n",bestLike[highestLike],highestLike);
  fprintf(flog,"best like %f after %i runs!\n",bestLike[highestLike],highestLike);
  
  
  /////////////////////////////////////////////////////////// done - make output and clean /////////////////////////////////  
  // Print F and Q in files
  

  FILE *fp=openFile(outfiles,".qopt");
  // nBoot + 1, because first value is estimated Q
  printDouble(Q,nBoot+nConv,nPop,highestLike,nConv,ref.populations,fp);
  fclose(fp);
  fprintf(stderr,"FIRST row of .qopt file is BEST estimated Q, rest are nBoot bootstrapping Qs\n");
  fprintf(flog,"FIRST row of .qopt file is BEST estimated Q, rest are nBoot bootstrapping Qs\n");

  
  // only if certain flag, print adjusted freqs
  if(printFreq>0){
    gzFile fpGz=openFileGz(outfiles,".fopt.gz");
    printDoubleGz(F_1stRun,ref.refSites,nPop,ref.id,ref.populations,fpGz);
    gzclose(fpGz);
  }
  double nop;
  
  if(doAdjust>0 and method>0){
    emAccelV2(d, NULL, 0, NULL, NULL, F_new, Q_new[0],nop,NULL, 0,1,0,0,0);
  } else if(method > 0){
    emAccelUnadjustedV2(d, NULL, 0, NULL, NULL, Q_new[0],nop,NULL, 0, 1,0,0,0);   
  }

  //deallocate memory 
  dalloc(F,d.nSites); 
  dalloc(F_1stRun,d.nSites);  
  dalloc(F_new,d.nSites);  
  dalloc(Q,nBoot+nConv); 
  dalloc(Q_new,nBoot+nConv);   
  dalloc(F_org,d.nSites);
  dalloc(F_orgOrg,d.nSites);
  
  delete [] N;
  delete [] bestLike;
 
  d.idMap.clear();
  dallocBeagle(d);
  dallocBeagle(dOrg);
  dallocRefPanel(ref,nPop);

  for(std::vector<char*>::iterator it = tmpRef.begin(); it != tmpRef.end(); ++it){
    free(*it);
  } 
  tmpRef.clear();

  for(int i=0;1&&i<dumpedFiles.size();i++){
    free(dumpedFiles[i]);
  }
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  

  // print to log file
  fprintf(flog, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(flog, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fclose(flog); 
  return 0;
    
 }

