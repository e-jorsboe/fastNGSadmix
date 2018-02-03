/*
  log:
  g++ fastNGSadmix.cpp -lz -lpthread  -O3 -o fastNGSadmix

  log: (with readplink function)
  g++ fastNGSadmix.cpp readplinkV3.cpp -lz -lpthread  -O3 -o fastNGSadmix


  debug:
  g++ fastNGSadmix.cpp -lz -lpthread -ggdb -O3 -o fastNGSadmix

  debug: (with readplink function)
  g++ fastNGSadmix.cpp readplinkV3.cpp -lz -lpthread -ggdb -O3 -o fastNGSadmix

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
#include "readplinkV3.h"


//this is the max number of bytes perline
#define LENS 100000
//if we catch signal then quit program nicely
int SIG_COND =1;

double errTolMin=1e-5;
double errTolStart=0.05;
//frequencies and admixture coef cannot be less than this or more than 1-this
double errTol=errTolStart;

// checks if string is number
int validDouble(std::string someString){
  int isNumber = 0;
  int hasPoint = 0;
  for(int i = 0; i<someString.length(); i++){
    char s = someString[i];
    if(s =='.'){
      if(hasPoint){
	isNumber = 0;
	break;
      }
      hasPoint = 1;
      isNumber = 1;
    } else {
      isNumber = (s=='0' or s=='1' or s=='2' or s=='3' or s=='4' or s=='5' or s=='6' or s=='7' or s=='8' or s=='9');
    }

    // leaves loop if any not numbers
    if(isNumber==0){
      break;
    }
  }
  return(isNumber);
  
}



void minus1d(std::vector<double> &fst, std::vector<double> &sec,size_t x, std::vector<double> &res){
  for(size_t i=0;i<x;i++){
      res[i] = fst[i]-sec[i];
    }
}

void minus(std::vector< std::vector<double> >  &fst,std::vector< std::vector<double> > &sec,size_t x,size_t y,std::vector< std::vector<double> > &res){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++){
      res[i][j] = fst[i][j]-sec[i][j];
    }
  }
}

double sumSquare(std::vector< std::vector<double> >  &mat,size_t x,size_t y){
  double tmp=0;
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++){
      tmp += mat[i][j]*mat[i][j];
    }
  }
  return tmp;
}

double sumSquare1d(std::vector<double >  &mat,size_t x){
  double tmp=0;
  for(size_t i=0;i<x;i++){
      tmp += mat[i]*mat[i];
  }
  return tmp;
}

double calcThres(std::vector<double>  &d1, std::vector<double>  &d2, int x){
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
void map2domainQ(std::vector<double> &Q, int nPop){  
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
void map2domainF(std::vector< std::vector<double> > &F, int nSites, int nPop){
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

std::vector<std::string> dumpedFiles;
FILE *openFile(const char* a,const char* b){
  if(0){
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  }
  std::string c1(a, strlen(a));
  std::string c2(b, strlen(b));

  std::string c = c1 + c2;
  fprintf(stderr,"\t-> Dumping file: %s\n",c.c_str());
  if(0&&fexists(c.c_str())){
    fprintf(stderr,"File: %s exists will exist\n",c.c_str());
    fflush(stderr);
    exit(0);
  }

  dumpedFiles.push_back(c);

  FILE *fp = fopen(c.c_str(),"w");
  if (not(fp)){
    fprintf(stderr,"File: %s cannot be created specify valid path\n",c.c_str());
    fflush(stderr);
    exit(0);
  } 
  
  
  return fp;
}

gzFile openFileGz(const char* a,const char* b){
  if(0){
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  }

  std::string c1(a, strlen(a));
  std::string c2(b, strlen(b));

  std::string c = c1 + c2;
  
  fprintf(stderr,"\t-> Dumping file: %s\n",c.c_str());
  if(0&&fexists(c.c_str())){
    fprintf(stderr,"File: %s exists will exist\n",c.c_str());
    fflush(stderr);
    exit(0);
  }
  
  dumpedFiles.push_back(c);
  gzFile fp = gzopen(c.c_str(),"w");

  return fp;
}




// 0 indexed
double getGeno(const std::vector<double> &g, int row, int col){
  if(col>2){
    fprintf(stderr,"only one individaul in beagle file, thus only 3 columns \n");
    exit(0);
  }
  return(g[3*row+col]);
}


double getFreq(const std::vector<double> &freq, int pops, int row, int col){
  // vector is coded so first pops entries are row 0 (0pops...pops-1)
  // then from pops...2*pops-1 is row 1 and so on
  if(col>pops){
    fprintf(stderr,"only %i pops in ref panel \n",pops);
    exit(0);
  }
  return(freq[pops*row+col]);
}


//some struct with all the data from the beagle file
typedef struct{
  std::vector<double> genos;  
  std::vector<char> major;
  std::vector<char> minor;
  // for snp ids, chr_pos
  std::vector<std::string> id;
  int nSites;
  int nInd;
  // map of ids in beagle file for finding overlap with ref
  // id is like this: chr_pos 
  std::map <std::string,int> idMap;
  
}bgl;



// refPanel struct for reading in refPanel with header
typedef struct{
  std::vector<std::string> id;
  std::vector<int> chr;
  std::vector<int> pos;
  std::vector<std::string> name;
  std::vector<char> A0;
  std::vector<char> A1;
  // change this to std::vector<double>, create function that can return freq
  // has to know nSites and nPops
  std::vector<double> freqs;
  int refSites;
  int pops;
  std::vector<std::string> populations;
  // has map of column to keep in ref for calculations
  // is coded so key is old column number, from inputted ref
  // value is new column (in ref for analysis) number + 1 (cause has to be above 0)
  std::map <int,int> colsToKeep;
  std::map <std::string,int> popsToKeep;
  
}refPanel;
 




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


bgl readBeagle(const char* fname, const std::map <std::string,int> &overlap) {
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
  ret.nSites = overlap.size();

  int refIndex = 0;
  int bglIndex = 0;

  std::string dummyID;
  char dummyChar;

  ret.id.assign(ret.nSites,dummyID);
  ret.major.assign(ret.nSites,dummyChar);
  ret.minor.assign(ret.nSites,dummyChar);
  ret.genos.assign(ret.nSites*3,0);
    
  while(gzgets(fp,buf,LENS)!=NULL){

    // puts id of all sites in map for fast lookup
    char* bglID = strtok(buf,delims);
    std::string bglIDstring(bglID, strlen(bglID));
    // because allele might be coded 0,1,2,3    
    char A0 = intToChar(strtok(NULL,delims)[0]);
    char A1 = intToChar(strtok(NULL,delims)[0]);
    
    if( tolower(A0) < tolower(A1) ){
      bglIDstring = bglIDstring + "_" + A0 + "_" + A1;
    } else{
      bglIDstring = bglIDstring + "_" + A1 + "_" + A0;
    }
    // this keeps track of which position in beagle file a site is (+1 to be able to have .count() return TRUE)
    ret.idMap[bglIDstring] = refIndex+1;
    
    //then loop over the vector and parsing every line      
    if(overlap.count(bglIDstring)>0){
      
      ret.id.at(refIndex)=bglIDstring;
      ret.major.at(refIndex)=A0;
      ret.minor.at(refIndex)=A1;
      double tmpS = 0.0;     
      for(int i=0;i<ret.nInd*3;i++){
	double gl = atof(strtok(NULL,delims));
	ret.genos.at(refIndex*3+i)=gl;
	if(gl<0){
	  fprintf(stderr,"Likelihoods must be positive\n");
	  fprintf(stderr,"site %d ind %d geno %d has value %f\n",bglIndex,int(i*1.0/3),i%3,getGeno(ret.genos,refIndex,i));
	  exit(0);
	}
	tmpS+=gl;
	
	if(i==2 and !(tmpS>0)){
	  fprintf(stderr,"The sum of likelihoods for a genotypes must be positive\n");
	  fprintf(stderr,"individual %d site %d has sum %f\n",i,bglIndex,tmpS);
	  exit(0);
	}
      }
      
      // counts which line of overlapping sites between bgl and ref panel
      refIndex++;  
    }
    // counts which line of beagle file
    bglIndex++;
  }
    
  //clean up filepointer
  gzclose(fp); 
  return ret;
}

// read in plink file and converts to a bgl struct (beagle file)
bgl readPlinkToBeagle(const char* plinkName, const std::map <std::string,int> &overlap) {

  plink pl = readplink(plinkName);
  if(pl.fam.individuals > 1){
    fprintf(stderr,"More than one individual in input plink file - should only be one! \n");
    exit(0);
  }
  bgl b;
  b.nSites = overlap.size();
  // can only have one individual is this program
  b.nInd = 1; 
  int beagleIndex = 0;

  char dummyChar;
  std::string dummyID;
  b.id.assign(b.nSites,dummyID);
  b.major.assign(b.nSites,dummyChar);
  b.minor.assign(b.nSites,dummyChar);
  b.genos.assign(b.nSites*3,0);
  
  for(int i=0;i<(pl.y);i++){
    // bim id has chr_pos_A0_A1 ID (A0, A1 ordered alphabetically)
    if(overlap.count(pl.bim.id[i])<1){

      continue;
    }
    if(beagleIndex+1 > overlap.size()){
      fprintf(stderr,"Duplicated sites - multiple sites with same position in plink file!! \n");
      fprintf(stderr,"Use plink2 with --list-duplicate-vars suppress-first ids-only - and then --exclude \n");
      exit(0);
    }
    
    if(pl.d[0][i]==0){
      b.genos.at(beagleIndex*3)=0.0;
      b.genos.at(beagleIndex*3+1)=0.0;
      b.genos.at(beagleIndex*3+2)=1.0;
      
    } else if(pl.d[0][i]==1){
      b.genos.at(beagleIndex*3)=0.0;
      b.genos.at(beagleIndex*3+1)=1.0;
      b.genos.at(beagleIndex*3+2)=0.0;
    
    } else if(pl.d[0][i]==2){
      b.genos.at(beagleIndex*3)=1.0;
      b.genos.at(beagleIndex*3+1)=0.0;
      b.genos.at(beagleIndex*3+2)=0.0;

    }

    b.major.at(beagleIndex)=pl.bim.major[i];
    b.minor.at(beagleIndex)=pl.bim.minor[i];
    b.id.at(beagleIndex)=pl.bim.id[i];
    // stores id of all overlapped sites from plink file,
    // has index of which index in beagle file + 1
    b.idMap[pl.bim.id[i]] = beagleIndex +1;
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

refPanel readRefPanel(const char* fname, bgl b, const std::map <std::string,int> &includedPops, int nPop, const std::map <std::string,int> &overlap) {
  const char *delims = "\t \n";
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  char buf[LENS];
  refPanel ref;

  int totalSites = 0;
  int ncols = 0;
  
  // keeps track of which new column (index = newCol-1) has to be above 0 for lookup in map
  int newCol = 1;  
  //find number of columns
  while(gzgets(fp,buf,LENS)!=NULL){
    if(totalSites==0){      
      char* columnID = strtok(buf,delims);      
      while(columnID!=NULL){
	ncols++;
	// first 6 columns not freqs and has to at most K new columns included in ref
	if(ncols>6){
	  // if in supplied populations or if no populations supplied include
	  if(includedPops.count(columnID)>0 or includedPops.empty()){
	    // prints out which populations chosen
	    fprintf(stderr,"Chosen pop %s\n",columnID);
	    std::string columnIDstring(columnID, strlen(columnID));
	    // for which columns to keep
	    ref.populations.push_back(columnIDstring);
	    // keep track of which new column it will be
	    ref.popsToKeep[columnIDstring] = newCol;
	    // so that can translate from org column (where 7th column is first freq column) to new column (that is +1 here for lookup purposes)
	    ref.colsToKeep[ncols-7] = newCol;
	    newCol++;
	  } 
	}
	columnID = strtok(NULL,delims);	
      }
      
    }
    
    
    totalSites++;   
  }

  gzclose(fp); 
    
  if(ncols<7){
    // has to have at least 7 columns 
    fprintf(stderr,"Too few cols, ncols=%d\n",ncols);
    exit(0);
  }

  ref.pops = ref.popsToKeep.size();
  ref.refSites = b.nSites;

  std::string dummyID;
  char dummyChar;

  ref.freqs.assign(ref.pops*ref.refSites,0);
  ref.id.assign(ref.refSites,dummyID);
  ref.chr.assign(ref.refSites,0);
  ref.pos.assign(ref.refSites,0);
  
  ref.name.assign(ref.refSites,dummyID);
  // ref has A,C,G,T alleles
  ref.A0.assign(ref.refSites,dummyChar);
  ref.A1.assign(ref.refSites,dummyChar);
   
  gzFile fp1 = NULL;
  fp1=gzopen(fname,"r");
  
  // for keeping track of which index in new ref with
  int refIndex = 0;
  int refSite = 0;
  
  while(gzgets(fp1,buf,LENS)!=NULL){

    // looking at id value chr_pos for detecting overlap
    char* id = strtok(buf,delims);
    std::string stringID(id,strlen(id));
    
    int refChr = atoi(strtok(NULL,delims));
    int refPos = atoi(strtok(NULL,delims));
    
    char* name = strtok(NULL,delims);
    std::string stringName(name,strlen(name));

    // ref has A,C,G,T alleles
    char A0 = intToChar(strtok(NULL,delims)[0]);
    char A1 = intToChar(strtok(NULL,delims)[0]);

    // sorting also with alleles
    if( tolower(A0) < tolower(A1) ){
      stringID = stringID + "_" + A0 + "_" +A1;
    } else{
      stringID = stringID + "_" + A1 + "_" + A0;
    }
    
    // check if site is in overlap with beagle file
    // otherwise continues to next site in ref
    if(overlap.count(stringID) > 0){

      // index of site in overlap - and thereby also in beagle file (only overlapped sites)
      int overlapIndex = b.idMap[stringID]-1;

      // the idMap of the beagle struct has index sites was placed on in beagle file (index+1)
      ref.id.at(overlapIndex)=stringID;

      ref.chr.at(overlapIndex)=refChr;
      ref.pos.at(overlapIndex)=refPos;

      ref.name.at(overlapIndex)=stringName;
      // ref has A,C,G,T alleles
      ref.A0.at(overlapIndex)=A0;
      ref.A1.at(overlapIndex)=A1;
      
      //    reading in ref freqs
      for(int i=0;i<(ncols-6);i++){
	// check if org column to keep and thereby pop to keep in ref
      if(ref.colsToKeep.count(i)>0){
	// if bgl 1_1 A B GL(AA) GL(AB) GL(BB) Then ref 1 1 rs1 B A 1-f(B)
	// minor is last allele in beagle file
	if(ref.A0[overlapIndex]==b.minor[overlapIndex]){
	  // new col has to be - 1 for right index
	  //	  ref.freqs[refIndex][ref.colsToKeep[i]-1] = 1 - atof(strtok(NULL,delims));
	  ref.freqs.at(ref.pops*overlapIndex+(ref.colsToKeep[i]-1)) = 1 - atof(strtok(NULL,delims));
	} else{
	  //	  ref.freqs[refIndex][ref.colsToKeep[i]-1] = atof(strtok(NULL,delims));
	  ref.freqs.at(ref.pops*overlapIndex+(ref.colsToKeep[i]-1)) = atof(strtok(NULL,delims));
	}
	if(getFreq(ref.freqs,ref.pops,overlapIndex,ref.colsToKeep[i]-1)<0){
	  fprintf(stderr,"Frequencies must be positive\n");
	  fprintf(stderr,"site %d, pop %d, has value %f\n",refSite,i,getFreq(ref.freqs,ref.pops,overlapIndex,ref.colsToKeep[i]-1));
	  exit(0);
	}
      } else{
	// it has to skip cols that are not to be read in and move to next column which will be checked
	strtok(NULL,delims);
      }
      
      }
      refIndex++; 
    }
    refSite++;
    // only here if site was included in new ref
    
  }

  if(refIndex!=overlap.size()){
    fprintf(stderr,"refIndex %i, overlap sites %lu\n",refIndex,overlap.size());
    fprintf(stderr,"Seems like there are duplicate ids in input or reference panel!\n");
    exit(0);
  }

  gzclose(fp1); 
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
void readDouble1d(std::vector <double> &d,int nPop,const char*fname, std::map<std::string,int> popsToKeep){
  fprintf(stderr,"Opening nInd file: %s with nPop=%d\n",fname,nPop);
  const char*delims=" \n\t";
  gzFile fp = NULL;
  if((fp=gzopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];

  d.assign(nPop,0);
  
  std::vector<char*> tmp;
  // keeps track of org index of pop to keep and new index of pop to keep
  std::map<int,int> toKeep;

  gzgets(fp,buf,LENS);
    
  // looking at id value chr_pos for detecting overlap
  char* word = strtok(buf,delims);
  std::string stringID(word,strlen(word));
  int orgCol = 0;
  
  while(word!=NULL){
    std::string stringWord(word,strlen(word));
    if(popsToKeep.count(stringWord)>0){	
      // creates map of <nInd index, ref index> so can map from nInd order to ref order
      // so if first element in nInd is second in ref
      // the d array with nInd values has second element equal to nInd first value
      toKeep[orgCol] = popsToKeep[stringWord];	
    }
    word = strtok(NULL,delims); 
    orgCol++;
  }
     
  // first goes through the first line with names of pops, finds which should be included  
  // keeps track of which value at in nInd file, newCol has to be index+1, because has to be > 0 for lookup

  if(toKeep.size()!=popsToKeep.size()){
    fprintf(stderr,"nInd and ref panel do not have same size!\n");
    exit(0);
  }
  
  int index = 0;
  gzgets(fp,buf,LENS);
  word = strtok(buf,delims);        
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
  
  gzclose(fp);
}

void printDouble(const std::vector< std::vector<double> > &ret,size_t x,size_t y, int highestLike, int nConv, const std::vector<std::string> &populations ,FILE *fp){
  for(size_t i=0;i<x;i++){
    if(i==0){
      for(size_t j=0;j<y;j++){
	fprintf(fp,"%s ",populations[j].c_str());
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


void printDoubleGz(const std::vector< std::vector<double> > &ret, size_t x, size_t y, const std::vector<std::string> &id, const std::vector<std::string> &populations ,gzFile fp){

 for(size_t i=0;i<x;i++){
    if(i==0){
      gzprintf(fp,"id ");
      for(size_t j=0;j<y;j++){
	gzprintf(fp,"%s ",populations[j].c_str());
      }
      gzprintf(fp,"\n");
    }
    gzprintf(fp,"%s ",id[i].c_str());
    for(size_t j=0;j<y;j++){
      gzprintf(fp,"%.4f ",1-ret[i][j]);
    }
    gzprintf(fp,"\n");
  }
}

// calculate log(likelihood) from likelihood function
double likelihood(const std::vector<double> &Q, const std::vector< std::vector<double> > &F,int nSites, int nPop, const std::vector<double> &genos, int ploidy){
  double prod_ind = 0.0;
  for(int j = 0; j < nSites; j++) {
    double freq = 0.0;
    for(int k = 0; k < nPop; k++) {
      freq += (F[j][k])*Q[k];
    }
    // has to be like this, as I sort freqs
    // prior to running this program
    double f = freq;

    // ploidy will be either 1 or 2!
    if(ploidy==1){
      double sum = getGeno(genos,j,0)*f;
      sum +=  getGeno(genos,j,1)*(1-f);
      prod_ind += log(sum);
      
    } else if(ploidy==2){
      double sum = getGeno(genos,j,0)* f * f;
      sum +=  getGeno(genos,j,1)*2*f*(1-f);
      sum +=  getGeno(genos,j,2)*(1-f)*(1-f);
      prod_ind += log(sum);
    } 
  }
  return prod_ind;
}

// does bootstrapping sampling nSites random sites with replacement
void bootstrap(const std::vector<double> &genosOrg, std::vector<double> &genos, const std::vector< std::vector<double> > &F_orgOrg, std::vector< std::vector<double> > &F_org, std::vector< std::vector<double> > &F, int nPop, int nSites, int ploidy) {
  for(int j=0;j<nSites;j++){
    // generate random int from 0 to (nSites-1)
    int row = std::rand() % nSites;
    if(ploidy==1){
      genos[3*j+0] = getGeno(genosOrg,row,0);
      genos[3*j+1] = getGeno(genosOrg,row,1);

    } else if(ploidy==2){        
      genos[3*j+0] = getGeno(genosOrg,row,0);
      genos[3*j+1] = getGeno(genosOrg,row,1);
      genos[3*j+2] = getGeno(genosOrg,row,2);
    }
    for(int k = 0; k < nPop; k++) {
      F[j][k] = F_orgOrg[row][k];
      F_org[j][k] = F_orgOrg[row][k];
    }
  }
}

// em algorithm not adjusting F at every step
void emUnadjusted(std::vector<double> &Q, std::vector< std::vector<double> > &F, int nSites, int nPop, const std::vector<double> &genos, std::vector<double> &Q_1, int ploidy) {
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

      if(ploidy==1){
	double pp0=(1-fpart)*getGeno(genos,j,1);
	double pp1=fpart*getGeno(genos,j,0);
	double sum=pp0+pp1;
	expGG =(pp1)/sum;
      } else if(ploidy==2){
	  
	// pre GL (sites x 3) * (adjusted freq)
	// for calculating H range 0-2, this is the expected genotype
	double pp0=(1-fpart)*(1-fpart)*getGeno(genos,j,2);
	double pp1=2*(1-fpart)*fpart*  getGeno(genos,j,1);
	double pp2=fpart*fpart*        getGeno(genos,j,0);
	double sum=pp0+pp1+pp2;
	expGG =(pp1+2*pp2)/sum;
      }    	  
      
    }
    for(int k=0;k<nPop;k++){
      sumAG[k] += expGG/(fpart) * (Q[k] * F[j][k]); 
      sumBG[k] += (ploidy-expGG)/fpartInv * (Q[k] * (1-F[j][k]));
    }
    
  }
  for(int k=0;k<nPop;k++){        
    Q_1[k]=(sumAG[k] + sumBG[k])/(ploidy*nSites*1.0);        
  }
  
  map2domainQ(Q_1,nPop);
}

// em algorithm adjusting F at every step
void em(std::vector<double> &Q, std::vector< std::vector<double> > &F, int nSites, const std::vector<double> &nInd, int nPop, const std::vector<double> &genos, std::vector< std::vector<double> > &F_1, std::vector<double> &Q_1, std::vector< std::vector<double> > &F_org, int ploidy) {
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
    double sum=0;
    for(int k=0;k<nPop;k++){ 
      // admixture adjusted freq, for each pop
      fpart += F[j][k] * Q[k];
      fpartInv += (1-F[j][k]) * Q[k];
          
      if(ploidy==1){
	double pp0=(1-fpart)*getGeno(genos,j,1);
	double pp1=fpart*getGeno(genos,j,0);
	sum=pp0+pp1;
	expGG = (pp1)/sum;
      } else if(ploidy==2){	
	// pre GL (sites x 3) * (adjusted freq)
	// for calculating H range 0-2, this is the expected genotype
	double pp0=(1-fpart)*(1-fpart)*getGeno(genos,j,2);
	double pp1=2*(1-fpart)*fpart*  getGeno(genos,j,1);
	double pp2=fpart*fpart*        getGeno(genos,j,0);
	sum=pp0+pp1+pp2;
	expGG = (pp1+2*pp2)/sum;
      }
    }

    for(int k=0;k<nPop;k++){
      // similar to (H/(q*f))*q, for jth marker      
      sumAG[k] = (expGG) / (fpart) * (Q[k]*F[j][k]);
      sumBG[k] = (ploidy-expGG) / fpartInv * (Q[k]*(1-F[j][k]));
      sumA[k] += sumAG[k];
      sumB[k] += sumBG[k];
      sumAG[k] += nInd[k]*ploidy*F_org[j][k];
      sumBG[k] += ploidy*nInd[k]-(ploidy*nInd[k]*F_org[j][k]);
           
    }
    
    for(int k=0;k<nPop;k++){
      // adjust with ref panel, so we have input + ref expected number of alleles
      F_1[j][k]=sumAG[k]/(sumAG[k]+sumBG[k]);
      
    }      
  }

  for(int k=0;k<nPop;k++){  
    Q_1[k]=(sumA[k] + sumB[k])/(ploidy*nSites*1.0);
    
  }
  
  map2domainQ(Q_1,nPop);
  map2domainF(F_1,nSites,nPop);
}


int emAccelUnadjustedV2(const std::vector<double> &genos, const std::vector<double> &nInd, int nPop, std::vector< std::vector<double> > &F, std::vector<double> &Q, std::vector<double> &Q_new, int nit, int boot, int Qconv, double Qtol, double tol, int nSites, int ploidy){
 
  double stepMin = 1;
  double stepMax0 = 1;
  static double stepMax = stepMax0;
  double mstep = 4;
  double objfnInc = 1;
  //we make these huge struc
  static std::vector<double> Q_em1;
  static std::vector<double> Q_diff1;
  static std::vector<double> Q_em2;
  static std::vector<double> Q_diff2;
  static std::vector<double> Q_diff3;
  static std::vector<double> Q_tmp;
  static std::vector<double> Q_tmpDiff;
  
  if(Q_em1.empty()){
    Q_em1.assign(nPop,0);
    Q_diff1.assign(nPop,0);
    Q_em2.assign(nPop,0);
    Q_diff2.assign(nPop,0);
    Q_diff3.assign(nPop,0);
    Q_tmp.assign(nPop,0);
    Q_tmpDiff.assign(nPop,0);

  }
  // first EM run
  emUnadjusted(Q, F, nSites, nPop, genos, Q_em1, ploidy);

  minus1d(Q_em1,Q,nPop,Q_diff1);
  double sr2 = sumSquare1d(Q_diff1,nPop);
  // checks if convergence
  if(sqrt(sr2)<tol or (calcThres(Q,Q_em1,nPop) < Qtol and Qconv>0)){
    //fprintf(stderr,"like is %f\n",likelihood(Q_new, F, nSites, nPop,genos));    
    return 0;  
  }
  // second EM run
  emUnadjusted(Q_em1, F, nSites, nPop, genos, Q_em2, ploidy);
  minus1d(Q_em2,Q_em1,nPop,Q_diff2);
  double sq2 = sumSquare1d(Q_diff2,nPop);
  // checks if convergence - a second time
  if(sqrt(sq2)<tol  or (calcThres(Q_em1,Q_em2,nPop) < Qtol and Qconv>0)){
    //fprintf(stderr,"like is %f\n",likelihood(Q_new, F, nSites, nPop,genos));
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
    emUnadjusted(Q_new, F, nSites, nPop,genos,Q_tmp, ploidy);    
    minus1d(Q_tmp,Q_new,nPop,Q_tmpDiff);
    // adopted from squarem2 from SQUAREM package
    double res = sumSquare1d(Q_tmpDiff,nPop);
    double parnorm = (1/std::sqrt(nPop))*sumSquare1d(Q_tmpDiff,nPop);
    double kres = 1 + parnorm + sq2;
    if(res <= kres){
      Q_new.swap(Q_tmp);
    } else{
      Q_new.swap(Q_em2);
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
      double lnew = likelihood(Q_new, F, nSites, nPop,genos, ploidy);
      if(lnew!=lnew){
	fprintf(stderr,"likelihood is nan, probably because dividing by 0, go fix ref panel or input!\n");
	exit(0);
      }      
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
int emAccelV3(const std::vector<double> &genos, const std::vector<double> &nInd, int nPop, std::vector< std::vector<double> > &F, std::vector<double> &Q, std::vector< std::vector<double> > &F_new, std::vector<double> &Q_new, std::vector< std::vector<double> > &F_org, int nit, int boot, int Qconv, double Qtol, double tol, int nSites, int ploidy){

  double stepMin = 1;
  double stepMax0 = 1;
  static double stepMax = stepMax0;
  double mstep = 4;
  double objfnInc = 1;

  //we make these huge structures static such that we just allocate them the first time
  static std::vector<double> Q_em1;
  static std::vector<double> Q_diff1;
  static std::vector<double> Q_em2;
  static std::vector<double> Q_diff2;
  static std::vector<double> Q_diff3;
  static std::vector<double> Q_tmp;
  static std::vector<double> Q_tmpDiff;
  
  static std::vector< std::vector<double> > F_em1(nSites, std::vector<double>(nPop));
  static std::vector< std::vector<double> > F_diff1(nSites, std::vector<double>(nPop));
  static std::vector< std::vector<double> > F_em2(nSites, std::vector<double>(nPop));
  static std::vector< std::vector<double> > F_diff2(nSites, std::vector<double>(nPop));
  static std::vector< std::vector<double> > F_diff3(nSites, std::vector<double>(nPop));
  static std::vector< std::vector<double> > F_tmp(nSites, std::vector<double>(nPop));
  static std::vector< std::vector<double> > F_tmpDiff(nSites, std::vector<double>(nPop));

  if(Q_em1.empty()){  
    Q_em1.assign(nPop,0);
    Q_diff1.assign(nPop,0);
    Q_em2.assign(nPop,0);
    Q_diff2.assign(nPop,0);    
    Q_diff3.assign(nPop,0);      
    Q_tmp.assign(nPop,0);        
    Q_tmpDiff.assign(nPop,0);
  }
 
  // first EM run
  em(Q, F, nSites, nInd, nPop,genos, F_em1, Q_em1, F_org, ploidy);
  minus(F_em1,F,nSites,nPop,F_diff1);
  minus1d(Q_em1,Q,nPop,Q_diff1);
  double sr2 = sumSquare1d(Q_diff1,nPop) + sumSquare(F_diff1,nSites,nPop);
  // checks if convergence
  if(sqrt(sr2)<tol or (calcThres(Q,Q_em1,nPop) < Qtol and Qconv>0)){
    //fprintf(stderr,"like is %f\n",likelihood(Q_new, F_new, nSites, nPop,genos));    
    return 0;
  }
  // second EM run
  em(Q_em1, F_em1, nSites, nInd, nPop,genos, F_em2, Q_em2, F_org, ploidy);
  minus(F_em2,F_em1,nSites,nPop,F_diff2);
  minus1d(Q_em2,Q_em1,nPop,Q_diff2);
  double sq2 = sumSquare1d(Q_diff2,nPop) + sumSquare(F_diff2,nSites,nPop);
  // checks if convergence - a second time
  if(sqrt(sq2)<tol or (calcThres(Q_em1,Q_em2,nPop) < Qtol and Qconv>0)){
    //fprintf(stderr,"like is %f\n",likelihood(Q_new, F_new, nSites, nPop,genos));
    return 0;
  }
  minus(F_diff2,F_diff1,nSites,nPop,F_diff3);
  minus1d(Q_diff2,Q_diff1,nPop,Q_diff3);
  double sv2 = sumSquare1d(Q_diff3,nPop) + sumSquare(F_diff3,nSites,nPop);
  double alpha = sqrt(sr2/sv2);  
  // makes sure alpha does not go below 1 and above stepMax
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<nSites;i++){
    for(size_t j=0;j<nPop;j++){
      // based on ngsAdmix approach
      F_new[i][j] = F[i][j]+2*alpha*F_diff1[i][j]+alpha*alpha*F_diff3[i][j];
      F_tmp[i][j] = 1.0;
    }
  }
  map2domainF(F_new,nSites,nPop);
  for(size_t i=0;i<nPop;i++){
    Q_new[i] = Q[i]+2*alpha*Q_diff1[i]+alpha*alpha*Q_diff3[i];
    Q_tmp[i] = 1.0;
  }
  map2domainQ(Q_new,nPop);
  // if alpha not too close (0.01 close) to 1 
  if (fabs(alpha - 1) > 0.01){
    // we estimate new Q and F, with our inferred Q and F via alpha
    em(Q_new, F_new, nSites, nInd, nPop,genos,F_tmp,Q_tmp,F_org,ploidy);
    minus(F_tmp,F_new,nSites,nPop,F_tmpDiff);
    minus1d(Q_tmp,Q_new,nPop,Q_tmpDiff);
    double res = sumSquare1d(Q_tmpDiff,nPop) + sumSquare(F_tmpDiff,nSites,nPop);
    double parnorm = (1/std::sqrt(nPop))*sumSquare1d(Q_tmpDiff,nPop) + (1/std::sqrt(nSites*nPop))*sumSquare(F_tmpDiff,nSites,nPop);
    double kres = 1 + parnorm + sq2;
    if(res <= kres){
      Q_new.swap(Q_tmp);
      F_new.swap(F_tmp);
    } else{
      Q_new.swap(Q_em2);
      F_new.swap(F_em2);
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
      double lnew = likelihood(Q_new, F_new, nSites, nPop,genos,ploidy);
      if(lnew!=lnew){
	fprintf(stderr,"likelihood is nan, probably because dividing by 0, go fix ref panel or input!\n");
	exit(0);
      }
      fprintf(stderr,"iter[%d] like=%f alpha=%f ",nit,lnew,alpha);
      for(int i=0;i<nPop;i++){	      
	fprintf(stderr,"Q=%f, ",Q_new[i]);
      }
      fprintf(stderr,"\n");
    }
  }
  
  return 1;
}


// finds overlap between input and ref panel based on id chr_pos_A0_A1 (A0,A1 two alleles sorted), assumes no duplicate sites in terms of position and alleles
std::map <std::string,int> findOverlapV3(const char* lname, const char* plinkName, const char* fname, FILE* flog, const std::map <std::string,int> &includedPops, char* pops, double maf){
  std::map <std::string,int> inputSites;
  static std::map <std::string,int> overlap;
  const char *delims = "\t \n";
  // if beagle file input
  if(plinkName==NULL){
    int beagleIndex=0;
    gzFile fp1 = NULL;
    if(Z_NULL==(fp1=gzopen(lname,"r"))){
      fprintf(stderr,"Error opening file: %s\n",lname);
      exit(0);
    }
    char buf1[LENS];
    while(NULL!=gzgets(fp1,buf1,LENS)){
      if(beagleIndex>0){
	char* bglID = strtok(buf1,delims);
	std::string bglString(bglID,strlen(bglID));
	char A0 = intToChar(strtok(NULL,delims)[0]);
	char A1 = intToChar(strtok(NULL,delims)[0]);
	if( tolower(A0) < tolower(A1)){
	  bglString=bglString + "_" + A0 + "_" + A1;
	} else{
	  bglString=bglString + "_" + A1 + "_" + A0;
	}
	if(inputSites.count(bglString)>0){
	  fprintf(stderr,"Duplicate sites in beagle file: %s - Go fix!\n",bglString.c_str());
	  exit(0);
	} else{
	  inputSites[bglString]=1;
	}
      }
      beagleIndex++;
    }    
    gzclose(fp1);
    // if plink file input
  } else{
    plink pl_tmp = readplink(plinkName);
    // reads all plink sites into map, checks for duplicates!
    for(int s=0;s<(pl_tmp.y);s++){
      if(inputSites.count(pl_tmp.bim.id[s])>0){
	fprintf(stderr,"Duplicate sites in plink file: %s - Go fix!\n",pl_tmp.bim.id[s].c_str());
	exit(0);
      }
      else if(pl_tmp.d[0][s]==3){
	continue;
      }
      inputSites[pl_tmp.bim.id[s]]=1;  
    }
    kill_plink(pl_tmp);
  }
  fprintf(stderr,"Input has this many sites without missing data %zu\n",inputSites.size());
  fprintf(flog,"Input has this many sites without missing data %zu\n",inputSites.size());
  // reads ref panel
  gzFile fp2 = NULL;
  char buf2[LENS];
  if(Z_NULL==(fp2=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }

  int refIndex = 0;
  int invalidSites = 0;
  std::map <int,int> colsToRead;
  while(NULL!=gzgets(fp2,buf2,LENS)){
    // looking at id value chr_pos for detecting overlap
    // reading in which columns are the selected pops
    if(refIndex==0){
      int whichCol=0;
      char* refPop = strtok(buf2,delims);
      while(refPop!=NULL){
	std::string refPopString(refPop,strlen(refPop));
	int allPops = toupper(pops[0])=='A' and toupper(pops[1])=='L' and toupper(pops[2])=='L' and pops[3]=='\0';
	if(includedPops.count(refPopString)>0 or allPops){	  
	  colsToRead[whichCol] = 1;
	}
	whichCol++;
	refPop = strtok(NULL,delims);
      }
    
    } else{      
      
      char* id = strtok(buf2,delims);
      std::string refStringID(id,strlen(id));

      // in order to get past the pre-allele stuff
      strtok(NULL,delims);
      strtok(NULL,delims);
      strtok(NULL,delims);

      // in order to get alleles   
      char A0 = intToChar(strtok(NULL,delims)[0]);
      char A1 = intToChar(strtok(NULL,delims)[0]);
      
      if( tolower(A0) < tolower(A1)){
	refStringID=refStringID + "_" + A0 + "_" + A1;
      } else{
	refStringID=refStringID + "_" + A1 + "_" + A0;
      }

      int colBeingRead = 5;
      int skipLine = 0;
    
      while(id!=NULL){	
	std::string refString(id,strlen(id));
	// check if a freq column
	if(colBeingRead>=6){
	  // checks if one of selected columns
	  if(colsToRead.count(colBeingRead)>0){
	    // checks if not double NA for instance
	    if(not validDouble(refString)){
	      skipLine = 1;
	      // checks if freq below maf threshold
	    } else  if(validDouble(refString) and (atof(id)<maf or atof(id)>1-maf)){
	      skipLine = 1;
	      
	    }
	    if(skipLine){
	      invalidSites++;
	    }
	  }
	}
	id = strtok(NULL,delims);
	colBeingRead++;	
      }
      
      if(overlap.count(refStringID)>0){
	fprintf(stderr,"Duplicate site in ref panel: %s - Go fix!\n",refStringID.c_str());     
	exit(0);
	// check if site is in beagle or plink file
	// otherwise continues to next site in ref
      } else {

	if(inputSites.count(refStringID) > 0 and not skipLine){
	  overlap[refStringID] = 1;	
	}
      }
      
    }
    refIndex++;

  }
  
  // because ref index also counts header
  fprintf(stderr,"Ref has this many sites %i\n",refIndex-1);
  fprintf(flog,"Ref has this many sites %i\n",refIndex-1);

  fprintf(stderr,"This many sites in ref are either not-valid-number or below maf in any of the chosen pops %i\n",invalidSites);
  fprintf(flog,"This many sites in ref are either not-valid-number or below maf in any of the chosen pops %i\n",invalidSites);
 
  // starts at 1 to avoid header   
  if(overlap.size()==0){
    fprintf(stderr,"No overlapping sites where found!!\n");
    exit(0);
  }
  gzclose(fp2);
  // cleaning
 
  return(overlap);
}
void info(){
  
  fprintf(stderr,"Arguments:\n");
  fprintf(stderr,"\t-likes Beagle likelihood filename\n");
  fprintf(stderr,"\t-plink Plink file in the binary bed format\n");
  fprintf(stderr,"\t-Nname Number of individuals in each reference populations\n");
  fprintf(stderr,"\t-fname Ancestral population frequencies\n");
  fprintf(stderr,"\t-whichPops Which populations from the ref panel to include in analysis, denotes number of populations (nPop) for admixture estimation\n \t If 'all' all pops in ref are analyzed, must be comma seperated (pop1,pop2,..)\n");
 
  fprintf(stderr,"Optional:\n");
  fprintf(stderr,"\t-haploid Raise this flag if haploid organism being analyzed, first two cols of beagle file will be used - does not work for plink! Only write '-haploid'.\n");
  fprintf(stderr,"\t-out Prefix for output files\n"); 
  fprintf(stderr,"\t-printFreq print admixture adjusted allele frequencies of reference panel + input individual (1: yes, 0: no (default))\n"); 

  fprintf(stderr,"Setup:\n");
  fprintf(stderr,"\t-doAdjust Adjusts the frequencies in the reference populations with the input (1: yes (default), 0: no)\n");
  fprintf(stderr,"\t-seed Seed for initial guess in EM and for bootstrap\n"); 
  fprintf(stderr,"\t-method If 0 no acceleration of EM algorithm (1: yes (default), 0: no)\n");
  fprintf(stderr,"\t-maf Filters away sites with lower maf in any of analyzed pops, default 0\n");

  fprintf(stderr,"Stop chriteria:\n"); 
  fprintf(stderr,"\t-Qconv Stopping criteria based on change in Q (works best when using doAdjust) (1: yes, 0: no (default))\n"); 
  fprintf(stderr,"\t-Qtol Tolerance value for stopping criteria based on change in Q (0.001 (default))\n"); 
  fprintf(stderr,"\t-tol Tolerance for convergence - can only be set for the unaccelerated EM algorithm (EM: 1e-5, EMAcc: 1e-7)\n"); 
  fprintf(stderr,"\t-maxiter Maximum number of EM iterations\n"); 
  fprintf(stderr,"\t-boot Number of bootstrapping iterations, default 0, can at most be 10000, .qopt FIRST row BEST estimated Q, rest bootstraps!!\n"); 
  fprintf(stderr,"\t-conv Number of convergence iterations, each with random starting point, to check if has converged, default 1, can at most be 10\n");
  fprintf(stderr,"\t-randomBoot if 1 takes random Q starting points for each bootstrap, instead of converged upon estimate, default 0\n");

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
  char* pops = '\0';
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
  double maf = 0.00;
  int randomBoot = 0;
  // ploidy only works for 1 and 2 - only beagle files
  int ploidy = 2;
  
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
    // which populations in ref panel to be ananlyzed, must agree with nPop
    else if(strcmp(*argv,"-whichPops")==0 || strcmp(*argv,"-pops")==0) pops=*++argv;
    else if(strcmp(*argv,"-haploid")==0 || strcmp(*argv,"-h")==0 ) ploidy=1;
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
    else if(strcmp(*argv,"-maf")==0) maf=atof(*++argv);
    else if(strcmp(*argv,"-randomBoot")==0) randomBoot=atoi(*++argv); 
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

  if(pops=='\0'){
    fprintf(stderr,"Please supply which populations to be analyzed - 'all' for all pops: -Nname\n");
    info();
 
  }
  // max 10000 bootstraps 10 conv runs and 
  nBoot = std::min(std::max(nBoot,0),10000);
  nConv = std::min(std::max(nConv,1),10);
  
  //out put files
  FILE *flog=openFile(outfiles,".log");

  fprintf(stderr,"Input: likes=%s plink=%s Nname=%s fname=%s outfiles=%s\n",lname,plinkName,Nname,fname,outfiles);
  fprintf(stderr,"Setup: seed=%d method=%d\n",seed,method);
  fprintf(stderr,"Ploidy of %i has been chosen\n\n",ploidy);
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
  fprintf(flog,"Ploidy of %i has been chosen\n\n",ploidy);
  
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
  if(not (toupper(pops[0])=='A' and toupper(pops[1])=='L' and toupper(pops[2])=='L' and pops[3]=='\0')){
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
    overlap = findOverlapV3(lname, NULL, fname,flog,includedPops,pops,maf);
    d=readBeagle(lname,overlap);
    dOrg=readBeagle(lname,overlap);
  } else if(plinkName!=NULL){
    // finds overlapping sites, then reads plink file
    if(ploidy==1){
      fprintf(stderr,"Haploid analysis does not work with plink files!\n");
      fprintf(flog,"Haploid analysis does not work with plink files!\n");
      info();
    }
    overlap = findOverlapV3(NULL, plinkName, fname,flog,includedPops,pops,maf);
    d=readPlinkToBeagle(plinkName, overlap);
    dOrg=readPlinkToBeagle(plinkName, overlap);  
  }

  if(maf>0){
    fprintf(stderr,"Overlap: of %zu sites between input and ref, after maf filter of %f\n",overlap.size(),maf);
    fprintf(flog,"Overlap: of %zu sites between input and ref, after maf filter of %f\n",overlap.size(),maf);  
  } else{
    fprintf(stderr,"Overlap: of %zu sites between input and ref\n",overlap.size());
    fprintf(flog,"Overlap: of %zu sites between input and ref\n",overlap.size());  
  }
  
  refPanel ref;  
  // if nPop == 0 then reads in all of them refs
  ref = readRefPanel(fname,d,includedPops,nPop,overlap);
  if(toupper(pops[0])=='A' and toupper(pops[1])=='L' and toupper(pops[2])=='L' and pops[3]=='\0'){
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

  std::vector< std::vector<double> > F(d.nSites, std::vector<double>(nPop));
  std::vector< std::vector<double> > F_new(d.nSites, std::vector<double>(nPop));
  // F_org is for storing initial freqs from ref panel
  std::vector< std::vector<double> > F_org(d.nSites, std::vector<double>(nPop));
  // for sampling from when doing bootstrap
  std::vector< std::vector<double> > F_orgOrg(d.nSites, std::vector<double>(nPop));
  // for when having to write adjusted freqs
  std::vector< std::vector<double> > F_1stRun(d.nSites, std::vector<double>(nPop));
  // gets freq values from ref struct, where read into
  for(int i=0;i<d.nSites;i++){
    for(int k=0;k<nPop;k++){
      double f = getFreq(ref.freqs,ref.pops,i,k);
      F[i][k] = f;
      F_org[i][k] = f;
      F_orgOrg[i][k] = f;
    }   
  }

  // because has to have conv values for converge runs, and then nBoot bootstrapped values
  std::vector< std::vector<double> > Q(nBoot+nConv, std::vector<double>(nPop));
  std::vector< std::vector<double> > Q_new(nBoot+nConv, std::vector<double>(nPop));
  std::vector<double> sum;
  sum.assign(nBoot+nConv,0);
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
  
  std::vector<double> N;
  N.assign(nPop,0);
  // reading nInd, where colsToKeep from ref to read in the same columns as in ref
  readDouble1d(N,nPop,Nname,ref.popsToKeep);
  fprintf(flog,"Opening nInd file: %s with nPop=%d\n",fname,nPop); 
  for(int i=0;i<nPop;i++){
    fprintf(flog,"Chosen pop %s\n",ref.populations[i].c_str());
    fprintf(flog,"N = %f\n",N[i]);
    fprintf(stderr,"N = %f\n",N[i]);
    // being printed in function to stderr
  }
  fprintf(stderr,"\n");
  fprintf(flog,"\n");

  std::vector<double> bestLike; bestLike.assign(nConv,0);
  int highestLike = 0;
  // initial likelihood   
  double lold = likelihood(Q[0], F_org, d.nSites, nPop, d.genos, ploidy);
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
      bootstrap(dOrg.genos,d.genos,F_orgOrg,F_org,F,nPop,d.nSites,ploidy); 
      fprintf(stderr,"At this bootstrapping: %i out of: %i\n",b-(nConv-1),nBoot);
      fprintf(flog,"At this bootstrapping: %i out of: %i\n",b-(nConv-1),nBoot);
    }
    for(nit=1;SIG_COND and nit<maxIter;nit++) {	
      if(doAdjust==0){
	if(method==0){
	  // unadjusted, unaccelerated EM
	  emUnadjusted(Q[b], F, d.nSites, nPop,d.genos,Q_new[b], ploidy);
	} else{
	  // unadjusted, accelerated EM	  
	  if(0==emAccelUnadjustedV2(d.genos, N, nPop, F, Q[b], Q_new[b], nit, b, Qconv, Qtol, tol, d.nSites, ploidy)){
	    if(b<nConv){
	      // stores all likelihoods so max can be found
	      bestLike[b] = likelihood(Q[b], F_org, d.nSites, nPop,d.genos, ploidy);
	      fprintf(stderr,"like after %f\n",bestLike[b]);
	    }
	    break;
	  }
	}
      } else{   
	if(method==0){
	  // adjusted, unaccelerated EM
	  em(Q[b], F, d.nSites, N, nPop,d.genos,F_new,Q_new[b],F_org, ploidy);
	} else{
	  // adjusted, accelerated EM
	  if(0==emAccelV3(d.genos, N, nPop, F, Q[b], F_new, Q_new[b],F_org, nit, b, Qconv, Qtol, tol, d.nSites, ploidy)){
	    if(b<nConv){
	      double tmpLike =  likelihood(Q_new[b], F_new, d.nSites, nPop,d.genos,ploidy);
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
	F.swap(F_new);
      }
      // swaps adresses of Q and Q_new, as F_new has the results from the iteration just paseed
      Q.swap(Q_new);
      
      //stopping criteria, for EM unaccelerated
      if((nit%10)==0 and method == 0){
	double lik = likelihood(Q[b], F, d.nSites, nPop,d.genos,ploidy);	
	if(likeLast!=likeLast and lik!=lik){
	  fprintf(stderr,"likelihood is nan, probably because dividing by 0, go fix ref panel or input!");
	  exit(0);	  
	}
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
      if((nit+1)==maxIter and b<nConv){
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
    if(b==(nConv-1) and randomBoot==0){
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
  

  FILE *fpQ=openFile(outfiles,".qopt");
  // nBoot + 1, because first value is estimated Q
  printDouble(Q,nBoot+nConv,nPop,highestLike,nConv,ref.populations,fpQ);
  fclose(fpQ);
  fprintf(stderr,"FIRST row of .qopt file is BEST estimated Q, rest are nBoot bootstrapping Qs\n");
  fprintf(flog,"FIRST row of .qopt file is BEST estimated Q, rest are nBoot bootstrapping Qs\n");

  
  // only if certain flag, print adjusted freqs
  if(printFreq>0){
    gzFile fpGz=openFileGz(outfiles,".fopt.gz");
    printDoubleGz(F_1stRun,ref.refSites,nPop,ref.id,ref.populations,fpGz);
    gzclose(fpGz);
  }
   
  /*
  for(int i=0;1&&i<dumpedFiles.size();i++){
    free(dumpedFiles[i]);
    }*/
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  

  // print to log file
  fprintf(flog, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(flog, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fclose(flog); 
  return 0;
    
 }

