
/*
  log:
  g++ fastNGSadmix.cpp -lz -lpthread  -O3 -o fastNGSadmix

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

//This is taken from here:
//http://blog.albertarmea.com/post/47089939939/using-pthread-barrier-on-mac-os-x
#ifdef __APPLE__

#ifndef PTHREAD_BARRIER_H_
#define PTHREAD_BARRIER_H_

#include <pthread.h>
#include <errno.h>


/////////////////////////////////////////////////////////////////////
// globale ting: functioner + variable(ikke god stil)

typedef int pthread_barrierattr_t;
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int tripCount;
} pthread_barrier_t;


int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count)
{
    if(count == 0)
    {
        errno = EINVAL;
        return -1;
    }

if(pthread_mutex_init(&barrier->mutex, 0) < 0)
    {
        return -1;
    }
    if(pthread_cond_init(&barrier->cond, 0) < 0)
    {
        pthread_mutex_destroy(&barrier->mutex);
        return -1;
    }

    barrier->tripCount = count;
    barrier->count = 0;

    return 0;
}
int pthread_barrier_destroy(pthread_barrier_t *barrier)
{
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
    return 0;
}

int pthread_barrier_wait(pthread_barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->mutex);
    ++(barrier->count);
    if(barrier->count >= barrier->tripCount)
    {
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return 1;
    }
    else
    {
        pthread_cond_wait(&barrier->cond, &(barrier->mutex));
        pthread_mutex_unlock(&barrier->mutex);
        return 0;
    }
}

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

//this struct contains nescerray information needed for threading
typedef struct{
  int ID; 
  int threadNumber; 
  int nThreads; 
  double **F_1;
  double **Q_1;
  double **F;
  double **Q;
  int start; 
  int stop;
  int startI; 
  int stopI;
  double ***prodA;
  double ***prodB;
  double nPop;
  double **genos;
  int nInd;
  int nSites;
  double lres;//this is the likelihood for a block of data. total likelihood is sum of lres.
}pars;

//to make life simple we make stuff relating to the threading global
pthread_t *threads = NULL;
pars * myPars= NULL; //den ovenfor definerede struct type


double **allocDouble(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}

double ***allocDouble3(size_t x,size_t y,size_t z){
  double ***ret= new double**[x];
  for(size_t i=0;i<x;i++){
    ret[i] = new double*[y];
    for(size_t j=0;j<y;j++)
      ret[i][j] = new double[z];
  }
  return ret;
}



void dallocDouble3(double ***ret,size_t x,size_t y){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++)
      delete[] ret[i][j];
  }

  for(size_t i=0;i<x;i++){
    delete ret[i];
  }

  delete ret;
}

void minus(double **fst,double **sec,size_t x,size_t y,double **res){
  //  fprintf(stderr,"x=%lu y=%lu\n",x,y);
  for(size_t i=0;i<x;i++)
    for(size_t j=0;j<y;j++){
      //  fprintf(stderr,"i=%lu j=%lu\n",i,j);
      res[i][j] = fst[i][j]-sec[i][j];
    }
}

double sumSquare(double **mat,size_t x,size_t y){
  double tmp=0;
  for(size_t i=0;i<x;i++)
    for(size_t j=0;j<y;j++)
      tmp += mat[i][j]*mat[i][j];
  return tmp;
}


double sumSquareMinus(double **mat1,double **mat2,size_t x,size_t y){
  double tmp=0;
  for(size_t i=0;i<x;i++)
    for(size_t j=0;j<y;j++){
      double tmp2 = mat1[i][j]-mat2[i][j];
      tmp += tmp2*tmp2;
    }
  return tmp;
}


//same as above but swapped.... to lazy to change code
void dalloc(double **ret,size_t x){
  for(size_t i=0;i<x;i++)
    delete [] ret[i] ;
  delete [] ret;
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
  // fprintf(stderr,"start=%f\n",start);
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
      //  fprintf(stderr,"%f %f %f\n",W0,W1,W2);
      if(0&&std::isnan(sum)){
	//fprintf(stderr,"PRE[%d]: W %f\t%f\t%f sum=%f\n",i,W0,W1,W2,sum);
	exit(0);
      }
    }

    p=sum/keepInd;
    // fprintf(stderr,"it=%d\tp=%f\tsum=%f\tkeepInd=%d\n",it,p,log(sum),keepInd);
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }



  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
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
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=(loglike[i*3+0])*(1-p)*(1-p);
      W1=(loglike[i*3+1])*2*p*(1-p);
      W2=(loglike[i*3+2])*p*p;
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //fprintf(stderr,"p=%f W %f\t%f\t%f sum=%f loglike: %f\n",p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
    }
    p=-999;
    // exit(0);
  }
  
  return(p);
}


void map2domainQEmil(double* Q, int K){
  
  double sum=0;
  for(int k=0;k<K;k++){
    if(Q[k]<errTol){
      Q[k] = errTol;
    }
    if(Q[k]>(1-errTol)){
      Q[k] = 1-errTol;
    }
    sum+=Q[k];
  }
  for(int k=0;k<K;k++){
    Q[k]=Q[k]/sum;
  }
}

// errTol 0.1

void map2domainFEmil(double** F, int nSites, int K){
  for(int s=0;s<nSites;s++)
    for(int k=0;k<K;k++){
      if(F[s][k]<errTol){
	if(s==1){
	  fprintf(stderr,"Error and freqs is: %f and error %f\n,",F[s][k],errTol);
	}
	F[s][k] = errTol;
      }
      if(F[s][k]>1-errTol){
	if(s==1){
	  fprintf(stderr,"1-Error and freqs is: %f and error %f\n,",F[s][k],1-errTol);
	}
	F[s][k] = 1-errTol;
      }
    }
  
}





///////////////////////////////////////
// Trying to do faster version


void getExpGandFstarEmilFast(double* Q, double** F, int nSites, double* nInd, int K,double **genos, double **F_1, double* Q_1, double** F_org){

  map2domainFEmil(F, nSites, K);
  double* hjSq = new double[nSites];
  double* hjinvSq = new double[nSites]; 
  double* aSum = new double[K];
  double* bSum = new double[K];

  double** sumG = allocDouble(nSites,K);
  double** a = allocDouble(nSites,K);
  double** b = allocDouble(nSites,K);

  double rowTop = 0.0;
  double rowPre = 0.0;
  double ksRowSum = 0.0;
  double ksInvRowSum = 0.0;

 // set to zero, because refects values of RAM created from
 for(int i=0;i<nSites;i++){
   hjSq[i]=0;
   hjinvSq[i]=0;
   
   for(int k=0 ; k<K ; k++){
     sumG[i][k] = 0;
     a[i][k] = 0;
     b[i][k] = 0;
     
   }
 }

 for(int i=0;i<nSites;i++){  
   // calculates top which genotype likelihood times prob of genotype, first col times 2, then 1, 0
   for(int k=0;k<K;k++){
     // this should be original refFreqs
     sumG[i][k] = F_org[i][k] * (nInd[k]*2);
     hjSq[i] += (F[i][k] * Q[k]);
     hjinvSq[i] += (1-F[i][k]) * Q[k];

   }
   rowTop = (genos[i][0]*hjSq[i]*hjSq[i])*2 + (genos[i][1]*2*hjSq[i]*hjinvSq[i]);
   rowPre = (genos[i][0] * hjSq[i] * hjSq[i] + genos[i][1] * 2 * hjSq[i] * hjinvSq[i] + genos[i][2] * hjinvSq[i] * hjinvSq[i]);
   
   for(int k=0 ; k<K ; k++){
     
     ksRowSum+=F[i][k] * Q[k];
     ksInvRowSum+= (1-F[i][k]) * Q[k];
   }
   if(ksRowSum < 10e-7){
     for(int k=0 ; k<K ; k++){
       a[i][k] = sumG[i][k];
     }
   }
   else{
     for(int k=0 ; k<K ; k++){
       a[i][k] = (rowTop / rowPre) * ((F[i][k] * Q[k]) / ksRowSum) + sumG[i][k];
     }
   }
   if(ksInvRowSum < 10e-7){
     for(int k=0 ; k<K ; k++){
       b[i][k] = (nInd[k]*2) - sumG[i][k];
     }
   }
   else{
     for(int k=0 ; k<K ; k++){
       b[i][k] = (2 - (rowTop / rowPre)) * ((1-F[i][k]) * Q[k] / ksInvRowSum) + ((nInd[k]*2) - sumG[i][k]);
     }
   }
   for(int k=0 ; k<K ; k++){
     F_1[i][k] = a[i][k] / (a[i][k] + b[i][k]);
   }
   
   hjSq[i]=0.0;
   hjinvSq[i]=0.0;
   rowTop = 0.0;
   rowPre = 0.0;
   ksRowSum = 0.0;
   ksInvRowSum = 0.0;
   
 }
 
 // do rowSums for each row and then calculate a & b, and if all_ks has rowsum of 0 set a 0, if all_ks_inv has rowsum of 0 set b 0
 map2domainFEmil(F_1, nSites, K);
 
for(int i=0 ; i<nSites ; i++){ 
  for(int k=0;k<K;k++){
   
    hjSq[i] += (F_1[i][k] * Q[k]);   
    hjinvSq[i] += (1-F_1[i][k]) * Q[k];
   
  }
  rowTop = (genos[i][0]*hjSq[i]*hjSq[i])*2 + (genos[i][1]*2*hjSq[i]*hjinvSq[i]);
  rowPre = (genos[i][0] * hjSq[i] * hjSq[i] + genos[i][1] * 2 * hjSq[i] * hjinvSq[i] + genos[i][2] * hjinvSq[i] * hjinvSq[i]);

 for(int k=0 ; k<K ; k++){
    
   ksRowSum+=(F_1[i][k] * Q[k]);
   ksInvRowSum+=(1-F_1[i][k]) * Q[k];
 }
 if(ksRowSum < 10e-7){
     for(int k=0 ; k<K ; k++){
       a[i][k] = 0.0;
     }
 }
 else{
   for(int k=0 ; k<K ; k++) {
     a[i][k] = (rowTop / rowPre) * (F_1[i][k] * Q[k] / ksRowSum);
     
     }
 }
 if(ksInvRowSum < 10e-7){
   for(int k=0 ; k<K ; k++){
     b[i][k] = 0.0;
   }
 }
 else{
   for(int k=0 ; k<K ; k++){
     b[i][k] = (2 - (rowTop / rowPre)) * ((1-F_1[i][k]) * Q[k] / ksInvRowSum);
     
   }
 }
 
 for(int k=0 ; k<K ; k++){
 
   // sum a and b first and then assignt Q
   aSum[k]+=a[i][k];
   bSum[k]+=b[i][k];     
     
   }
 
   // fprintf(stderr,"Q=%f\n,",Q_1[k]);
   rowTop = 0.0;
   rowPre = 0.0;
   ksRowSum = 0.0;
   ksInvRowSum = 0.0;
   
 }
 for(int k=0 ; k<K ; k++){
     Q_1[k] = (aSum[k] + bSum[k]) / (2.0*nSites);
 }
map2domainQEmil(Q_1, K);

dalloc(sumG,nSites);
dalloc(a,nSites);
dalloc(b,nSites);
delete[] hjSq;
delete[] hjinvSq;
delete[] aSum;
delete[] bSum;


}









/////////////////////////////
// taken function from R and pasted in here
// deleted again! and pasted yet again


void getExpGandFstarEmil2(double* Q, double** F, int nSites, double* nInd, int K,double **genos, double **F_1, double* Q_1, double** F_org){

  map2domainFEmil(F, nSites, K);
  double** sumG = allocDouble(nSites,K);
  double* hjSq = new double[nSites];
  double* hjinvSq = new double[nSites]; 
  double** all_ks = allocDouble(nSites,K);
  double** all_ks_inv = allocDouble(nSites,K);
  double** pre = allocDouble(nSites,3); // genotypes based on HWE
  
  double** a = allocDouble(nSites,K);
  double** b = allocDouble(nSites,K);
 
 // set to zero, because refects values of RAM created from
 for(int i=0;i<nSites;i++){
   hjSq[i]=0;
   hjinvSq[i]=0;
   
   for(int k=0 ; k<K ; k++){
     all_ks_inv[i][k] = 0;
     all_ks[i][k] = 0;
     sumG[i][k] = 0;
     a[i][k] = 0;
     b[i][k] = 0;
     if(k<3){
       pre[i][k] = 0;
     }
   }
 }

 double q = 0.0;
 double qInv = 0.0;
 
 for(int i=0;i<nSites;i++){
   for(int k=0;k<K;k++){
     // this should be original refFreqs
     sumG[i][k] = F_org[i][k] * (nInd[k]*2);
     q = (F[i][k] * Q[k]);
     qInv = (1-F[i][k]) * Q[k];
     hjSq[i] += q;    
     hjinvSq[i] += qInv;
     all_ks[i][k] = q;
     all_ks_inv[i][k] = qInv;
     q = 0.0;
     qInv = 0.0;
   }
 }  
 
 for(int i=0;i<nSites;i++){  
   // calculates top which genotype likelihood times prob of genotype, first col times 2, then 1, 0
   pre[i][0] = genos[i][0] * hjSq[i] * hjSq[i];
   pre[i][1] = genos[i][1] * 2 * hjSq[i] * hjinvSq[i];
   pre[i][2] = genos[i][2] * hjinvSq[i] * hjinvSq[i];
   hjSq[i]=0.0;
   hjinvSq[i]=0.0;
   
}
 
 // do rowSums for each row and then calculate a & b, and if all_ks has rowsum of 0 set a 0, if all_ks_inv has rowsum of 0 set b 0
 double rowTop = 0.0;
 double rowPre = 0.0;
 double ksRowSum = 0.0;
 double ksInvRowSum = 0.0;
 
for(int i=0 ; i<nSites ; i++){
  for(int k=0 ; k<K ; k++){
    if(k<3){
      rowTop+=pre[i][k] * (2-k);
      rowPre+=pre[i][k];
    }
    ksRowSum+=all_ks[i][k];
    ksInvRowSum+=all_ks_inv[i][k];
  }
  if(ksRowSum < 10e-7){
    for(int k=0 ; k<K ; k++){
       a[i][k] = sumG[i][k];
    }
  }
  else{
    for(int k=0 ; k<K ; k++){
      a[i][k] = (rowTop / rowPre) * (all_ks[i][k] / ksRowSum) + sumG[i][k];
    }
  }
  if(ksInvRowSum < 10e-7){
    for(int k=0 ; k<K ; k++){
      b[i][k] = (nInd[k]*2) - sumG[i][k];
    }
  }
  else{
    for(int k=0 ; k<K ; k++){
      b[i][k] = (2 - (rowTop / rowPre)) * (all_ks_inv[i][k] / ksInvRowSum) + ((nInd[k]*2) - sumG[i][k]);
    }
  }
  for(int k=0 ; k<K ; k++){
    F_1[i][k] = a[i][k] / (a[i][k] + b[i][k]);
    if(i==0){
      fprintf(stderr,"F adj=%f\n,",F_1[i][k]);

    }
  }
  
  rowTop = 0.0;
  rowPre = 0.0;
  ksRowSum = 0.0;
  ksInvRowSum = 0.0;
 }
 map2domainFEmil(F_1, nSites, K);
 
 q = 0.0;
 qInv = 0.0;
 
 for(int i=0;i<nSites;i++){
   for(int k=0;k<K;k++){
     
     q = (F_1[i][k] * Q[k]);
     qInv = (1-F_1[i][k]) * Q[k];
     hjSq[i] += q;    
     hjinvSq[i] += qInv;
    
     all_ks[i][k] = q;
     all_ks_inv[i][k] = qInv;
     q = 0.0;
     qInv = 0.0;
   }
 }  

 for(int i=0;i<nSites;i++){
   
   // calculates top which genotype likelihood times prob of genotype, first col times 2, then 1, 0
   pre[i][0] = genos[i][0] * hjSq[i] * hjSq[i];
   pre[i][1] = genos[i][1] * 2 * hjSq[i] * hjinvSq[i];
   pre[i][2] = genos[i][2] * hjinvSq[i] * hjinvSq[i];
    
 }

 rowTop = 0.0;
 rowPre = 0.0;
 ksRowSum = 0.0;
 ksInvRowSum = 0.0; 
 
 for(int i=0 ; i<nSites ; i++){
   for(int k=0 ; k<K ; k++){
     if(k<3){
       rowTop+=pre[i][k] * (2-k);
       rowPre+=pre[i][k];
     }
    
     ksRowSum+=all_ks[i][k];
     ksInvRowSum+=all_ks_inv[i][k];
   }
   if(ksRowSum < 10e-7){
     for(int k=0 ; k<K ; k++){
       a[i][k] = 0.0;
     }
   }
   else{
     for(int k=0 ; k<K ; k++) {
       a[i][k] = (rowTop / rowPre) * (all_ks[i][k] / ksRowSum);
 
     }
   }
   if(ksInvRowSum < 10e-7){
     for(int k=0 ; k<K ; k++){
       b[i][k] = 0.0;
     }
   }
   else{
     for(int k=0 ; k<K ; k++){
       b[i][k] = (2 - (rowTop / rowPre)) * (all_ks_inv[i][k] / ksInvRowSum);
 
     }
   }
   
   rowTop = 0.0;
   rowPre = 0.0;
   ksRowSum = 0.0;
   ksInvRowSum = 0.0;
 }
 double aSum = 0.0;
 double bSum = 0.0;
 for(int k=0 ; k<K ; k++){
   for(int i=0; i<nSites;i++){
 
     // sum a and b first and then assignt Q
     aSum+=a[i][k];
     bSum+=b[i][k];     
     
   }
   Q_1[k] = (aSum + bSum) / (2.0*nSites);
   // fprintf(stderr,"Q=%f\n,",Q_1[k]);
   
   aSum = 0.0;
   bSum = 0.0;
 }

map2domainQEmil(Q_1, K);

dalloc(sumG,nSites);
dalloc(all_ks,nSites);
dalloc(all_ks_inv,nSites);
dalloc(pre,nSites);
dalloc(a,nSites);
dalloc(b,nSites);
delete[] hjSq;
delete[] hjinvSq;

}

////////////////////////

void getExpGandFstarEmil(double* Q, double** F, int nSites, double* nInd, int K,double **genos, double **F_1, double* Q_1, double** F_org){
  
  map2domainFEmil(F, nSites, K);
  double sumAG[K];
  double sumBG[K];
  double sumAGadj[K]; // for after freqs adjusted with input
  double sumBGadj[K];
  // do this for each site
  for(int k=0;k<K;k++){
    sumAGadj[k]=0;
    sumBGadj[k]=0;
  }
  for(int j=0;j<nSites;j++){   
    for(int k=0;k<K;k++){ //time killar
      sumAG[k]=0;
      sumBG[k]=0;
      
    }
    
    double fpart=0;
    double fpartInv=0;
    double expGG=0;
    for(int k=0;k<K;k++){ //time killar
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
      // K is number of ancestral populations

    }
    for(int k=0;k<K;k++){
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
    for(int k=0;k<K;k++){ //time killar
      // admixture adjusted freq, for each pop
         
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
      // K is number of ancestral populations
    }
    for(int k=0;k<K;k++){ //time killar
      sumAGadj[k] += expGG/(fpartAdj) * (Q[k] * F_1[j][k]); //proteckMe
      sumBGadj[k] += (2-expGG)/fpartAdjInv * (Q[k] * (1-F_1[j][k])); //proteckMe
    }
    
  }
  for(int k=0;k<K;k++){ //time killar
    Q_1[k]=(sumAGadj[k] + sumBGadj[k])/(2.0*nSites);
  }
  map2domainQEmil(Q_1,K);
  // should this be checked before using F_1??
  map2domainFEmil(F_1,nSites,K);
}



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
  char **ids;
  int nSites;
  int nInd;
  char **keeps; //matrix containing 0/1 indicating if data or missing
  int *keepInd; //keepInd[nSites] this is the number if informative samples
  float *mafs;
}bgl;

//utility function for cleaning up out datastruct
void dalloc(bgl &b){
  for(int i=0;i<b.nSites;i++){
    delete [] b.genos[i];
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
  const char *delims = "\t \n";
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  
  bgl ret;
  char buf[LENS];

  //find number of columns
  gzgets(fp,buf,LENS);
  strtok(buf,delims);
  int ncols=1;
  while(strtok(NULL,delims))
    ncols++;
  if(0!=(ncols-3) %3 ){
    fprintf(stderr,"ncols=%d\n",ncols);
    exit(0);
  }
  ret.nInd = (ncols-3)/3;//this is the number of samples
  
  //read every line into a vector
  std::vector<char*> tmp;
  while(gzgets(fp,buf,LENS))
    tmp.push_back(strdup(buf));
  
  //now we now the number of sites
  ret.nSites=tmp.size();
  ret.major= new char[ret.nSites];
  ret.minor= new char[ret.nSites];
  ret.ids = new char*[ret.nSites];
  ret.genos= new double*[ret.nSites];

  //then loop over the vector and parsing every line
  for(int s=0;SIG_COND&& (s<ret.nSites);s++){
    ret.ids[s] = strdup(strtok(tmp[s],delims));
    ret.major[s] =strtok(NULL,delims)[0];
    ret.minor[s] =strtok(NULL,delims)[0];
    ret.genos[s]= new double[3*ret.nInd];
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
  //  keeps=ret.keeps;
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
      //      fprintf(stderr,"i=%d j=%d\n",i,j);
      if(neg)
	d[i][j] = -atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
  fclose(fp);
}

void readDoubleGZ(double **d,int x,int y,const char*fname,int neg){
  fprintf(stderr,"opening : %s with x=%d y=%d\n",fname,x,y);
  const char*delims=" \n";
  gzFile fp = NULL;
  if((fp=gzopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000;
  char buf[lens];
  for(int i=0;i<x;i++){
      if(NULL==gzgets(fp,buf,lens)){
      fprintf(stderr,"Increase buffer\n");
      exit(0);
    }
    if(neg)
      d[i][0] = -atof(strtok(buf,delims));
    else
      d[i][0] = atof(strtok(buf,delims));
    for(int j=1;j<y;j++){
      //      fprintf(stderr,"i=%d j=%d\n",i,j);
      if(neg)
	d[i][j] = -atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
  gzclose(fp);
}

// Emil
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
    //      fprintf(stderr,"i=%d j=%d\n",i,j);   
    d[j] = atof(strtok(NULL,delims));
  }
  fclose(fp);
}



void printDoubleEmil(double *ret,size_t x,FILE *fp){
  for(size_t i=0;i<x;i++){
      fprintf(fp,"%.20f ",ret[i]);
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


int printer =0;


double likelihoodEmil(double* Q, double** F,int nSites, int K,double **genos){
  // F is sites times npop
  // Q is nInd times npop
  
  double prod_ind = 0.0;
  for(int j = 0; j < nSites; j++) {
    // extracts array of pointers - genos array of array of pointers
    double *gg = genos[j];
    double freq = 0.0;
    for(int k = 0; k < K; k++) {
      freq += (F[j][k])*Q[k];
    }
    double f = 1 - freq;
    double sum = gg[0] * f * f;
    sum += gg[1]*2*f*(1-f);
    sum += gg[2]*(1-f)*(1-f);
    prod_ind += log(sum);
  }
  
  
  
  return prod_ind;
}


//threaded
double likelihood(double** Q, double** F,int nSites, int nInd,int K,double **genos){
  // F is sites times npop
  // Q is nInd times npop
  double prod_sit = 0.0;
  for(int i = 0; i < nInd; i++) {
    double prod_ind = 0.0;
    for(int j = 0; j < nSites; j++) {
      double *gg = genos[j]+i*3;
      double freq = 0.0;
      for(int k = 0; k < K; k++) 
	freq += (F[j][k])*Q[i][k];
      double f =1- freq;
      double sum = gg[0] * f * f;
      sum += gg[1]*2*f*(1-f);
      sum += gg[2]*(1-f)*(1-f);
      prod_ind += log(sum);
    }
    prod_sit += prod_ind;     
  }

  #ifdef CHECK
    checkFQ(F,Q,nSites,nInd,K,__FUNCTION__);
  #endif
  return prod_sit;
}





// log likelihood for single site
// fills it directly into the structure
void likelihood_thd(double** Q, double** F,int nSites,int startI, int stopI,int K,double **genos,double &prod_sit){
  // F is sites times npop
  // Q is nInd times npop
  prod_sit =0;
  for(int i = startI; i < stopI; i++) {
    double prod_ind = 0.0;
    for(int j = 0; j < nSites; j++) {
      double *gg = genos[j]+i*3;
      double freq = 0.0;
      for(int k = 0; k < K; k++) 
	freq += (F[j][k])*Q[i][k];
      double f =1- freq;
      double sum = gg[0] * f * f;
      sum += gg[1]*2*f*(1-f);
      sum += gg[2]*(1-f)*(1-f);
      prod_ind += log(sum); //collect all site for a single individual
    }
    prod_sit += prod_ind;     
  }
  //checkFQ(F,Q,nSites,nInd,K,__FUNCTION__);
}
void *lkWrap(void *a){
  pars *p = (pars *)a; 
  likelihood_thd(p->Q,p->F,p->nSites,p->startI,p->stopI,p->nPop,p->genos,p->lres);
  return NULL;
}

// total log likelihood
double like_tsk(double **Q,double **F,int nThreads){
  for(int i=0;i<nThreads;i++){
    myPars[i].Q=Q;
    myPars[i].F=F;
    if( pthread_create(&threads[i],NULL,lkWrap,&myPars[i]))// here loglike is calculated in threads
      fprintf(stderr,"Problems starting threads\n");
  }
  for(int i=0;i<nThreads;i++)
    if(pthread_join(threads[i], NULL))
      fprintf(stderr,"problems joining\n");
  double res=0;
  for(int i=0;i<nThreads;i++){
    //    fprintf(stderr,"lres[%d]=%f\n",i,myPars[i].lres);
    res += myPars[i].lres;
  }
  return res;
}

void emEmilFast(double* Q, double** F, int nSites, double* nInd, int K,double **genos,double **F_1,double *Q_1, double** F_org) {

   if(Q==NULL){//cleanup
    
     return;
  }

   getExpGandFstarEmilFast(Q,F,nSites,nInd,K,genos,F_1,Q_1,F_org);

}


void emEmil2(double* Q, double** F, int nSites, double* nInd, int K,double **genos,double **F_1,double *Q_1, double** F_org) {

   if(Q==NULL){//cleanup
    
     return;
  }

   getExpGandFstarEmil2(Q,F,nSites,nInd,K,genos,F_1,Q_1,F_org);

}


// genos - genotype likelihoods
void emEmil(double* Q, double** F, int nSites, double* nInd, int K,double **genos,double **F_1,double *Q_1, double** F_org) {

   if(Q==NULL){//cleanup
    
     return;
  }

   getExpGandFstarEmil(Q,F,nSites,nInd,K,genos,F_1,Q_1,F_org);

}



pthread_barrier_t barr;
int dumpOld =0;

//Q,F are the old F_1,Q_1 are the next startpoint

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

  fprintf(stderr,"Setup:\n"); 
  fprintf(stderr,"\t-seed Seed for initial guess in EM\n"); 
  fprintf(stderr,"\t-P Number of threads\n"); 
  fprintf(stderr,"\t-method If 0 no acceleration of EM algorithm\n"); 
  fprintf(stderr,"\t-misTol Tolerance for considering site as missing\n");

  fprintf(stderr,"Stop chriteria:\n"); 
  fprintf(stderr,"\t-tolLike50 Loglikelihood difference in 50 iterations\n"); 
  fprintf(stderr,"\t-tol Tolerance for convergence\n"); 
  fprintf(stderr,"\t-dymBound Use dymamic boundaries (1: yes (default) 0: no)\n"); 
  fprintf(stderr,"\t-maxiter Maximum number of EM iterations\n"); 


  fprintf(stderr,"Filtering\n"); 
  fprintf(stderr,"\t-minMaf Minimum minor allele frequency\n"); 
  fprintf(stderr,"\t-minLrt Minimum likelihood ratio value for maf>0\n"); 
  fprintf(stderr,"\t-minInd Minumum number of informative individuals\n");


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

      //      fprintf(stderr,"diiffs=%f\n",fabs(d1[i][j]-d2[i][j]));
    if(fabs(d1[i]-d2[i])>diff){
      diff=fabs(d1[i]-d2[i]);
    }
  }
  return diff;
}


float calcThres(double **d1,double **d2, int x,int y){
  // finds the largest difference between 2 arrays
  // arrays has dimention x times y
  float diff=0;
  for(int i=0;i<x;i++)
    for(int j=0;j<y;j++){
      //      fprintf(stderr,"diiffs=%f\n",fabs(d1[i][j]-d2[i][j]));
      if(fabs(d1[i][j]-d2[i][j])>diff)
	diff=fabs(d1[i][j]-d2[i][j]);
    }
  return diff;
}

int whichMax(double *g){
  // equality signs such that when some are equal
  // the lowest genotype is preferred
  if(g[0]>=g[1]&&g[0]>=g[2])
    return 0;
  if(g[1]>g[0]&&g[1]>=g[2])
    return 1;
  else
    return 2;
}

void printLikes(bgl &d){
  // to write likelihoods for debugging
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

//void modLikesMinMaf(bgl &d,float minMaf){
void filterMinMaf(bgl &d,float minMaf){
  //  fprintf(stderr,"WARNING filtering minMaf=%f \n",minMaf);
  int posi =0;
  for(int s=0;s<d.nSites;s++){
    //    fprintf(stderr,"minmaf=%f mafs[%d]=%f lik=%f\n",minMaf,s,d.mafs[s],lik);
    if(d.mafs[s]>minMaf&&d.mafs[s]<1-minMaf){
      d.genos[posi] = d.genos[s];
      d.major[posi] = d.major[s];
      d.minor[posi] = d.minor[s];
      d.ids[posi] = d.ids[s];
      d.keeps[posi] = d.keeps[s];
      d.keepInd[posi] = d.keepInd[s];
      d.mafs[posi] = d.mafs[s];
      posi++;   
    }else{
      //      fprintf(stderr,"skippping\n");
    }
  }
  d.nSites=posi;
}

//void modLikesMiss(bgl &d,int minInd){
void filterMiss(bgl &d,int minInd){
  //  fprintf(stderr,"WARNING filtering mis=%d \n",minInd);
  int posi =0;
  for(int s=0;s<d.nSites;s++){
    if(d.keepInd[s]>minInd){
      d.genos[posi] = d.genos[s];
      d.major[posi] = d.major[s];
      d.minor[posi] = d.minor[s];
      d.ids[posi] = d.ids[s];
      d.keeps[posi] = d.keeps[s];
      d.keepInd[posi] = d.keepInd[s];
      d.mafs[posi] = d.mafs[s];
      posi++;   
    }else{
      // fprintf(stderr,"skippping\n");
    }
  }
  d.nSites=posi;
}

//void modLikesMinLrt(bgl &d,float minLrt){
void filterMinLrt(bgl &d,float minLrt){
  //  fprintf(stderr,"WARNING filtering minlrt=%f \n",minLrt);
  int posi =0;
  for(int s=0;s<d.nSites;s++){
    float lik=likeFixedMinor(d.mafs[s],d.genos[s],d.nInd,d.keeps[s]);
    float lik0=likeFixedMinor(0.0,d.genos[s],d.nInd,d.keeps[s]);
    //    fprintf(stderr,"minlrt=%f mafs[%d]=%f lik=%f lik0=%f 2*(lik0-lik)=%f\n",minLrt,s,d.mafs[s],lik,lik0,2.0*(lik0-lik));
    if(2.0*(lik0-lik)>minLrt){
      d.genos[posi] = d.genos[s];
      d.major[posi] = d.major[s];
      d.minor[posi] = d.minor[s];
      d.ids[posi] = d.ids[s];
      d.keeps[posi] = d.keeps[s];
      d.keepInd[posi] = d.keepInd[s];
      d.mafs[posi] = d.mafs[s];
      posi++;   
    }else{
      //  fprintf(stderr,"skippping\n");
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
  int minInd = 0;
  int printInfo = 0;
  float minMaf =0.05;
  float minLrt =0;
  const char* lname = NULL;
  const char* fname = NULL;
  const char* qname = NULL;
  const char* Nname = NULL;
  const char* outfiles = NULL;
  int nPop = 3;
  int seed =time(NULL);
  int nThreads = 1;
  // float tolLike50=0.1;
  float tolLike50=0.01; // changed by Emil
  // reading arguments
  argv++;
  while(*argv){
    // GL in the shape of beagle file
    if(strcmp(*argv,"-likes")==0 || strcmp(*argv,"-l")==0) lname=*++argv; //name / char arrays
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
    // ??
    else if(strcmp(*argv,"-P")==0) nThreads=atoi(*++argv); 
    else if(strcmp(*argv,"-printInfo")==0) printInfo=atoi(*++argv); 
    else if(strcmp(*argv,"-method")==0 || strcmp(*argv,"-m")==0) method=atoi(*++argv); 
    // different stop chriteria
    else if(strcmp(*argv,"-tolLike50")==0||strcmp(*argv,"-lt50")==0) tolLike50=atof(*++argv); //float/double - atof - char array to double/float
    // do I need those when I have only one individual??
    else if(strcmp(*argv,"-tol")==0||strcmp(*argv,"-t")==0) tol=atof(*++argv);
    else if(strcmp(*argv,"-maxiter")==0 || strcmp(*argv,"-i")==0) maxIter=atoi(*++argv); 
    // different filterings 
    else if(strcmp(*argv,"-misTol")==0 || strcmp(*argv,"-mt")==0) misTol=atof(*++argv);
    else if(strcmp(*argv,"-minMaf")==0||strcmp(*argv,"-maf")==0) minMaf=atof(*++argv);
    // min likelihood ratio value for maf>0
    else if(strcmp(*argv,"-minLrt")==0||strcmp(*argv,"-lrt")==0) minLrt=atof(*++argv);
    else if(strcmp(*argv,"-minInd")==0||strcmp(*argv,"-mis")==0) minInd=atoi(*++argv);
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
  if(lname==NULL){
    fprintf(stderr,"Please supply beagle file: -likes");
    info();
  }
  if(outfiles==NULL){
    fprintf(stderr,"Will use beagle fname as prefix for output\n");
    outfiles=lname;
  }

  //out put files
  FILE *flog=openFile(outfiles,".log");
  FILE *ffilter=openFile(outfiles,".filter");

  fprintf(stderr,"Input: lname=%s nPop=%d, fname=%s qname=%s outfiles=%s\n",lname,nPop,fname,qname,outfiles);
  fprintf(stderr,"Setup: seed=%d nThreads=%d method=%d\n",seed,nThreads,method);
  fprintf(stderr,"Convergence: maxIter=%d tol=%f tolLike50=%f dymBound=%d\n",maxIter,tol,tolLike50,dymBound);
  fprintf(stderr,"Filters: misTol=%f minMaf=%f minLrt=%f minInd=%d\n",misTol,minMaf,minLrt,minInd);

  fprintf(flog,"Input: lname=%s nPop=%d, fname=%s qname=%s outfiles=%s\n",lname,nPop,fname,qname,outfiles);
  fprintf(flog,"Setup: seed=%d nThreads=%d method=%d\n",seed,nThreads,method);
  fprintf(flog,"Convergence: maxIter=%d tol=%f tolLike50=%f dymBound=%d\n",maxIter,tol,tolLike50,dymBound);
  fprintf(flog,"Filters: misTol=%f minMaf=%f minLrt=%f minInd=%d\n",misTol,minMaf,minLrt,minInd);

  if(dymBound==0){
    errTolStart = errTolMin;
    errTol = errTolMin;
  }
    
  clock_t t=clock();//how long time does the run take
  time_t t2=time(NULL);

 
  //read BEAGLE likelihood file  
  // made into object to give it additional info
  bgl d=readBeagle(lname);
  fprintf(stderr,"Input file has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);
  fprintf(flog,"Input file has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);
  
  // printing name of first SNP out - First line of beagle treated as header
  fprintf(stderr,"first SNP has id %s\n",d.ids[0]);

  // prints out 3 first GL
  for(int i=0;i<3;i++){
    fprintf(stderr,"col %d has geno %f\n",i,d.genos[0][i]);
  }

  
  // filter sites based on MAF - 
  if(minMaf!=0.0)
    filterMinMaf(d,minMaf);
  if(minLrt!=0.0)
    filterMinLrt(d,minLrt);
  if(minInd!=0)
    filterMiss(d,minInd);
   if(printInfo)
    printKeepSites(d,ffilter);

  #ifdef DO_MIS
  keeps = d.keeps;
  #endif

  //  printLikes(d);
  fprintf(stderr,"Input file has dim (AFTER filtering): nsites=%d nind=%d\n",d.nSites,d.nInd);
  fprintf(flog,"Input file has dim (AFTER filtering): nsites=%d nind=%d\n",d.nSites,d.nInd);
  fflush(stderr);
  
  //set seed
  srand(seed);

  //unknown parameters
  // we only have a vector of Qs and then a matrix of freqs
  // so allocates a an array of pointers, each pointing to a second array
  // this gives matrix like structure of F
  double **F = allocDouble(d.nSites,nPop);
  double **F_org = allocDouble(d.nSites,nPop);
  double *Q = new double[nPop];
  double *N = new double[nPop];

  double **F_new = allocDouble(d.nSites,nPop);
  double *Q_new = new double[nPop];
  double *N_new = new double[nPop];

  //get start values
  readDoubleGZ(F,d.nSites,nPop,fname,0);
  readDoubleGZ(F_org,d.nSites,nPop,fname,0);
  // putting in an intial guess - matrix of Qs - change to vector
  
  double sum=0;
  for(int k=0;k<nPop;k++){
    break;
    // Q[k]=rand()*1.0/RAND_MAX;
    // Q[k]=1.0/K
    sum+=Q[k];
  }
  for(int k=0;k<nPop;k++){
    Q[k]=1.0/nPop;
    sum+=Q[k];
  }

  for(int k=0;k<nPop;k++) {
    // to make sure that proportions sum to 1
    Q[k]= Q[k]/sum; 
    Q_new[k]=Q[k];
  }
  
  
  if(Nname==NULL){
    fprintf(stderr,"Please supply number of individauls file: -Nname");
    info();
  }
  // reading initial Qs - matrix of Qs - change to vector
  readDouble1d(N,nPop,Nname);
  
  
  // Also managed to read in refFreqs
  for(int i=0;i<3;i++){
    fprintf(stderr,"Freq col %d has geno %f\n",i,F[0][i]);
    fprintf(stderr,"There are this many in pop%d %f\n",i,N[i]); 
    fprintf(stderr,"There is this ancestry in pop%d %f\n",i,Q[i]); 
  }

  //  double res =likelihood(Q, F, d.nSites, d.nInd, nPop,d.genos);    
  //  fprintf(stderr,"startres=%f\n",res);
  // square stuff beloq
  
  //  if(1) {
  
  //update the global stuff NOW
  //update the internal stuff in the pars for the threading
 

  double lold = -likelihoodEmil(Q, F, d.nSites, nPop,d.genos);

  fprintf(stderr,"iter[start] like is=%f\n",lold);

  //////////////////////////////////////// em ///////////////////////////////////  
  //below is the main looping trhought the iterations.
  // we have 4 possible ways, threading/nothreading line/noline
  int nit;
  double likeLast= lold;
  for(int nit=0; nit<100;nit++) {
    emEmil(Q, F, d.nSites, N, nPop,d.genos,F_new,Q_new,F_org);
    std::swap(Q,Q_new);
    std::swap(F,F_new);
    double lik = likelihoodEmil(Q, F, d.nSites, nPop,d.genos);
    likeLast=-lik;
  }  
  for(nit=1;SIG_COND&& nit<maxIter;nit++) {
    break;
    emEmil(Q, F, d.nSites, N, nPop,d.genos,F_new,Q_new,F_org);
    std::swap(Q,Q_new);
    std::swap(F,F_new);

    if((nit%3)==0 ){ //stopping criteria
      double lik = -likelihoodEmil(Q, F, d.nSites, nPop,d.genos);
      // thres is largest differense in admixture fractions
      fprintf(stderr,"iter[%d] like is=%f thres=%f\n",nit,lik,calcThresEmil(Q,Q_new,nPop));
      fprintf(stderr,"iter[%d] diff in likelihood is=%f",nit,std::abs(lik-likeLast));      
      //	fprintf(stderr,"iter[%d] like is=%f like old is=%f %f %f \n",nit,lik,likeLast, lik+likeLast, tolLike50);
      if(errTol>errTolMin){
	errTol=errTol/10;
	if(errTol<errTolMin)
	  errTol=errTolMin;
	//	fprintf(stderr,"50 changing errTol to %f\n",errTol);
      }
      else if(std::abs(lik-likeLast) < tolLike50){
	fprintf(stderr,"Convergence achived becuase log likelihooditer difference for 50 iteraction is less than %f\n",tolLike50);
	if(lik+likeLast<-1){
	  fprintf(stderr,"Convergence achived because log likelihooditer difference was NEGATIVE\n");
	}
	break;
      }
      
      likeLast=-lik;
      
    }
    
  }

  
  lold = likelihoodEmil(Q, F, d.nSites, nPop,d.genos);
  fprintf(stderr,"best like=%f after %d iterations\n",lold,nit);
  /////////////////////////////////////////////////////////// done - make output and clean /////////////////////////////////  
  // Print F and Q in files
  FILE *fp=openFile(outfiles,".qopt");
  printDoubleEmil(Q,nPop,fp);
  fclose(fp);
  
  gzFile fpGz=openFileGz(outfiles,".fopt.gz");
  printDoubleGz(F,d.nSites,nPop,fpGz);
  gzclose(fpGz);
  
  
  
  //deallocate memory 

  dalloc(F_new,d.nSites);
  delete [] Q_new;
  dalloc(d);
  delete [] threads;
  delete [] myPars;
  
  // if(nThreads==1){
  //  em(NULL, NULL, 0, 0, 0,NULL,NULL,NULL);
  // }


  for(int j = 0; j < d.nSites; j++) {
    delete[] F[j];
  }
  delete[] F;
  
  delete[] Q;
  for(int i=0;1&&i<dumpedFiles.size();i++){
    //    fprintf(stderr,"dumpedfiles are: %s\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  


  // print to log file

  fprintf(flog, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(flog, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(flog,"best like=%f after %d iterations\n",lold,nit);
  fclose(flog); 
  
  fclose(ffilter);
  return 0;
  
  
 }

