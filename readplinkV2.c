/*
  log:
  g++ readplinkV2.cpp -lz -lpthread  -O3 -o readplinkV2


  debug:
  g++ readplinkV2.cpp -lz -lpthread -ggdb -O3 -o readplinkV2


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "readplinkV2.h"


void dallocFamFile(famFile &f){
  for(int i=0;i<f.individuals;i++){
    free(f.familyID[i]);
    free(f.individualID[i]);
    free(f.maternalID[i]);
    free(f.paternalID[i]);
  }
  delete [] f.sex;
  delete [] f.phenotype;
  delete [] f.familyID;
  delete [] f.individualID;
  delete [] f.maternalID;
  delete [] f.paternalID;
  
}

void dallocBimFile(bimFile &b){
  for(int i=0;i<b.sites;i++){
    free(b.rs[i]);
    free(b.id[i]);
  }

  delete [] b.chr;
  delete [] b.rs;
  delete [] b.recombrate;
  delete [] b.position;
  delete [] b.minor;
  delete [] b.major;
  delete [] b.id;
  
}


/*
  output will be genotype[ind][nsites]
  g \in {0,1,2,3},counts of allele2, 3=NA/missing
*/
//modified from snpMatrix by clayton and hintak leung 2007
unsigned char **readbed(const char* file, int nrow,int ncol) {
  int i;
  const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
  const unsigned char mask = '\x03';


  FILE *in = fopen(file, "r");
  if (!in){
    fprintf(stderr,"Couldn't open input file: %s\n", file);
    exit(0);
  }
  unsigned char start[3];
  if (fread(start, 1, 3, in)!=3){
    fprintf(stderr,"Failed to read first 3 bytes");
    exit(0);
  }
  if (start[0]!='\x6C' || start[1]!='\x1B'){
    fprintf(stderr,"Input file does not appear to be a .bed file (%X, %X)", 
	   start[0], start[1]);
    exit(0);
  }
  /* Create output object */

  
  unsigned char** results =(unsigned char**) calloc(nrow,sizeof(unsigned char*));
  for(i=0;i<nrow;i++)
    results[i] =(unsigned char*) calloc(ncol,sizeof(unsigned char));

  /* Read in data */

  int snp_major = start[2];
  int part=0, ij=0, j=0;i=0;
  while (1) {
    unsigned char byte;
    if (!part) {
      if (feof(in) || !fread(&byte, 1, 1, in)) {
	printf("Unexpected end of file reached");
	exit(0);
      }
      part = 4;
    }
    unsigned char code = byte & mask;
    byte = byte >> 2;
    part--;
    unsigned char tmp = recode[code];
    if(tmp==0)
      results[i][j] = 3;
    else
      results[i][j] =tmp-1;

    assert(results[i][j]>=0 && results[i][j]<=3);

    if (snp_major) {
      ij++;
      i++;
      if (i==nrow) {
	i = part = 0;
	j++;
	if (j==ncol)
	  break;
      }
    }	
    else {
      ij += nrow;
      j++;
      if (j==ncol){
	j = part = 0;
	i++;
	if (i==nrow)
	  break;
	ij = i;
      }
    }
  }
  fclose(in);
  return results;
}

famFile readFamFile(const char* fname) {

   FILE* pFile;
   char buffer [10000];
   const char *delims = "\t \n";
   pFile = fopen (fname , "r");
   if (pFile == NULL) perror ("Error opening file");
   std::vector<char*> tmp;
   while(1){
       if (fgets(buffer,100, pFile) == NULL) break;
       tmp.push_back(strdup(buffer));
   }
   int lines = tmp.size();
   
   famFile fam;
   fam.individuals = lines;
   fam.familyID = new char*[lines];
   fam.individualID = new char*[lines];
   fam.paternalID = new char*[lines];
   fam.maternalID = new char*[lines];
   fam.sex = new int[lines];
   fam.phenotype = new double[lines];
   
   
   //fprintf(stdout,"Vector %s, %s, l=%i\n",tmp[0],tmp[1],tmp.size());
   for(int i=0; i<lines;i++){
           
    fam.familyID[i] = strdup(strtok(tmp[i],delims));
    fam.individualID[i] = strdup(strtok(NULL,delims));
    fam.maternalID[i] = strdup(strtok(NULL,delims));
    fam.paternalID[i] = strdup(strtok(NULL,delims));
    fam.sex[i] = atoi(strtok(NULL,delims));
    fam.phenotype[i] = atof(strtok(NULL,delims));
      
    }
    
   fclose (pFile);
   
   return(fam);
}


bimFile readBimFile(const char* fname) {

   FILE* pFile;
   char buffer [10000];
   const char *delims = "\t \n";
   pFile = fopen (fname , "r");
   if (pFile == NULL) perror ("Error opening file");
   std::vector<char*> tmp;
   while(1){
       if (fgets(buffer,100, pFile) == NULL) break;
       tmp.push_back(strdup(buffer));
   }
   int lines = tmp.size();
   
   bimFile bim;
   bim.sites = lines;
   bim.chr = new int[lines];
   bim.rs = new char*[lines];
   bim.recombrate = new double[lines];
   bim.position = new int[lines];
   bim.minor = new char[lines];
   bim.major = new char[lines];
   // do chr_pos id here
   bim.id = new char*[lines];
   
   //fprintf(stdout,"Vector %s, %s, l=%i\n",tmp[0],tmp[1],tmp.size());
   for(int i=0; i<lines;i++){
           
     bim.chr[i] = atoi(strtok(tmp[i],delims));
     bim.rs[i] = strdup(strtok(NULL,delims));
     bim.recombrate[i] = atof(strtok(NULL,delims));
     bim.position[i] = atoi(strtok(NULL,delims));
     bim.minor[i] = strtok(NULL,delims)[0];
     bim.major[i] = strtok(NULL,delims)[0];
     
     char chr[3];
     char pos[10];
     
     sprintf(chr, "%d", bim.chr[i]);
     sprintf(pos, "%d", bim.position[i]);
     
     char *concat =(char*) malloc(strlen(pos)+5);
     strcpy(concat,chr);strcpy(concat+strlen(chr),"_");strcpy(concat+(strlen(chr)+1),pos);
     bim.id[i] = strdup(concat);
     free(concat);
     
     
   }
   
   fclose (pFile);
   
   return(bim);
}


 

int nlines(const char *fname){
  FILE *in =NULL;
  if(!((in=fopen(fname,"r")))){
      fprintf(stderr,"Problem opening file: %s\n",fname);
      return 0;
  }
  int c;
  int n=0;
  while ( (c=fgetc(in)) != EOF ) 
    if ( c == '\n' )  n++;
  fclose(in);
  return n;
}

void kill_plink(plink *p){
  int i;
  for(i=0;i<p->x;i++){
    free(p->d[i]);
  }
  dallocBimFile(p->bim);
  dallocFamFile(p->fam);
  free(p->d);
  free(p);
  p=NULL;
}

plink *readplink(const char *str){
  plink *p =(plink*) malloc(sizeof( plink));

  char *fname =(char*) malloc(strlen(str)+5);
  strcpy(fname,str);strcpy(fname+strlen(fname),".bim");

  
  int ncol = nlines(fname);
  bimFile bim;
  bim = readBimFile(fname);
  
  
  strcpy(fname+strlen(str),".fam");
  famFile fam;
  fam = readFamFile(fname);
  int nrow = nlines(fname);
  //  nrow=3;
  if(ncol==0||nrow==0)
    return 0;
  fprintf(stderr,"Done reading file: \'%s\' with dim ncol:%d\tnrow:%d\n",str,ncol,nrow);
  strcpy(fname+strlen(str),".bed");
  unsigned char **dat = readbed(fname,nrow,ncol);
  p->x=nrow;
  p->y=ncol;
  p->d=dat;
  p->bim=bim;
  p->fam=fam;
  free(fname);
  return p;
}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  int i,j;
  if(argc==1){
    fprintf(stderr,"Supply prefix for plink binary files\n");
    return 0;
  }
  char *fname =(char*) malloc(strlen(argv[1])+5);
  strcpy(fname,argv[1]);strcpy(fname+strlen(fname),".bim");
  int ncol = nlines(fname);
  
  strcpy(fname+strlen(argv[1]),".fam");
  int nrow = nlines(fname);

  if(ncol==0||nrow==0)
    return 0;
  fprintf(stderr,"ncol:%d\tnrow:%d\n",ncol,nrow);
  strcpy(fname+strlen(argv[1]),".bed");
  unsigned char **dat = readbed(fname,nrow,ncol);
#if 1
  for(i=0;i<ncol;i++){
    for(j=0;j<nrow;j++)
      fprintf(stdout,"%d ",dat[j][i]);
    fprintf(stdout,"\n");
  }
#endif

  for(i=0;i<nrow;i++)
    free(dat[i]);
  free(dat);
  free(fname);
  return 0;
}
#endif
