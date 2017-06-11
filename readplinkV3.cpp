/*
  log:
  g++ readplinkV3.cpp -lz -lpthread  -O3 -o readplinkV3


  debug:
  g++ readplinkV3.cpp -lz -lpthread -ggdb -O3 -o readplinkV3


*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "readplinkV3.h"
#include <string>


void dallocFamFile(famFile &f){
  delete [] f.sex;
  delete [] f.phenotype;
  
}

void dallocBimFile(bimFile &b){
  
  delete [] b.chr;
  delete [] b.recombrate;
  delete [] b.position;
  delete [] b.minor;
  delete [] b.major;
  
  
}

// bed file reader modified from snpMatrix by Clayton and Hintak Leung 2007
/*
  output will be genotype[ind][nsites]
  g \in {0,1,2,3},counts of allele2, 3=NA/missing
*/
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

// reads fam file
famFile readFamFile(const char* fname, int lines) {

   FILE* pFile;
   char buffer [10000];
   const char *delims = "\t \n";
   pFile = fopen (fname , "r");
   if (pFile == NULL) {
     perror ("Error opening fam file");
     exit(0);
   }

   
   famFile fam;
   fam.individuals = lines;
   fam.sex = new int[lines];
   fam.phenotype = new double[lines];
   int famIndex=0;

   while(fgets(buffer,100, pFile) != NULL){
     
     char* familyID = strtok(buffer,delims);
     std::string familyIDstring(familyID, strlen(familyID));
     fam.familyID.push_back(familyIDstring);
       
     char* individualID = strtok(NULL,delims);
     std::string individualIDstring(individualID, strlen(individualID));
     fam.individualID.push_back(individualIDstring);
     
     char* paternalID = strtok(NULL,delims);
     std::string paternalIDstring(paternalID, strlen(paternalID));
     fam.paternalID.push_back(paternalIDstring);
     
     char* maternalID = strtok(NULL,delims);
     std::string maternalIDstring(maternalID, strlen(maternalID));
     fam.maternalID.push_back(maternalIDstring);
     
     fam.sex[famIndex] = atoi(strtok(NULL,delims));
     fam.phenotype[famIndex] = atof(strtok(NULL,delims));
     famIndex++;
				    
   } 

   
   fclose (pFile);
   
   return(fam);
}

// read bim file with chr_pos ids
bimFile readBimFile(const char* fname,int lines) {

  
  FILE* pFile;
  char buffer [10000];
  const char *delims = "\t \n";
  pFile = fopen (fname , "r");
  if (pFile == NULL) {
    perror ("Error opening fam file");
    exit(0);
  }
  
   
  bimFile bim;
  bim.sites = lines;
  bim.chr = new int[lines];
  
  bim.recombrate = new double[lines];
  bim.position = new int[lines];
  bim.minor = new char[lines];
  bim.major = new char[lines];
  // chr_pos id here

  int bimIndex=0;
  while(fgets(buffer,100, pFile) != NULL){

    char* chrChar = strtok(buffer,delims);
    bim.chr[bimIndex] = atoi(chrChar);
    
    char* rs = strtok(NULL,delims);
    std::string rsString(rs, strlen(rs));
    bim.rs.push_back(rsString);
    
    bim.recombrate[bimIndex] = atof(strtok(NULL,delims));
    char* posChar = strtok(NULL,delims);
    bim.position[bimIndex] = atoi(posChar);
    bim.minor[bimIndex] = strtok(NULL,delims)[0];
    bim.major[bimIndex] = strtok(NULL,delims)[0];
     
    std::string chrString(chrChar,strlen(chrChar));
    std::string posString(posChar,strlen(posChar));

    std::string bimID = chrString + "_" + posString;
    // for concatenating strings now also with alleles
    if( tolower(bim.minor[bimIndex]) < tolower(bim.major[bimIndex]) ){
      bimID = bimID + "_" + bim.minor[bimIndex] + "_" + bim.major[bimIndex];
    } else{
      bimID = bimID + "_" + bim.major[bimIndex] + "_" + bim.minor[bimIndex];
    }
    
    bim.id.push_back(bimID);     
    bimIndex++;
    
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




void kill_plink(plink &p){

  for(int i=0;i<p.x;i++){
    free(p.d[i]);
  }
  dallocBimFile(p.bim);
  dallocFamFile(p.fam);
  free(p.d);

}



plink readplink(const char *str){
  plink p;

  char *fname =(char*) malloc(strlen(str)+5);
  strcpy(fname,str);strcpy(fname+strlen(fname),".bim");
  int ncol = nlines(fname);

  bimFile bim;
  bim = readBimFile(fname,ncol);
  
  strcpy(fname+strlen(str),".fam");
  int nrow = nlines(fname);

  famFile fam;
  fam = readFamFile(fname,nrow);
  
  //  nrow=3;
  if(ncol==0||nrow==0){
    fprintf(stderr,"Invalid binary plink file! \n");
    exit(0);
  }
  strcpy(fname+strlen(str),".bed");
  unsigned char **dat = readbed(fname,nrow,ncol);
  p.x=nrow;
  p.y=ncol;
  p.d=dat;
  p.bim=bim;
  p.fam=fam;
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
