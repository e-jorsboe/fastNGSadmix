#pragma once
#include <stdlib.h>

typedef struct{
  char** familyID; 
  char** individualID; 
  char** paternalID; 
  char** maternalID; 
  int* sex;
  double* phenotype;
  int individuals;
  
}famFile;

typedef struct{
  int* chr; 
  char** rs; 
  double* recombrate; 
  int* position; 
  int* sex;
  char* minor;
  char* major;
  char** id;
  int sites;
  
}bimFile;


typedef struct{
  size_t x;
  size_t y;
  unsigned char **d;
  famFile fam;
  bimFile bim;
}plink;

bimFile readBimFile(const char* fname);
unsigned char **readbed(const char* file, int nrow,int ncol);
plink *readplink(const char *str);
void kill_plink(plink *p);
