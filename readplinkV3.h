#pragma once
#include <stdlib.h>
#include <vector>
#include <string>


typedef struct{
  std::vector<std::string> familyID; 
  std::vector<std::string> individualID; 
  std::vector<std::string> paternalID; 
  std::vector<std::string> maternalID; 
  int* sex;
  double* phenotype;
  int individuals;
  
}famFile;

typedef struct{
  int* chr; 
  std::vector<std::string> rs; 
  double* recombrate; 
  int* position; 
  int* sex;
  char* minor;
  char* major;
  std::vector<std::string> id;
  int sites;
  
}bimFile;


typedef struct{
  size_t x;
  size_t y;
  unsigned char **d;
  famFile fam;
  bimFile bim;
}plink;

famFile readFamFile(const char* fname, int lines);
bimFile readBimFile(const char* fname, int lines);
unsigned char **readbed(const char* file, int nrow,int ncol);
plink readplink(const char *str);
void kill_plink(plink &p);
