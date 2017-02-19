# fastNGSadmix

Program for infering admixture proportions and doing PCA with a single NGS sample. Inferences based on reference panel.

For how to use program go to:
http://www.popgen.dk/software/
 
Installation:
=====

git clone https://github.com/e-jorsboe/fastNGSadmix.git;

cd fastNGSadmix; make

For the R files the snpStats package is required, it can be obatined thus:
source("https://bioconductor.org/biocLite.R")
biocLite("snpStats")

It has been tested for R version 3.2.x and later.

The R scripts and the C++ program has been tested on a 8 GB linux system,
however if one wants to create larger reference panels (and thereby genotypes),
doing it on a server with more RAM would be advisable.
As the R scripts will take up a lot of RAM in that case.
