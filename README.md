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
