# fastNGSadmix

Program for infering admixture proportions and doing PCA with a single NGS sample. Inferences based on reference panel.

For how to use program go to:
http://www.popgen.dk/software/index.php/FastNGSadmix

The program has been published in Bioinformatics:

Emil Jørsboe, Kristian Hanghøj, Anders Albrechtsen;
fastNGSadmix: Admixture proportions and principal component analysis of a single NGS sample,
Bioinformatics, btx474, https://doi.org/10.1093/bioinformatics/btx474

Installation:
=====

git clone https://github.com/e-jorsboe/fastNGSadmix.git;

cd fastNGSadmix; make

For the R files the snpStats package is required, it can be obatined thus:

source("https://bioconductor.org/biocLite.R");
biocLite("snpStats")

It has been tested for R version 3.2.x and later.

The R scripts and the C++ program has been tested on a 8 GB linux system,
however if one wants to create larger reference panels (and thereby genotypes),
doing it on a server with more RAM would be advisable.
As the R scripts will take up a lot of RAM in that case.

iAdmix (Bansal, 2015):
=====

I have provided some script for running iAdmix with genotype likelihoods directly.

First download the iAdmix software.

Then a .beagle file can be converted to the iAdmix genotype likelihoods file format via beagle2GL2.R:

Rscript beagle2GL2.R test.beagle

And then run with the .GL2 file in iAdmix using runancestryV2.py:

python runancestryV2.py --freq Ref.txt --GL test.GL2 --out test --path .
