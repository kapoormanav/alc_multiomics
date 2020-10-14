#!/bin/bash

for file in {1..22}

do

cat <<EOF> run.matrix.eqtl.chr$file.R

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/fairfax.chr$file.CD14.47231.414.b.qced.txt", sep="");
snps_location_file_name = paste(base.dir, "/hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/peer_expression/fairfax.chr$file.CD14.47231.414.b.qced.location.txt", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "/hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/CD14.47231.414.b.PEER_20.expression.txt", sep="");
gene_location_file_name = paste(base.dir, "/hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/peer_expression/gene.location.grch37_ensembl_95.update.01.07.19.txt", sep="");

# Covariates file name
# Set to character() for no covariates
#covariates_file_name = paste(base.dir, "/sc/orga/scratch/kapoom02/RNASEQ_AA/covariate_aa.txt", sep="");
covariates_file_name = character()
# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1.0;
pvOutputThreshold_tra = 0.000000005;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;


## Load genotype data

snps = SlicedData$new();
snps\$fileDelimiter = "\t";      # the TAB character
snps\$fileOmitCharacters = "NA"; # denote missing values;
snps\$fileSkipRows = 1;          # one row of column labels
snps\$fileSkipColumns = 1;       # one column of row labels
snps\$fileSliceSize = 1000;      # read file in slices of 2,000 rows
snps\$LoadFile("/hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/fairfax.chr$file.CD14.47231.414.b.qced.txt");

## Load gene expression data

gene = SlicedData$new();
gene\$fileDelimiter = "\t";      # the TAB character
gene\$fileOmitCharacters = "NA"; # denote missing values;
gene\$fileSkipRows = 1;          # one row of column labels
gene\$fileSkipColumns = 1;       # one column of row labels
gene\$fileSliceSize = 1000;      # read file in slices of 2,000 rows
gene\$LoadFile("/hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/peer_expression/CD14.47231.414.b.PEER_20.expression.txt");

for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);

## Load covariates

cvrt = SlicedData$new();
cvrt\$fileDelimiter = " ";      # the TAB character
cvrt\$fileOmitCharacters = "NA"; # denote missing values;
cvrt\$fileSkipRows = 1;          # one row of column labels
cvrt\$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt\$LoadFile("/sc/orga/scratch/kapoom02/RNASEQ_AA/covariate_aa.txt");
}

## Run the analysis

snpspos = read.table("/hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/fairfax.chr$file.CD14.47231.414.b.qced.location.txt", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("/hpc/users/kapoom02/6/from_scratch_SMR/Fairfax_hrc/EQTL_analysis/peer_expression/gene.location.grch37_ensembl_95.update.01.07.19.txt", header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me\$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
cis_out <- me\$cis\$eqtls
write.table(cis_out, "cis.eqtls.fairfax.chr$file.CD14.47231.414.b.qced.07.01.19.txt", row.names=FALSE, quote=FALSE)
#cat('Detected distant eQTLs:', '\n');
#trans_out <- me\$trans\$eqtls
#write.table(trans_out, "trasn.eqtl.fairfax.chr1.CD14.47231.414.b.qced.txt")
#show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

EOF

done
