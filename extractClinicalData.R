#!/usr/bin/env Rscript
# extractClinicalData.R

## PURPOSE: Extracts and parses clinical data from a given GEO Series accession
## number (e.g. GSExxx) using GEOquery; provides a baseline tsv file that can
## be further edited depending on the specific dataset
##
## USAGE: Rscript extractClinicalData.R [GSExxx] [output path - optional]
##
## [GSExxx]: The GEO Series accession number beginning with GSE followed by the
##           number referencing the specific series (e.g. GSE31210)
## [output name]: The path containing the name of the output file
##
## Example 1: Rscript extractClinicalData.R GSE31210 GSE31210_clinical 
## Example 2: Rscript extractClinicalData.R GSE31210

library(GEOquery)
library(stringr)
args <- commandArgs(trailingOnly=TRUE)

accession <- args[1]
if (length(args)==1) {
    # default output file
    output <- args[1]
} else {
    output <- args[2]
}

gse <- getGEO(accession)
phenoData <- gse[[1]]@phenoData@data
ch1 <- phenoData[, grepl(':ch1', colnames(phenoData))]
for (col in 1:ncol(ch1)) {
    colnames(ch1)[col] <- gsub(':ch1', '', colnames(ch1)[col])
}
ch1 <- cbind(data.frame(phenoData$geo_accession, phenoData$title, phenoData$type), ch1)
colnames(ch1)[1] <- 'sample'
colnames(ch1)[2] <- 'patient'
colnames(ch1)[3] <- 'type'
write.table(ch1, file=paste0(output, '.tsv'), row.names=FALSE, sep='\t', quote=FALSE)
