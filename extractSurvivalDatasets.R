#!/usr/bin/env Rscript
# extractSurvivalDatasets.R

## PURPOSE: Finds lung, colon, prostate, breast, and pancreatic cancer 
## expression datasets from NCBI GEO using GEOmetadb
##
## USAGE: Rscript extractSurvivalDatasets.R [output name] [keywords]
##
## [output name]: The path containing the name of the output file
## [keywords]: A string consisting of a list of keywords containing names 
##             of the cancer type, where each keyword is separated by a |
##
##             Or, a default cancer type, where keywords have already been
##             determined, can be specified as the second argument in place
##             of a list of keywords. Below is a list of defaults that can
##             be entered in place of keywords:
##                  
##                  lung_cancer
##                  colon_cancer
##                  prostate_cancer
##                  breast_cancer
##                  pancreatic_cancer
##
## Example 1: Rscript extractSurvivalDatasets.R colon_cancer_survival "colon cancer|colon carcinoma|colon adenocarcinoma|colonic cancer|colonic carcinoma|colonic adenocarcinoma|rectal cancer|rectal carcinoma|rectal adenocarcinoma|CRC|COAD"
## Example 2: Rscript extractSurvivlDatasets.R lung_cancer_survival lung_cancer

args = commandArgs(trailingOnly=TRUE)

output = args[1]
keywords = args[2]

library(GEOmetadb)
library(stringr)

# extract data from NCBI GEO
# ** note: this will create a large file (approximately 7.28 GB) in the working directory
if (!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
con <- dbConnect(SQLite(), 'GEOmetadb.sqlite')


# extract all sample members with characteristics containing survival keywords
samples_survival <- dbGetQuery(con, paste("SELECT DISTINCT *", "FROM gsm", "WHERE characteristics_ch1 LIKE '%rfs%'", "OR characteristics_ch1 LIKE '%dfs%'", "OR characteristics_ch1 LIKE '%surv%'", "OR characteristics_ch1 LIKE '%dead%'", "OR characteristics_ch1 LIKE '%death%'", "OR characteristics_ch1 LIKE '%os[_]%'", "OR characteristics_ch1 LIKE '%pfs%'", "OR characteristics_ch1 LIKE '%dss%'", "OR characteristics_ch1 LIKE '%alive%'", "OR characteristics_ch1 LIKE '%os[:]%'", "OR characteristics_ch1 LIKE '%vital status%'", "OR characteristics_ch1 LIKE '%dfi%'", "OR characteristics_ch1 LIKE '%pfi%'", "AND organism_ch1='Homo sapiens'"))

# rename column from title to sample_title to distinguish gsm title from gse title
names(samples_survival)[2] <- "sample_title"

# gseWithTitle contains all gse accession numbers, the title, and the type of study
gseWithTitle <- dbGetQuery(con, paste("SELECT DISTINCT gse, title, type", "FROM gse"))

# gplWithTitle contains all gpl accession numbers, the title, and the manufacturer
gplWithTitle <- dbGetQuery(con, paste("SELECT DISTINCT gpl, title, manufacturer", "FROM gpl"))

# rename column from title to gpl_title to distinguish from gse title
names(gplWithTitle)[2] <- "gpl_title"

# gse2gsm contains all gse_gsm pairings
gse2gsm <- dbGetQuery(con, paste("SELECT *", "FROM gse_gsm"))
# filter for all the gse_gsm pairs with survival data
surv_gse_gsm = gse2gsm[gse2gsm$gsm %in% samples_survival$gsm, ]
# merge with gseWithTitle to include title and type
surv_gse_gsm = merge(surv_gse_gsm, gseWithTitle)
# put gse in a data frame to count number of samples per gse, name cols
gse_freq = as.data.frame(table(surv_gse_gsm$gse))
colnames(gse_freq) = c("gse", "sample_count")
# order by decreasing count, filter for gse with at least 50 samples
gse_freq = gse_freq[order(gse_freq$sample_count, decreasing = T), ]
gse_freq = gse_freq[gse_freq$sample_count >= 50, ]
# merge with gseWithTitle to include title and type
gse_title = merge(gse_freq, gseWithTitle)

# gse2gpl contains all gse_gpl pairings
gse2gpl <-dbGetQuery(con, paste("SELECT *", "FROM gse_gpl"))
# filter for all the gse_gpl pairs with survival data
surv_gse_gpl = gse2gpl[gse2gpl$gse %in% gse_freq$gse, ]
# merge with gplWithTitle to include gpl title and manufacturer
surv_gse_gpl = merge(surv_gse_gpl, gplWithTitle)
# merge with gse_count to add gpl info
gse_info = merge(gse_title, surv_gse_gpl)

# Filter for only expression data
gse_expression = subset(gse_info, grepl("Expression profiling by array|Expression profiling by high throughput sequencing|Non-coding RNA profiling by array|Non-coding RNA profiling by high throughput sequencing", type))
# Get rid of \t to avoid confusion when switching to tsv
gse_expression$type = gsub("\t", " ", gse_expression$type)

# filter for cancer type using default or custom keywords
if (keywords == "lung_cancer") {
    survival_cancer = subset(gse_expression, grepl("lung cancer|lung carcinoma|lung adenocarcinoma|SCLC|lung squamous|LUAD|LUSC", title, ignore.case=TRUE))
} else if (keywords == "colon_cancer") {
    survival_cancer = subset(gse_expression, grepl("colon cancer|colon carcinoma|colon adenocarcinoma|colonic cancer|colonic carcinoma|colonic adenocarcinoma|rectal cancer|rectal carcinoma|rectal adenocarcinoma|CRC|COAD|READ", title, ignore.case=TRUE))
} else if (keywords == "prostate_cancer") {
    survival_cancer = subset(gse_expression, grepl("prostate cancer|prostatic carcinoma|prostatic adenocarcinoma|prostate carcinoma|PRAD|prostate adenocarcinoma", title, ignore.case=TRUE))
} else if (keywords == "breast_cancer") {
    survival_cancer = subset(gse_expression, grepl("breast cancer|breast carcinoma|BRCA|breast adenocarcinoma", title, ignore.case=TRUE))
} else if (keywords == "pancreatic_cancer") {
    survival_cancer = subset(gse_expression, grepl("pancreatic cancer|pancreatic carcinoma|pancreatic adenocarcinoma|PAAD", title, ignore.case=TRUE))
} else {
    survival_cancer = subset(gse_expression, grepl(keywords, title, ignore.case=TRUE))
}

# write cancer survival data to tsv file
write.table(survival_cancer, file=paste0(output,".tsv"), sep='\t', row.names=FALSE, quote=TRUE)
