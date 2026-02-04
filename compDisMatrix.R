
# install.packages("argparse")
# install.packages("corrr")
# install.packages("ggcorrplot")
# install.packages("FactoMineR")
library("FactoMineR")
library(argparse)
library(cluster)
library(proxy)
library('corrr')
source("myRfunc.R", keep.source = TRUE)

# Create a parser object
parser <- ArgumentParser(description = "handle input filenames")

# Define arguments
parser$add_argument("-f","--full",type = "character", default = "null", help = "full input file name")
parser$add_argument("-c","--core", type = "character", default = "null", help = "core input file name")

# Parse the arguments
args <- parser$parse_args()

#interactive:
# args$full="out_create_dismatrix_bray_curtis_ward.D_clsnum_NOWHITE_nocore_noclr_norelab.rds"
# args$core="out_create_dismatrix_bray_curtis_ward.D_clsnum_NOWHITE_iscore_noclr_norelab.rds"

args$full="out_create_dismatrix_mahalanobis_ward.D_clsnum_NOWHITE_nocore_noclr_norelab.rds"
args$core="out_create_dismatrix_mahalanobis_ward.D_clsnum_NOWHITE_iscore_noclr_norelab.rds"

# print the arguments
cat("full filename:", args$full, "\n")
cat("core filename:", args$core, "\n")

#load files
fullmatrix <- readRDS(args$full)
corematrix <- readRDS(args$core)

# single_sample_analisys(args$full$results$dendo, 3,args$full$results$distmatrix)

clusnum<-2
debug<-TRUE

coresil<-single_sample_analisys(corematrix$results$dendo,clusnum,corematrix$results$distmatrix, debug)

# fullsil<-single_sample_analisys(fullmatrix$results$dendo,clusnum,fullmatrix$results$distmatrix, debug)

