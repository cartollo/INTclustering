#define variables
isdebug<-0
clusnum<-3

# install.packages("argparse")
# install.packages("corrr")
# install.packages("ggcorrplot")
# install.packages("FactoMineR")
# install.packages("drclust")
# install.packages("coca")
# library(drclust)
library(coca)
library(argparse)
library(cluster)
library(proxy)
library("FactoMineR")
library(ggcorrplot)
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
args$core<-"out_create_dismatrix_euclidean_ward.D_clsnum_NOWHITE_iscore_isclr_isrelab.rds"
args$full<-"out_create_dismatrix_euclidean_ward.D_clsnum_NOWHITE_isclr_isrelab.rds"
# args$core="out_create_dismatrix_bray_curtis_ward.D_clsnum_NOWHITE_iscore_noclr_norelab.rds"

if(args$full!="null"){ isfull<-TRUE}else{isfull<-FALSE}
if(args$core!="null"){ iscore<-TRUE}else{iscore<-FALSE}

if(isfull){
  fullfilename<-args$full
  fullsil<-single_sample_analisys(fullfilename,clusnum)
}
# if(iscore){
#   corefilename<-args$core
#   coresil<-single_sample_analisys(corefilename,clusnum) 
# }

  



