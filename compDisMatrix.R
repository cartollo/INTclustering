
# install.packages("argparse")
# install.packages("corrr")
# install.packages("ggcorrplot")
# install.packages("FactoMineR")
# install.packages("drclust")
# install.packages("coca")
# install.packages("uwot")
# install.packages("umap")
# library(drclust)
# install.packages("Rtsne")
# install.packages("clue")
# install.packages("aricode")
# install.packages("dbscan")
library(Rtsne)
library(dbscan)
library(clue)
library(aricode)
library(pheatmap)
library(coca)
library(grid)
library(argparse)
library(cluster)
library(proxy)
library("FactoMineR")
library(ggcorrplot)
library('corrr')
library(uwot)
library(FNN)
# library(umap)
source("myRfunc.R", keep.source = TRUE)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })


# Create a parser object
parser <- ArgumentParser(description = "handle input filenames")

# Define arguments
parser$add_argument("-f","--full",type = "character", default = "null", help = "full input file name")
parser$add_argument("-c","--core", type = "character", default = "null", help = "core input file name")
parser$add_argument("-n","--clusnum", type = "integer", default = 3, help = "number of clusters")
parser$add_argument("-o","--outfolder", type = "character", default = ".", help = "output folder")
parser$add_argument("-i","--interactive", type = "integer", default = 0, help = "interactive mode (1 for interactive, 0 for command line)")

# Parse the arguments
args <- parser$parse_args()
clusnum<-args$clusnum

#interactive:
if(args$interactive==1){
  clusnum<-3
  args$full="results/euclidean_ward.D_NOWHITE_CZM_isclr/out_create_dismatrix_euclidean_ward.D_NOWHITE_CZM_iscore_isclr.rds"
  args$core<-"results/euclidean_ward.D_NOWHITE_CZM_isclr/out_create_dismatrix_euclidean_ward.D_NOWHITE_CZM_isclr.rds"
}

# args$full<-"out_create_dismatrix_euclidean_ward.D_clsnum_NOWHITE_isclr_isrelab.rds"
# args$core="out_create_dismatrix_euclidean_ward.D2_NOWHITE_CZM_iscore_isclr_isrelab.rds"

if(args$full!="null"){ isfull<-TRUE}else{isfull<-FALSE}
if(args$core!="null"){ iscore<-TRUE}else{iscore<-FALSE}

set.seed(123)

if(iscore){
  corefilename<-args$core
  coreres<-single_sample_analisys(filename=corefilename,clusnum=clusnum, noprint = FALSE, txtoutfilename=txtoutfilename, newfile=TRUE, folder=args$outfolder, iscore=TRUE, plot_result=TRUE, iskmeans=TRUE) 
  txtoutfilename<-coreres$txtoutfilename
}


if(isfull){
  fullfilename<-args$full
  fullres<-single_sample_analisys(filename=fullfilename,clusnum=clusnum, noprint=FALSE, txtoutfilename=txtoutfilename, newfile=FALSE, folder=args$outfolder, iscore=FALSE, plot_result=TRUE, iskmeans=TRUE)
}


if(isfull && iscore){
  compare_twoclustering(first = fullres$dendo_cut,second = coreres$dendo_cut, firstfilename=fullres$filename, secondfilename=coreres$filename, txtoutfilename=txtoutfilename, folder=args$outfolder)
}


dev.off()



