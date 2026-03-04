
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
parser$add_argument("-f","--first",type = "character", default = "null", help = "first input file name")
parser$add_argument("-s","--second", type = "character", default = "null", help = "second input file name")
parser$add_argument("-n","--clusnum", type = "integer", default = 3, help = "number of clusters")
parser$add_argument("-o","--outfolder", type = "character", default = ".", help = "output folder")
parser$add_argument("-i","--interactive", type = "integer", default = 0, help = "interactive mode (1 for interactive, 0 for command line)")

# Parse the arguments
args <- parser$parse_args()
clusnum<-args$clusnum
isinteractive<-args$interactive
#da commentare o scommentare
# isinteractive<-TRUE

#interactive:
if(isinteractive==TRUE){
  clusnum<-3
  args$first="results/pearson_average_NOWHITE_CZM_isclr/out_create_dismatrix_pearson_average_NOWHITE_CZM_isclr.rds"
  args$second<-"results/pearson_average_NOWHITE_CZM_isclr/out_create_dismatrix_pearson_average_NOWHITE_CZM_iscore_isclr.rds"
}

# args$full<-"out_create_dismatrix_euclidean_ward.D_clsnum_NOWHITE_isclr_isrelab.rds"
# args$core="out_create_dismatrix_euclidean_ward.D2_NOWHITE_CZM_iscore_isclr_isrelab.rds"

if(args$first!="null"){ isfirst<-TRUE}else{isfirst<-FALSE}
if(args$second!="null"){ issecond<-TRUE}else{issecond<-FALSE}

set.seed(123)

if(isfirst){
  fullfilename<-args$first
  firstres<-single_sample_analisys(filename=fullfilename,clusnum=clusnum, noprint=FALSE, txtoutfilename=txtoutfilename, newfile=TRUE, folder=args$outfolder, plot_result=TRUE, iskmeans=TRUE)
  txtoutfilename<-firstres$txtoutfilename
}
if(isfirst) newfile<-FALSE else newfile<-TRUE

if(issecond){
  secondfilename<-args$second
  secondres<-single_sample_analisys(filename=secondfilename,clusnum=clusnum, noprint = FALSE, txtoutfilename=txtoutfilename,newfile=newfile, folder=args$outfolder, plot_result=TRUE, iskmeans=TRUE) 
  if(!isfirst)  txtoutfilename<-secondres$txtoutfilename
}




if(isfirst && issecond){
  cat("Compare between full and core\n",file=txtoutfilename,append=TRUE)
  compare_twoclustering(first = firstres$dendo_cut,second = secondres$dendo_cut, firstfilename=firstres$filename, secondfilename=secondres$filename, txtoutfilename=txtoutfilename, folder=args$outfolder)
   for (i in seq_along(firstres$res_umapdistmatrix)) {
    compare_twoclustering(
      first = firstres$res_umapdistmatrix[[i]]$kmeans_end$km$cluster,
      second = secondres$res_umapdistmatrix[[i]]$kmeans_end$km$cluster, 
      firstfilename=paste0(firstres$filename, 
      "\nkmeans on umap kvalue=",firstres$res_umapdistmatrix[[i]]$kvalue, 
      " min_dist=", firstres$res_umapdistmatrix[[i]]$min_dist,
      " clusnum=",firstres$res_umapdistmatrix[[i]]$clusnum), 
      secondfilename=paste0(secondres$filename, 
      "\nkmeans on umap kvalue=",secondres$res_umapdistmatrix[[i]]$kvalue, 
      " min_dist=",secondres$res_umapdistmatrix[[i]]$min_dist, 
      " clusnum=",secondres$res_umapdistmatrix[[i]]$clusnum), txtoutfilename=txtoutfilename, 
      folder=args$outfolder,
      txtlabel=paste0("Compare kmeans on umap between first and second (kvalue=",firstres$res_umapdistmatrix[[i]]$kvalue, " min_dist=", firstres$res_umapdistmatrix[[i]]$min_dist, " clusnum=",firstres$res_umapdistmatrix[[i]]$clusnum, ")")
      )
  }
}


dev.off()



