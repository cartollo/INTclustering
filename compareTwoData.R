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
source("myRfunc.R", keep.source = TRUE)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })


# Create a parser object
parser <- ArgumentParser(description = "handle input filenames")

# Define arguments
parser$add_argument("-a","--acase",type = "character", default = "null", help = "first input file name")
parser$add_argument("-b","--bcase", type = "character", default = "null", help = "second input file name")
parser$add_argument("-n","--clusnum", type = "integer", default = 3, help = "number of clusters")
parser$add_argument("-o","--outfolder", type = "character", default = ".", help = "output folder")
parser$add_argument("-i","--interactive", type = "integer", default = 0, help = "interactive mode (1 for interactive, 0 for command line)")

# Parse the arguments
args <- parser$parse_args()
clusnum<-args$clusnum
isinteractive<-args$interactive
#da commentare o scommentare
isinteractive<-TRUE

#interactive:
if(isinteractive==TRUE){
  clusnum<-3
  args$acase="results/correlation_average_NOWHITE_CZM_isclr/out_create_dismatrix_correlation_average_NOWHITE_CZM_isclr_analysis.rds"
  args$bcase="results/correlation_average_NOWHITE_CZM_isclr/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_analysis.rds"
}

set.seed(123)

if(args$acase=="null" || args$bcase=="null")
  stop("Please provide both input files using -a and -b arguments.")

databasea <- readRDS(args$acase)
databaseb <- readRDS(args$bcase)

txtoutfilename<-"prova.txt"

outfilename=paste0(args$outfolder,"/",txtoutfilename)
cat("first file:",args$acase, "\nsecond file:",args$bcase,file=outfilename,append=TRUE)

dendo_cuta <- cutree(databasea$results$dendo, k = clusnum)
dendo_cutb <- cutree(databaseb$results$dendo, k = clusnum)
#compare dendo clusteering
res<-compare_twoclustering(
  first = dendo_cuta,
  second = dendo_cutb, 
  firstfilename=databasea$filename, 
  secondfilename=databaseb$filename, 
  txtoutfilename=txtoutfilename, 
  folder=args$outfolder, 
  txtlabel=paste0("compare hierarchical clustering")
)

#compare kmeans on umap 
for (i in seq_along(databasea$results$res_umapdistmatrix)) {
  res<-compare_twoclustering(
    first = databasea$results$res_umapdistmatrix[[i]]$kmeans_end$km$cluster,
    second = databaseb$results$res_umapdistmatrix[[i]]$kmeans_end$km$cluster, firstfilename=databasea$filename, 
    secondfilename=databaseb$filename, 
    txtoutfilename=txtoutfilename, 
    folder=args$outfolder, 
    txtlabel=paste0("compare kmeans on umap between\n ", databasea$config$rdsoutfilename, " and ", databaseb$config$rdsoutfilename, "\n (kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue, " min_dist=", databasea$results$res_umapdistmatrix[[i]]$min_dist, " clusnum=",databasea$results$res_umapdistmatrix[[i]]$clusnum, ")")
    )
}

# #compare kmeans on hierarchical clustering
# for (i in seq_along(databasea$results$res_umapdistmatrix)) {
#   compare_twoclustering(
# #    TODO: da sistemare qua per confronto tra robe fatte con metodi diverse
#     first = databasea$results$kmeans_end[i]$res_list$km$cluster,
#     second = databaseb$results$kmeans_end[i]$res_list$km$cluster, 
#     firstfilename=databasea$filename, 
#     secondfilename=databaseb$filename, 
#     txtoutfilename=databasea$txtoutfilename, 
#     folder=args$outfolder, 
#     txtlabel=paste0("compare hierarchical clustering between\n ", databasea$filename, "\n ", databaseb$filename)
#   )
# }


# dev.off()



