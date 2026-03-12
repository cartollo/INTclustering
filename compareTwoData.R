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
parser$add_argument("-f","--folder", type = "character", default = ".", help = "output folder")
parser$add_argument("-o","--outputfile", type = "character", default = "comparedOutput.pdf", help = "output file")
parser$add_argument("-i","--interactive", type = "integer", default = 0, help = "interactive mode (1 for interactive, 0 for command line)")

# Parse the arguments
args <- parser$parse_args()
clusnum<-args$clusnum
isinteractive<-args$interactive
#da commentare o scommentare
isinteractive<-FALSE
# isinteractive<-TRUE

#interactive:
if(isinteractive==TRUE){
  clusnum<-3
  args$acase<-"results/subsample_10_1000_3clus_shotgun/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_shotgun_seed998_analysis.RDS"
  args$bcase<-"results/subsample_10_1000_3clus_shotgun/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_shotgun_seed997_analysis.RDS"
  # args$bcase<-"results/correlation_average_NOWHITE_CZM_isclr_shotgun/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_shotgun_clusnum3_analysis.RDS"
}

set.seed(123)

if(args$acase=="null" || args$bcase=="null")
  stop("Please provide both input files using -a and -b arguments.")

databasea <- readRDS(args$acase)
databaseb <- readRDS(args$bcase)
pdf(file = paste0(args$folder,"/",args$outputfile))
dendo_cuta <- cutree(databasea$results$dendo, k = clusnum)
dendo_cutb <- cutree(databaseb$results$dendo, k = clusnum)
if(length(databasea$config$rdsoutfilename)>0){
  firstfilename<-databasea$config$rdsoutfilename
}else{
  firstfilename<-databasea$config$inputfile
}
if(length(databaseb$config$rdsoutfilename)>0){
  secondfilename<-databaseb$config$rdsoutfilename 
}else{
  secondfilename<-databaseb$config$inputfile 
}
#compare dendo clusteering
res<-compare_twoclustering(
  first = dendo_cuta,
  second = dendo_cutb, 
  firstfilename=firstfilename,
  secondfilename=secondfilename,
  folder=args$folder, 
  plot_result =  TRUE,
  txtlabel=paste0("compare hierarchical clustering")
)

#TODO: provv: 
mydendo_picturea<-dendo_picture(distinputmatrix=databasea$results$distinputmatrix, clus_method=databasea$config$clus_method, dendo=databasea$results$dendo, clusnum=clusnum, title="Dendrogram clustering a")
mydendo_pictureb<-dendo_picture(distinputmatrix=databaseb$results$distinputmatrix, clus_method=databaseb$config$clus_method, dendo=databaseb$results$dendo, clusnum=clusnum, title="Dendrogram clustering b")

dendo_cuta <- cutree(databasea$results$dendo, k = clusnum)
dendo_cutb <- cutree(databaseb$results$dendo, k = clusnum)
patientsa<-vector("list", length=clusnum)
patientsb<-vector("list", length=clusnum)
browser()

for(i in 1:clusnum){
  patientsa[[i]]<-sort(names(dendo_cuta)[dendo_cuta==i])
  patientsb[[i]]<-sort(names(dendo_cutb)[dendo_cutb==i])
  # linesa<-append(linesa,paste("cluster = ",i," patientes:",patientsa[i],sep=" "))
  # linesb<-append(linesb,paste("cluster = ",i," patientes:",patientsb[i],sep=" "))
}
linesa<-""
linesb<-""
for(i in 1:clusnum){
  linesa<-paste(linesa,patientsa[i], sep="\n")
  linesb<-paste(linesb,patientsb[i], sep="\n")
}
linesa<-paste0(linesa,"\n","check a caso: una sovvrapposizione dece essere di: ",length(intersect(names(dendo_cuta)[dendo_cuta == 3],names(dendo_cutb)[dendo_cutb == 1])),"\n") #provv
txt <- c(linesa)
grid::grid.newpage()
grid.text(
  # label=paste("pazienti a:",linesa,"pazienti b", linesb,sep="\n"),
  label=txt,
  x = 0.05, y = 0.95,
  just = c("left", "top"),
  gp = gpar(fontsize = 10)
)
txt <- c(linesb)
grid::grid.newpage()
grid.text(
  # label=paste("pazienti a:",linesa,"pazienti b", linesb,sep="\n"),
  label=txt,
  x = 0.05, y = 0.95,
  just = c("left", "top"),
  gp = gpar(fontsize = 10)
)


#compare kmeans on umap 
for (i in seq_along(databasea$results$res_umapdistmatrix)) {
  res<-compare_twoclustering(
    first = databasea$results$res_umapdistmatrix[[i]]$kmeans_end$km$cluster,
    second = databaseb$results$res_umapdistmatrix[[i]]$kmeans_end$km$cluster, 
    firstfilename=firstfilename, 
    secondfilename=secondfilename, 
    folder=args$folder, 
    plot_result = TRUE,
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
#     folder=args$folder, 
#     txtlabel=paste0("compare hierarchical clustering between\n ", databasea$filename, "\n ", databaseb$filename)
#   )
# }


dev.off()



