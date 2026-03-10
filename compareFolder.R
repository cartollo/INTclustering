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
library(fitdistrplus)
source("myRfunc.R", keep.source = TRUE)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
set.seed(123)


######################## define my parameters ##################
clusnum<-3
# isfam<-TRUE
isfam<-FALSE
isshotgun<-TRUE
# isshotgun<-FALSE
#ouput file
if(isfam){
  folder<-"results/subsample_10_1000_16S_4clus_fam"
  reffile<-"results/correlation_average_NOWHITE_CZM_isclr_isfam/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_isfam_analysis.RDS"
}else{
  folder<-"results/subsample_10_1000_5clus_genus"
  reffile<-"results/correlation_average_NOWHITE_CZM_isclr/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_analysis.rds"
  # reffile="results/subsample_10_1000_16S_genus/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_seed549_analysis.RDS"
}
if(isshotgun){
  folder<-"results/subsample_10_1000_3clus_shotgun"
  # reffile<-"results/correlation_average_NOWHITE_CZM_isclr_shotgun/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_shotgun_clusnum3_analysis.RDS"
  reffile<-"results/subsample_10_1000_3clus_shotgun/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_shotgun_seed54_analysis.RDS"
}
pdf(file = paste0(folder,"/","comparefolderouput_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_subvssub.pdf")) 

databasea <- readRDS(reffile)
txtoutfilename<-NULL
plot_result<-FALSE
files <- list.files(path = folder, pattern="*\\_analysis.RDS$",full.names = TRUE)
ami_dendoclus<-vector()
ari_dendoclus<-c()
overlap_dendoclus<-c()
ari_kmean<-matrix(,nrow=length(databasea$results$res_umapdistmatrix), ncol=length(files))
ami_kmean<-matrix(,nrow=length(databasea$results$res_umapdistmatrix), ncol=length(files))
overlap_kmean<-matrix(,nrow=length(databasea$results$res_umapdistmatrix), ncol=length(files))
filenum<-0
if(length(databasea$config$rdsoutfilename)>0){
  firstfilename<-databasea$config$rdsoutfilename
}else{
  firstfilename<-databasea$config$inputfile
}
for (bcase in files) {
  filenum<-filenum+1
  cat("Processing", bcase, "\n")
  data <- readRDS(bcase)
  databaseb <- readRDS(bcase)
  if(length(databaseb$config$rdsoutfilename)>0){
    secondfilename<-databaseb$config$rdsoutfilename 
  }else{
    secondfilename<-databaseb$config$inputfile 
  }  
  if(length(txtoutfilename)>0)
    cat("first file:",reffile, "\nsecond file:",bcase,file=txtoutfilename,append=TRUE)
  dendo_cuta <- cutree(databasea$results$dendo, k = clusnum)
  dendo_cutb <- cutree(databaseb$results$dendo, k = clusnum)
  #compare dendo clusteering
  res<-compare_twoclustering(
    first = dendo_cuta,
    second = dendo_cutb, 
    firstfilename=firstfilename, 
    secondfilename=secondfilename, 
    txtoutfilename=txtoutfilename, 
    folder=folder, 
    plot_result=plot_result,
    txtlabel=paste0("compare hierarchical clustering")
  )

  #fill stuff
  ami_dendoclus<-append(ami_dendoclus,res$AMI)
  ari_dendoclus<-append(ari_dendoclus,res$ARI)
  overlap_dendoclus<-append(overlap_dendoclus,res$overlap)

  #compare kmeans on umap 
  for (i in seq_along(databasea$results$res_umapdistmatrix)) {
    res<-compare_twoclustering(
      first = databasea$results$res_umapdistmatrix[[i]]$kmeans_end$km$cluster,
      second = databaseb$results$res_umapdistmatrix[[i]]$kmeans_end$km$cluster, firstfilename=firstfilename, 
      secondfilename=secondfilename, 
      txtoutfilename=txtoutfilename, 
      folder=folder, 
      plot_result=plot_result,
      txtlabel=paste0("compare kmeans on umap between\n ", databasea$config$rdsoutfilename, " and ", databaseb$config$rdsoutfilename, "\n (kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue, " min_dist=", databasea$results$res_umapdistmatrix[[i]]$min_dist, " clusnum=",databasea$results$res_umapdistmatrix[[i]]$clusnum, ")")
      )
      ari_kmean[i,filenum]<-res$ARI
      ami_kmean[i,filenum]<-res$AMI
      overlap_kmean[i,filenum]<-res$overlap
  }
}

#plot stuff
fitmethod<-"gaus"
create_histo(x=ami_dendoclus, main="ami for hierarchical clustering", xlab="ami", fitmethod = fitmethod, xlim=c(-0.1,0.1), breaks=40)
create_histo(x=ari_dendoclus, main="ari for hierarchical clustering", xlab="ari", fitmethod = fitmethod,xlim=c(-0.1,0.1), breaks=40)
create_histo(x=overlap_dendoclus, main="overlap for hierarchical clustering", xlab="overlap", fitmethod = fitmethod, xlim=c(0,1), breaks=50)
for (i in seq_along(databasea$results$res_umapdistmatrix)) {
  create_histo(x=ami_kmean[i,],main=paste0("ami for kmeans on umap kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue," min_dist=",databasea$results$res_umapdistmatrix[[i]]$min_dist), xlab="ami", fitmethod = fitmethod,xlim=c(-0.1,0.1), breaks=40)
  create_histo(x=ari_kmean[i,],main=paste0("ari for kmeans on umap kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue," min_dist=",databasea$results$res_umapdistmatrix[[i]]$min_dist), xlab="ari", fitmethod = fitmethod,xlim=c(-0.1,0.1), breaks=40)
  create_histo(x=overlap_kmean[i,],main=paste0("overlap for kmeans on umap kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue," min_dist=",databasea$results$res_umapdistmatrix[[i]]$min_dist), xlab="overlap", fitmethod = fitmethod,xlim=c(0.,1), breaks=50)
}




dev.off()






