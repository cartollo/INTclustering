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
library(gplots)
source("myRfunc.R", keep.source = TRUE)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
set.seed(123)


######################## define my parameters ##################
debug<-0
clusnum<-3
# isfam<-TRUE
isfam<-FALSE
isshotgun<-TRUE
# isshotgun<-FALSE
#ouput file
if(isfam){
  folder<-"results/subsample_10_1000_correlation_average_16S_3clus_fam"
  reffile<-"results/correlation_average_NOWHITE_CZM_isclr_isfam/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_isfam_clusnum3_analysis.RDS"
}else{
  folder<-"results/subsample_10_1000_5clus_genus"
  reffile<-"results/correlation_average_NOWHITE_CZM_isclr/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_analysis.rds"
  # reffile="results/subsample_10_1000_16S_genus/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_seed549_analysis.RDS"
}
if(isshotgun){
  # folder<-"results/subsample_10_1000_shotgun_3clus_correlation_average/"
  folder<-"results/subsample_10_1000_spearman_average_3clus_shotgun"
  # folder<-"results/subsample_10_1000_euclidean_wardd2_3clus_shotgun"
  # reffile<-"results/euclidean_ward.D2_NOWHITE_CZM_isclr_shotgun/out_create_dismatrix_euclidean_ward.D2_NOWHITE_CZM_iscore_isclr_isrelab_shotgun_clusnum3_analysis.RDS"
  # reffile<-"results/correlation_average_NOWHITE_CZM_isclr_shotgun/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_shotgun_clusnum3_analysis.RDS"
  reffile<-"results/spearman_average_NOWHITE_CZM_isclr_3clus_shotgun/out_create_dismatrix_spearman_average_NOWHITE_CZM_iscore_isclr_isrelab_shotgun_clusnum3_analysis.RDS"
  # reffile<-"results/subsample_10_1000_correlation_average_3clus_shotgun/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_shotgun_seed99_analysis.RDS"
  # reffile<-"results/subsample_10_1000_3clus_shotgun/out_create_dismatrix_correlation_average_NOWHITE_CZM_iscore_isclr_shotgun_seed54_analysis.RDS"
}

outputfile<-paste0(folder,"/","comparefolderouput_create_dismatrix_spearman_average_NOWHITE_CZM_iscore_shotgun_clusnum_",clusnum,".pdf")
pdf(file = outputfile) 

databasea <- readRDS(reffile)
txtoutfilename<-NULL
plot_result<-FALSE
files <- list.files(path = folder, pattern="*\\_analysis.RDS$",full.names = TRUE)

#two file comparison stuff
ami_dendoclus<-c()
ari_dendoclus<-c()
overlap_dendoclus<-c() 
max_clus_overlap<-c() #valore di overlap per cluster con overlap massimo
min_clus_overlap<-c() #valore di overlap per cluster con overlap minimo
max_min_diff_clus_overlap<-c() #differnce between maximum and minimum overlap value
clus_size_diffmax<-c() #cluster con size massima - size minima
clus_size_max<-c() #size del cluster con size massima
clus_size_min<-c() #size del cluster con size minima
clus_size_median<-c() #size del cluster mediano
ari_kmean<-matrix(,nrow=length(databasea$results$res_umapdistmatrix), ncol=length(files))
ami_kmean<-matrix(,nrow=length(databasea$results$res_umapdistmatrix), ncol=length(files))
overlap_kmean<-matrix(,nrow=length(databasea$results$res_umapdistmatrix), ncol=length(files))

#single subfile analysis stuff
sub_overlap_kmeanwithdendo<-matrix(,nrow=length(databasea$results$res_umapdistmatrix), ncol=length(files))

filenum<-0
if(length(databasea$config$rdsoutfilename)>0){
  firstfilename<-databasea$config$rdsoutfilename
}else{
  firstfilename<-databasea$config$inputfile
}
for (bcase in files) {
  filenum<-filenum+1
  if(debug==1)
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
  total_counts<-sum(res$contingency_table)
  max_clus_overlap<-append(max_clus_overlap, max(diag(res$contingency_table))/total_counts)
  min_clus_overlap<-append(min_clus_overlap, min(diag(res$contingency_table))/total_counts)
  max_min_diff_clus_overlap<-append(max_min_diff_clus_overlap, (max(diag(res$contingency_table))-min(diag(res$contingency_table)))/total_counts)
  clus_size_diffmax<-append(clus_size_diffmax, max(colSums(res$contingency_table))-min(colSums(res$contingency_table)))
  clus_size_max<-append(clus_size_max, max(colSums(res$contingency_table)))
  clus_size_min<-append(clus_size_min, min(colSums(res$contingency_table)))
  clus_size_median<-append(clus_size_median, median(colSums(res$contingency_table)))

  #compare kmeans on umap 
  for (i in seq_along(databasea$results$res_umapdistmatrix)) {
    res<-compare_twoclustering(
      first = databasea$results$res_umapdistmatrix[[i]]$kmeans_end$km$cluster,
      second = databaseb$results$res_umapdistmatrix[[i]]$kmeans_end$km$cluster, 
      firstfilename=firstfilename, 
      secondfilename=secondfilename, 
      txtoutfilename=txtoutfilename, 
      folder=folder, 
      plot_result=plot_result,
      txtlabel=paste0("compare kmeans on umap between\n ", databasea$config$rdsoutfilename, " and ", databaseb$config$rdsoutfilename, "\n (kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue, " min_dist=", databasea$results$res_umapdistmatrix[[i]]$min_dist, " clusnum=",databasea$results$res_umapdistmatrix[[i]]$clusnum, ")")
    )
    ari_kmean[i,filenum]<-res$ARI
    ami_kmean[i,filenum]<-res$AMI
    overlap_kmean[i,filenum]<-res$overlap
    
    sub_overlap_kmeanwithdendo[i,filenum]<-  databaseb$results$res_umapdistmatrix[[i]]$compare_res$overlap
  }
}

#plot stuff
fitmethod<-NULL
# fitmethod<-"gaus"
create_histo(x=ami_dendoclus, main="ami for hierarchical clustering", xlab="ami", fitmethod = fitmethod, xlim=c(-0.1,1.1), breaks=120)
create_histo(x=ari_dendoclus, main="ari for hierarchical clustering", xlab="ari", fitmethod = fitmethod,xlim=c(-0.1,1.1), breaks=120)
create_histo(x=overlap_dendoclus, main="overlap for hierarchical clustering", xlab="overlap", fitmethod = fitmethod, xlim=c(-0.1,1.1), breaks=120)
create_histo(x=max_clus_overlap, main="maximum cluster overlap value for hierarchical clustering", xlab="max_clus_overlap", fitmethod = fitmethod, xlim=c(-0.1,1.1), breaks=120)
create_histo(x=min_clus_overlap, main="minimum cluster overlap value for hierarchical clustering", xlab="min_clus_overlap", fitmethod = fitmethod, xlim=c(-0.1,1.1), breaks=120)
create_histo(x=max_min_diff_clus_overlap, main="maximum - minimum cluster overlap value for hierarchical clustering", xlab="max_min_diff_clus_overlap", fitmethod = fitmethod, xlim=c(-0.1,1.1), breaks=120)
create_histo(x=clus_size_diffmax, main="maximum - minimum cluster size for hierarchical clustering", xlab="clus_size_diffmax", fitmethod = fitmethod, xlim=c(0,180), breaks=90)
create_histo(x=clus_size_max, main="cluster size maximum for hierarchical clustering", xlab="clus_size_max", fitmethod = fitmethod, xlim=c(0,180), breaks=90)
create_histo(x=clus_size_min, main="cluster size minimum for hierarchical clustering", xlab="clus_size_max", fitmethod = fitmethod, xlim=c(0,180), breaks=90)
if(clusnum>2)
  create_histo(x=clus_size_median, main="medan cluster size for hierarchical clustering", xlab="clus_size_median", fitmethod = fitmethod, xlim=c(0,180), breaks=90)

#facciamo un paio di 2d histos
hist2d(x=overlap_dendoclus,y=max_clus_overlap, xlab="overlap_dendoclus", ylab="max_clus_overlap")
hist2d(x=overlap_dendoclus,y=min_clus_overlap, xlab="overlap_dendoclus", ylab="min_clus_overlap")
hist2d(x=overlap_dendoclus,y=max_min_diff_clus_overlap, xlab="overlap_dendoclus", ylab="max_min_diff_clus_overlap")
hist2d(x=overlap_dendoclus,y=clus_size_diffmax, xlab="overlap_dendoclus", ylab="clus_size_diffmax")
hist2d(x=overlap_dendoclus,y=clus_size_max, xlab="overlap_dendoclus", ylab="clus_size_max")
hist2d(x=overlap_dendoclus,y=clus_size_min, xlab="overlap_dendoclus", ylab="clus_size_min")

if(clusnum>2)
  hist2d(x=overlap_dendoclus,y=clus_size_median, xlab="overlap_dendoclus", ylab="clus_size_median")

for (i in seq_along(databasea$results$res_umapdistmatrix)) {
  create_histo(x=ami_kmean[i,],main=paste0("ami for kmeans on umap kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue," min_dist=",databasea$results$res_umapdistmatrix[[i]]$min_dist), xlab="ami", fitmethod = fitmethod,xlim=c(-0.1,1.1), breaks=120)
  create_histo(x=ari_kmean[i,],main=paste0("ari for kmeans on umap kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue," min_dist=",databasea$results$res_umapdistmatrix[[i]]$min_dist), xlab="ari", fitmethod = fitmethod,xlim=c(-0.1,1.1), breaks=120)
  create_histo(x=overlap_kmean[i,],main=paste0("overlap for kmeans on umap kvalue=",databasea$results$res_umapdistmatrix[[i]]$kvalue," min_dist=",databasea$results$res_umapdistmatrix[[i]]$min_dist), xlab="overlap", fitmethod = fitmethod,xlim=c(-0.1,1.1), breaks=120)

  create_histo(x=sub_overlap_kmeanwithdendo[i,],main=paste0("overlap of kmeans on umap vs hierarchical=",databasea$results$res_umapdistmatrix[[i]]$kvalue," min_dist=",databasea$results$res_umapdistmatrix[[i]]$min_dist), xlab="overlap", fitmethod = NULL,xlim=c(-0.1,1.1), breaks=120)
}

cat("done, output file=",outputfile)

dev.off()






