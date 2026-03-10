
source("create_distmatrix.R", keep.source = TRUE)

#possible clus_methods:
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

# possible distances:
# "euclidean", "maximum","manhattan", "canberra","binary", "minkowski", "mahalanobis", "bray_curtis", "jaccard", "aitchison", "pearson", "spearman", "correlation"

#possible whitemethods:
#"ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky" or "NOWHITE"

#possible zeroimpmethods: 
#"GBM"=default,"SQ","BL","CZM"=used by us,"user" , "pseudocount", "skip"
#in teoria, se uzi CZM o dovresti impostare amche gli altri parametri quali frac e threshold

  
#su DRAP
# inputfile="/Users/ymac/myINT/myemicrain/mountemicrain/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData"
#local
inputfile<-"/Users/ymac/Library/CloudStorage/OneDrive-FONDAZIONEIRCCSISTITUTONAZIONALEDEITUMORI/File\ di\ Iacovacci\ Jacopo\ -\ EMICRAIN/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData"
# inputfile<-"/Users/ymac/Library/CloudStorage/OneDrive-FONDAZIONEIRCCSISTITUTONAZIONALEDEITUMORI/File\ di\ Iacovacci\ Jacopo\ -\ EMICRAIN/Datasets/Microlearner/HNC/microbiome/shotgun/ML_data_shotgun_HNC.RData"

for (i in 1:1000){
  label<-paste0("_seed", i)
  create_db(seed=i, clusnum=0, distance="correlation", clus_method="average", iscore=TRUE, isweighted=FALSE, isclr=TRUE, ispca=FALSE, ispcoa=FALSE, isrelab=FALSE, whitemethod="NOWHITE", zeroimpmethod="CZM", iscreatedatabase=FALSE, inputfile=inputfile, isfam=TRUE, subsample_size=163, outputfilelable = label, folderoutput="results/subsample_10_1000_16S_5clus_fam")
}
  

