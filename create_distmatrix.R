# user defined variables
# possible cluster methods:
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0){
  clusnum<-0
}else{ 
  clusnum <-as.integer(args[1])
} #clusternum dato da primo argomento

# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
clus_method <- "ward.D"
# possible distances:
# "euclidean", "maximum","manhattan", "canberra","binary", "minkowski", "mahalanobis", "bray_curtis", "jaccard", "aitchison"
distance <- "euclidean"

# iscore <- TRUE  #usare solo il core
iscore <- FALSE  #usare solo il core

isclr <- TRUE  #usare CLR
# isclr <- FALSE  #usare CLR
# isrelab<-TRUE  #usare la relative abundance
isrelab<-FALSE  #usare la relative abundance
#"ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky" or "NOWHITE"
whitemethod<-"NOWHITE"
#possible zeroimpmethods: "GBM"=default,"SQ","BL","CZM","user" , "pseudocount", "skip"
#in teoria, se uzi CZM o dovresti impostare amche gli altri parametri quali frac e threshold
zeroimpmethod<-"CZM"
#fare pca per massimizzare varianza 
ispcoa<-FALSE
# ispcoa<-TRUE
ispca<-FALSE
# ispca<-TRUE

iscreatedatabase<-FALSE

isdebug<-0
options(error=traceback)

#su DRAP
# inputfile="/Users/ymac/myINT/myemicrain/mountemicrain/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData"
#local
inputfile="/Users/ymac/Library/CloudStorage/OneDrive-FONDAZIONEIRCCSISTITUTONAZIONALEDEITUMORI/File\ di\ Iacovacci\ Jacopo\ -\ EMICRAIN/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData"

#do some checks
if(isclr && (distance=="bray_curtis" || distance=="jaccard")){
  print(paste("you chosed  with CLR with ",distance,", and this is nosense, please, be careful and try again"))
  quit(save = "no")  
}
if(!(distance == "euclidean" || distance == "maximum" || distance == "manhattan" || distance == "canberra" || distance == "binary" || distance == "minkowski") && whitemethod!="NOWHITE"){
  print(paste("you chosed  with whitening with ",distance,", and this is nosense, please, be careful and try again"))
  quit(save = "no")  
}
if(ispcoa==TRUE && (iscore==FALSE || ispca==TRUE)){
  print("please check iscore, ispca and ispccoa, coglione")
  quit(save="no")
}

# other variables
outputprefix <- paste("out_create_dismatrix", distance, clus_method,whitemethod,zeroimpmethod, sep = "_")
if (iscore) {outputprefix <- paste(outputprefix, "iscore", sep = "_")}
if (isclr) {outputprefix <- paste(outputprefix, "isclr", sep = "_")}
if (isrelab) { outputprefix <- paste(outputprefix, "isrelab", sep = "_")}
if(ispca){outputprefix <- paste(outputprefix, "ispca", sep = "_")}
if(ispcoa){outputprefix <- paste(outputprefix, "ispcoa", sep = "_")}

tobesaved<-list(
  config=list(
    clus_method=clus_method,
    distance=distance,
    iscore=iscore,
    isclr=isclr,
    ispca=ispca,
    ispcoa=ispcoa,
    isrelab=isrelab,
    whitemethod=whitemethod,
    zeroimpmethod=zeroimpmethod,
    inputfile=inputfile
  )
)

# print(sessionInfo())

# logfile <- paste0(outputprefix,".log")

# to do only once (or after each R update)
# install.packages("pheatmap")
# install.packages("zCompositions")
# install.packages("MatchIt")
# install.packages("ecodist")
# install.packages("vegan")
# install.packages("whitening")
# install.packages("tidyverse")
# install.packages("compositions")
# install.packages("chemometrics")
# install.packages("corrr")
# install.packages("ggcorrplot")
# install.packages("FactoMineR")
# install.packages("factoextra")
# install.packages("vegan")

# options
options(keep.source = TRUE)

# load libraries
library(pheatmap)
library("FactoMineR")
library(ggcorrplot)
library(factoextra)
library('corrr')
library(ecodist)
library(vegan)
library("whitening")
library(chemometrics)
source("myRfunc.R", keep.source = TRUE)

load(file = inputfile, verbose = TRUE)
# object ML.data

ML.dis.baseline.gen <- ML.data$OTU$dis$gen$bas
ML.val.baseline.gen <- ML.data$OTU$val$gen$bas
# Merge discovery and validation dataset
ML.dis.baseline.gen$genus <- rownames(ML.dis.baseline.gen)
ML.val.baseline.gen$genus <- rownames(ML.val.baseline.gen)
ML.baseline.gen <- merge(ML.dis.baseline.gen, ML.val.baseline.gen, by = "genus", all = TRUE)
rownames(ML.baseline.gen) <- ML.baseline.gen$genus
ML.baseline.gen$genus <- NULL # cancello colonna ridondante

ML.baseline.gen[is.na(ML.baseline.gen)] <- 0 # isnan allora mette 0 (ci sono vari invalid number in tabella)
# se una cella ha meno di 10 conteggi assegna 0
ML.baseline.gen[ML.baseline.gen < 10] <- 0
# righe/colonne che contengono solo 0 o NA vengono rimosse
ML.baseline.gen <- ML.baseline.gen[
  rowSums(ML.baseline.gen, na.rm = TRUE) > 0,
  colSums(ML.baseline.gen, na.rm = TRUE) > 0
]

if(iscreatedatabase){
  createdb$countsall<-ML.baseline.gen
}

# relative abundances
ML.baseline.gen.relab <- (ML.baseline.gen / rep(colSums(ML.baseline.gen), each = nrow(ML.baseline.gen))) * 100

if(iscreatedatabase){
  createdb$relaball<-ML.baseline.gen.relab
}

# Calcolo % di campioni in cui il genere ha ab. relativa ≥ 2%
# ed è presente in almeno il 10% dei casi
if(iscore){
  genus <- rowSums(ML.baseline.gen.relab >= 2) > (round(0.10 * ncol(ML.baseline.gen.relab)))
}else{
  genus <- rowSums(ML.baseline.gen.relab >= 1) > (round(0.05 * ncol(ML.baseline.gen.relab)))
}

if(iscreatedatabase){
  core_relabundance<-2
  core_presence<-0.1
  full_relabundance<-0
  full_presence<-0.05
  core_genus <- rowSums(ML.baseline.gen.relab >= core_relabundance) > (round(core_presence * ncol(ML.baseline.gen.relab)))
  full_genus <- rowSums(ML.baseline.gen.relab >= full_relabundance) > (round(full_presence * ncol(ML.baseline.gen.relab)))
  createdb$core_genus<-core_genus
  createdb$full_genus<-full_genus
  createdb$core_relabundance<-core_relabundance
  createdb$core_presence<-core_presence
  createdb$full_relabundance<-full_relabundance
  createdb$full_presence<-full_presence
}

# # Filtro la matrice rel ab con solo i genus "core"
# ML.baseline.gen.core <- ML.baseline.gen.relab[core_genus, ]


genus_list <- rownames(ML.baseline.gen.relab)[genus]
if (isdebug) {
  genus_list # stampo nomi genus core
}

# prendo solo conteggi
if (isrelab) {
  ML.baseline.gen.count <- ML.baseline.gen.relab[genus, ]
}else{
  ML.baseline.gen.count <- ML.baseline.gen[genus, ]
}

if(iscreatedatabase){
  createdb$core_counts<-ML.baseline.gen.relab[core_genus, ]
  createdb$core_counts_relab<-ML.baseline.gen.[core_genus, ]
  createdb$full_counts<-ML.baseline.gen.relab[full_genus, ]
  createdb$full_counts_relab<-ML.baseline.gen.[full_genus, ]
  tmp_jaccard<-ML.baseline.gen[core_genus, ]
  createdb$jaccard_core_input<-as.matrix(tmp_jaccard > 0.1) * 1
  tmp_jaccard<-ML.baseline.gen[full_genus, ]
  createdb$jaccard_full_input<-as.matrix(tmp_jaccard > 0.1) * 1
}

# solo 1 e 0 per dire presenza e assenza, assenza è quando è <10% numero che si può tunare
if(distance=="jaccard"){
  ML.baseline.gen.binary <- as.matrix(ML.baseline.gen.count > 0.1) * 1
}

#zero imputation
distinputmatrix<-zero_imputation(ML.baseline.gen.count,method=zeroimpmethod)
if(iscreatedatabase){
  createdb$core_counts_zeroed<-zero_imputation(createdb$core_counts,method=zeroimpmethod)
  createdb$core_counts_relab_zeroed<-zero_imputation(createdb$core_counts_relab,method=zeroimpmethod)
  createdb$full_counts_zeroed<-zero_imputation(createdb$full_counts,method=zeroimpmethod)
  createdb$full_counts_relab_zeroed<-zero_imputation(createdb$full_counts_relab,method=zeroimpmethod)

  createdb$zeroimpmethod<-zeroimpmethod
}

#check if there are zeroes:
if(any(distinputmatrix == 0)){
  print(paste("WARNING: distinputmatrix zero_imputation done with",zeroimpmethod," produced a matrix with ",sum(distinputmatrix == 0), "zeroes!"))
}
if(iscreatedatabase){
  if(any(createdb$core_counts_relab_zeroed == 0)){
    print(paste("WARNING: createdb$core_counts_relab_zeroed zero_imputation done with",zeroimpmethod," produced a matrix with ",sum(distinputmatrix == 0), "zeroes!"))
  }
  if(any(createdb$core_counts_zeroed == 0)){
    print(paste("WARNING: createdb$core_counts_zeroed zero_imputation done with",zeroimpmethod," produced a matrix with ",sum(distinputmatrix == 0), "zeroes!"))
  }
  if(any(createdb$full_counts_relab_zeroed == 0)){
    print(paste("WARNING: createdb$full_counts_relab_zeroed zero_imputation done with",zeroimpmethod," produced a matrix with ",sum(distinputmatrix == 0), "zeroes!"))
  }
  if(any(createdb$full_counts_zeroed == 0)){
    print(paste("WARNING: createdb$full_counts_zeroed zero_imputation done with",zeroimpmethod," produced a matrix with ",sum(distinputmatrix == 0), "zeroes!"))
  }
}

## trasformazione CLR
if (isclr) {
  if(isdebug){mybrowser("prima di isclr")}
  distinputmatrix<-apply(log(distinputmatrix), 2, function(x) x - mean(x))
  distinputmatrix <- t(scale(t(distinputmatrix), center = TRUE, scale = TRUE))
  # distinputmatrix<-clr(distinputmatrix)
 }

 if(iscreatedatabase){
  tmp_clr<-apply(log(createdb$core_counts_relab_zeroed), 2, function(x) x - mean(x))
  tmp_clr<-t(scale(t(tmp_clr), center = TRUE, scale = TRUE))
  createdb$core_counts_relab_zeroed_clr<-tmp_clr

  tmp_clr<-apply(log(createdb$core_counts_zeroed), 2, function(x) x - mean(x))
  tmp_clr<-t(scale(t(tmp_clr), center = TRUE, scale = TRUE))
  createdb$core_counts_zeroed_clr<-tmp_clr

  tmp_clr<-apply(log(createdb$full_counts_relab_zeroed), 2, function(x) x - mean(x))
  tmp_clr<-t(scale(t(tmp_clr), center = TRUE, scale = TRUE))
  createdb$full_counts_relab_zeroed_clr<-tmp_clr

  tmp_clr<-apply(log(createdb$full_counts_zeroed), 2, function(x) x - mean(x))
  tmp_clr<-t(scale(t(tmp_clr), center = TRUE, scale = TRUE))
  createdb$full_counts_zeroed_clr<-tmp_clr


 }

if(whitemethod!="NOWHITE"){
  tmp_distinputmatrix<-whiten(t(distinputmatrix),method=whitemethod)
  distinputmatrix<-t(tmp_distinputmatrix)
}

if(iscreatedatabase){
  tmp_whiten<-whiten(t(createdb$core_counts_relab_zeroed_clr),method=whitemethod)
  createdb$core_counts_relab_zeroed_clr_whiten<-t(tmp_whiten)
  tmp_whiten<-whiten(t(createdb$core_counts_zeroed_clr),method=whitemethod)
  createdb$core_counts_zeroed_clr_whiten<-t(tmp_whiten)
  tmp_whiten<-whiten(t(createdb$full_counts_relab_zeroed_clr),method=whitemethod)
  createdb$full_counts_relab_zeroed_clr_whiten<-t(tmp_whiten)
  tmp_whiten<-whiten(t(createdb$full_counts_zeroed_clr),method=whitemethod)
  createdb$full_counts_zeroed_clr_whiten<-t(tmp_whiten)
  createdb$whitemethod<-whitemethod
}

if(ispca){
  #questo è per fare pca, non pcoa
  distinputpca <- princomp(x=t(distinputmatrix), scores=TRUE)
  threshold<-0.05 #voglio tenere variabili che contenogno almeno tot percento di varianza

  var_pc <- distinputpca$sdev^2
  prop_var <- var_pc / sum(var_pc) #varianza di ogni nuova variabile
  cum_prop_var <- cumsum(prop_var) #cumulativa di varianza
  k <- which(prop_var >= threshold)
  distinputmatrix <- t(distinputpca$scores[, 1:length(k), drop = FALSE])
  if(isdebug){
    summary(distinputpca)
    fviz_eig(distinputpca, addlabels = TRUE)
    fviz_pca_var(distinputpca, col.var = "black")
    fviz_pca_var(distinputpca, col.var = "cos2",gradient.cols = c("black", "orange", "green"),repel = TRUE)
  }
}

if(iscreatedatabase){
 pcathreshold<-0.05
 createdb$pcathreshold<-pcathreshold
 createdb$core_counts_zeroed_clr_whiten_pca<-applypca(createdb$core_counts_zeroed_clr_whiten,pcathreshold)
 createdb$core_counts_relab_zeroed_clr_whiten_pca<-applypca(createdb$core_counts_relab_zeroed_clr_whiten,pcathreshold)
 createdb$full_counts_zeroed_clr_whiten_pca<-applypca(createdb$full_counts_zeroed_clr_whiten,pcathreshold)
 createdb$full_counts_relab_zeroed_clr_whiten_pca<-applypca(createdb$full_counts_relab_zeroed_clr_whiten,pcathreshold)
}


if(isdebug){mybrowser("prima di distanza")}

if (distance == "mahalanobis") {
  print("mahalanobis distace")
    distmatrix <- MatchIt::mahalanobis_dist(~., t(distinputmatrix))
} else if (distance == "bray_curtis") {
  print("bray_curtis distace")
  distmatrix <- bcdist(t(distinputmatrix))#distanze calcolate su righe
} else if (distance == "jaccard") {
  print("jaccard")
  distmatrix <- vegdist(t(ML.baseline.gen.binary), method = "jaccard")
} else if (distance == "aitchison") {
  print("aitchison")
  distmatrix <- aDist(t(distinputmatrix))
} else if (distance == "euclidean" || distance == "maximum" || distance == "manhattan" || distance == "canberra" || distance == "binary" || distance == "minkowski") {
  print(paste(distance, "distance"))
  distmatrix <- dist(t(distinputmatrix), method = distance) 
}

if(iscreatedatabase){
  createdb$mahalanobis_core_zeroed_clr <- MatchIt::mahalanobis_dist(~., t(core_counts_zeroed_clr))
  createdb$mahalanobis_full_zeroed_clr <- MatchIt::mahalanobis_dist(~., t(full_counts_zeroed_clr))
  createdb$mahalanobis_core_zeroed_clr_whiten <- MatchIt::mahalanobis_dist(~., t(core_counts_zeroed_clr_whiten))
  createdb$mahalanobis_full_zeroed_clr_whiten <- MatchIt::mahalanobis_dist(~., t(full_counts_zeroed_clr_whiten))
  createdb$mahalanobis_core_zeroed_clr_whiten_pca <- MatchIt::mahalanobis_dist(~., t(core_counts_zeroed_clr_whiten_pca))
  createdb$mahalanobis_full_zeroed_clr_whiten_pca <- MatchIt::mahalanobis_dist(~., t(full_counts_zeroed_clr_whiten_pca))

  createdb$braycurtis_core <- bcdist(t(createdb$core_counts))
  createdb$braycurtis_full <- bcdist(t(createdb$full_counts))

  createdb$jaccard_core <- bcdist(t(createdb$jaccard_core_input))
  createdb$jaccard_full <- bcdist(t(createdb$jaccard_full_input))

  createdb$aitchison_core_zeroed_clr <- aDist(t(core_counts_zeroed_clr))
  createdb$aitchison_full_zeroed_clr <- aDist(t(full_counts_zeroed_clr))
  createdb$aitchison_core_zeroed_clr_whiten <- aDist(t(core_counts_zeroed_clr_whiten))
  createdb$aitchison_full_zeroed_clr_whiten <- aDist(t(full_counts_zeroed_clr_whiten))
  createdb$aitchison_core_zeroed_clr_whiten_pca <- aDist(t(core_counts_zeroed_clr_whiten_pca))
  createdb$aitchison_full_zeroed_clr_whiten_pca <- aDist(t(full_counts_zeroed_clr_whiten_pca))

  createdb$euclidean_core_zeroed_clr <- dist(t(core_counts_zeroed_clr),method = "euclidean")
  createdb$euclidean_full_zeroed_clr <- dist(t(full_counts_zeroed_clr),method = "euclidean")
  createdb$euclidean_core_zeroed_clr_whiten <- dist(t(core_counts_zeroed_clr_whiten),method = "euclidean")
  createdb$euclidean_full_zeroed_clr_whiten <- dist(t(full_counts_zeroed_clr_whiten),method = "euclidean")
  createdb$euclidean_core_zeroed_clr_whiten_pca <- dist(t(core_counts_zeroed_clr_whiten_pca),method = "euclidean")
  createdb$euclidean_full_zeroed_clr_whiten_pca <- dist(t(full_counts_zeroed_clr_whiten_pca),method = "euclidean")

  createdb$manhattan_core_zeroed_clr <- dist(t(core_counts_zeroed_clr),method = "manhattan")
  createdb$manhattan_full_zeroed_clr <- dist(t(full_counts_zeroed_clr),method = "manhattan")
  createdb$manhattan_core_zeroed_clr_whiten <- dist(t(core_counts_zeroed_clr_whiten),method = "manhattan")
  createdb$manhattan_full_zeroed_clr_whiten <- dist(t(full_counts_zeroed_clr_whiten),method = "manhattan")
  createdb$manhattan_core_zeroed_clr_whiten_pca <- dist(t(core_counts_zeroed_clr_whiten_pca),method = "manhattan")
  createdb$manhattan_full_zeroed_clr_whiten_pca <- dist(t(full_counts_zeroed_clr_whiten_pca),method = "manhattan")
}

print(paste("distmatrix size: ncol=",ncol(distmatrix)," nrow=",nrow(distmatrix)))

if(isdebug){mybrowser("distanza scelta")}

if(ispcoa){
  pcoaed<-wcmdscale(distmatrix,k=10, eig=TRUE,add=TRUE, x.ret = TRUE)
}


tmp_distmatrix <- as.dist(distmatrix) # TODO non sono sicuro serva per euclideo etc. o per aitchison
# clustering gerarchico
dendo <- hclust(tmp_distmatrix, method = clus_method)

if(iscreatedatabase){

}

#save stuff
tobesaved$results$distinputmatrix<-distinputmatrix
tobesaved$results$distmatrix<-distmatrix
tobesaved$results$dendo<-dendo
tobesaved$metadata$time<-Sys.time()
tobesaved$metadata$sessioninfo<-sessionInfo()
outputname <- paste0(outputprefix, ".rds")
saveRDS(tobesaved,file=outputname)

if(clusnum!=0){
  outputname <- paste0(outputprefix,"_clusnum", clusnum,".pdf")
  dendo_cut <- cutree(dendo, k = clusnum) 
  dendo_picture(dendo_cut=dendo_cut,clusnum=clusnum, distinputmatrix=distinputmatrix, clus_method = clus_method, dendo=dendo)
 }




