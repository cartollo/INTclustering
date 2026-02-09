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
iscore <- TRUE  #usare solo il core
# iscore <- FALSE  #usare solo il core
isclr <- TRUE  #usare CLR
# isclr <- FALSE  #usare CLR
isrelab<-TRUE  #usare la relative abundance
# isrelab<-FALSE  #usare la relative abundance
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
isdebug<-0
options(error=traceback)

inputfile="/Users/ymac/myINT/myemicrain/mountemicrain/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData"
#su DRAP
# inputfile="/home/yun/EMICRAIN/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData"

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
outputprefix <- paste("out_create_dismatrix", distance, clus_method, "clsnum",whitemethod, sep = "_")
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

print(sessionInfo())

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

# relative abundances
ML.baseline.gen.relab <- (ML.baseline.gen / rep(colSums(ML.baseline.gen), each = nrow(ML.baseline.gen))) * 100

# Calcolo % di campioni in cui il genere ha ab. relativa ≥ 2%
# ed è presente in almeno il 10% dei casi
if(iscore){
  core_genus <- rowSums(ML.baseline.gen.relab >= 2) > (round(0.10 * ncol(ML.baseline.gen.relab)))
}else{
  core_genus <- rowSums(ML.baseline.gen.relab >= 0) > (round(0.05 * ncol(ML.baseline.gen.relab)))
}

# # Filtro la matrice rel ab con solo i genus "core"
# ML.baseline.gen.core <- ML.baseline.gen.relab[core_genus, ]


core_genus_list <- rownames(ML.baseline.gen.relab)[core_genus]
if (isdebug) {
  core_genus_list # stampo nomi genus core
}

# prendo solo conteggi
if (isrelab) {
  ML.baseline.gen.count <- ML.baseline.gen.relab[core_genus, ]
}else{
  ML.baseline.gen.count <- ML.baseline.gen[core_genus, ]
}

# solo 1 e 0 per dire presenza e assenza, assenza è quando è <10% numero che si può tunare
if(distance=="jaccard"){
  ML.baseline.gen.binary <- as.matrix(ML.baseline.gen.count > 0.1) * 1
}

#zero imputation
distinputmatrix<-zero_imputation(ML.baseline.gen.count,method=zeroimpmethod)
#check if there are zeroes:
if(any(distinputmatrix == 0)){
  print(paste("WARNING: zero_imputation done with",zeroimpmethod," produced a matrix with ",sum(distinputmatrix == 0), "zeroes!"))
}

## trasformazione CLR
if (isclr) {
  if(isdebug){mybrowser("prima di isclr")}
  distinputmatrix<-apply(log(distinputmatrix), 2, function(x) x - mean(x))
  distinputmatrix <- t(scale(t(distinputmatrix), center = TRUE, scale = TRUE))
  # distinputmatrix<-clr(distinputmatrix)
 }

if(whitemethod!="NOWHITE"){
  tmp_distinputmatrix<-whiten(t(distinputmatrix),method=whitemethod)
  distinputmatrix<-t(tmp_distinputmatrix)
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

print(paste("distmatrix size: ncol=",ncol(distmatrix)," nrow=",nrow(distmatrix)))

if(isdebug){mybrowser("distanza scelta")}

if(ispcoa){
  pcoaed<-wcmdscale(distmatrix,k=10, eig=TRUE,add=TRUE, x.ret = TRUE)
}


tmp_distmatrix <- as.dist(distmatrix) # TODO non sono sicuro serva per euclideo etc. o per aitchison

# clustering gerarchico
dendo <- hclust(tmp_distmatrix, method = clus_method)

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
  dendo_picture(dendo,clusnum,outputname,distmatrix, distinputmatrix)
 }




