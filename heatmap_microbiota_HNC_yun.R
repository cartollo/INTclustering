#user defined variables
# clus_method="ward.D"
clus_method<-"ward.D2"
clusnum<-3
distance<-"mahalanobis"
# distance<-"bray_curtis"
# distance<-"jaccard"

#other variables
outputname<-paste("danogram",distance,clus_method,clusnum,sep="_")
outputname<-paste0(outputname,".png")


## -----------------------------------------------------------------------------
#to do only once (or after each R update)
# install.packages("pheatmap", lib="/home/yun/myRlibrary")
# install.packages("zCompositions", lib = "/home/yun/myRlibrary")
# install.packages("MatchIt", lib = "/home/yun/myRlibrary")
# install.packages("ecodist", lib = "/home/yun/myRlibrary")
# install.packages("vegan", lib = "/home/yun/myRlibrary")
# install.packages("whitening", lib = "/home/yun/myRlibrary")

#load libraries
library(pheatmap)
library(ecodist)
library(vegan)
library("whitening")


## -----------------------------------------------------------------------------
# imputation and clr transformation
#è incompleto, perchè metodo GBM è lo stesso se ci scrivo qualsiasi altra cosa, per ora è implementato il bayesian-multiplicative
clr_transformation <- function(count_matrix, zero_imputation_method = "GBM") {
  if (zero_imputation_method == "pseudocount") { #sostituisce tutti i zeri con dei 1
    count_matrix[count_matrix == 0] <- 1
    clr_table <- count_matrix
  } else {
    clr_table <- zCompositions::cmultRepl( #Bayesian-multiplicative replacement.
      X = count_matrix,
      method = zero_imputation_method,
      z.delete = FALSE   
    )
  }
    clr_table <- apply(log(clr_table), 2, function(x) x - mean(x))
  return(clr_table)
}


## -----------------------------------------------------------------------------
# load data MicroLearner HNC 16S/ put your path data
#load("/Users/michelesalerno/Library/CloudStorage/OneDrive-FONDAZIONEIRCCSISTITUTONAZIONALEDEITUMORI/Microbiota_in_HNC_Causal/data/ML_data_16S_HNC.RData")
# load("https://istitutonazionaledeitumori-my.sharepoint.com/personal/jacopo_iacovacci_istitutotumori_mi_it/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fjacopo%5Fiacovacci%5Fistitutotumori%5Fmi%5Fit%2FDocuments%2FEMICRAIN%2FDatasets%2FMicrolearner%2FHNC%2Fmicrobiome%2F16S&viewid=afc0a39d%2D5525%2D4f3a%2D8065%2D539ae364ae41&sharingv2=true&fromShare=true&at=9&CT=1769081090548&OR=OWA%2DNT%2DMail&CID=47df9aba%2D3549%2Dfeb8%2D6056%2D2243fe4cc2b3&FolderCTID=0x012000B369C4E12EC92A4A99B23E04543EF70C")
# load("/home/yun/EMICRAIN/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData", verbose = TRUE)
load(file = "/Users/ymac/myINT/myemicrain/mountemicrain/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData", verbose = TRUE)
# load("/Users/ymac/myINT/myemicrain/mountemicrain/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData", verbose = TRUE)
# object ML.data



# con <- gzcon(file("/Users/ymac/myINT/myemicrain/mountemicrain/Datasets/Microlearner/HNC/microbiome/16S/ML_data_16S_HNC.RData", "rb"))
# load(con)
# close(con)

## -----------------------------------------------------------------------------
ML.dis.baseline.gen <- ML.data$OTU$dis$gen$bas;
ML.val.baseline.gen <- ML.data$OTU$val$gen$bas;

# Merge discovery and validation dataset
ML.dis.baseline.gen$genus <- rownames(ML.dis.baseline.gen)
ML.val.baseline.gen$genus <- rownames(ML.val.baseline.gen)
ML.baseline.gen <- merge(ML.dis.baseline.gen, ML.val.baseline.gen, by = "genus", all = TRUE)
rownames(ML.baseline.gen) <- ML.baseline.gen$genus
ML.baseline.gen$genus <- NULL #cancello colonna ridondante

ML.baseline.gen[is.na(ML.baseline.gen)] <- 0 #isnan allora mette 0 (ci sono vari invalid number in tabella)
#se una cella ha meno di 10 conteggi assegna 0
ML.baseline.gen[ML.baseline.gen<10] <- 0 
# righe/colonne che contengono solo 0 o NA vengono rimosse 
ML.baseline.gen <- ML.baseline.gen[rowSums(ML.baseline.gen, na.rm = TRUE) > 0,
                     colSums(ML.baseline.gen, na.rm = TRUE) > 0]


## -----------------------------------------------------------------------------
# relative abundances
ML.baseline.gen.relab <- (ML.baseline.gen / rep(colSums(ML.baseline.gen), each = nrow(ML.baseline.gen)))*100


## -----------------------------------------------------------------------------
# Calcolo % di campioni in cui il genere ha ab. relativa ≥ 2% 
#ed è presente in almeno il 10% dei casi
core_genus <- rowSums(ML.baseline.gen.relab >= 2) > (round(0.10 * ncol(ML.baseline.gen.relab)))
# Filtro la matrice rel ab con solo i genus "core"
ML.baseline.gen.core <- ML.baseline.gen.relab[core_genus, ]


## -----------------------------------------------------------------------------
core_genus_list <- rownames(ML.baseline.gen.relab)[core_genus]
core_genus_list #stampo nomi genus core


## -----------------------------------------------------------------------------
#prendo solo conteggi
ML.baseline.gen.core.count <- ML.baseline.gen[core_genus,]
#solo 1 e 0 per dire presenza e assenza, assenza è quando è <10%
ML.baseline.gen.core.binary <- as.matrix(ML.baseline.gen.core.count > 0.1) * 1

## -----------------------------------------------------------------------------
##trasformazione CLR
clr <- clr_transformation(ML.baseline.gen.core.count, zero_imputation_method = "CZM") 


## -----------------------------------------------------------------------------
#prendo ogni riga di clr (ogni batterio/etc.) ci sottraggo la media e lo divido per std dev.
#funzione t è trasposta, qui faccio due volte trasposta per sfruttare funzione scale (https://rdrr.io/r/base/scale.html)
ML.baseline.gen.z <- t(scale(t(clr), center = TRUE, scale = TRUE))


## -----------------------------------------------------------------------------
#calcolo distanza di mahalanobis
# x <- MatchIt::mahalanobis_dist(~ ., t(ML.baseline.gen.z))

#x <- lower.tri(x, diag = FALSE)
#x <- x[lower.tri(x, diag = F)]
#x

# x_lower <- x * lower.tri(x, diag = FALSE) #prendo solo la parte di matrice sottostante i zeri
#x_lower #lo printo


## -----------------------------------------------------------------------------
browser()
if(distance=="mahalanobis"){
  print("mahalanobis distace")
  x <- MatchIt::mahalanobis_dist(~ ., t(ML.baseline.gen.z))
}else if(distance=="bray_curtis"){
  print("bray_curtis distace")
  x<-bcdist(ML.baseline.gen.core.count)
}else if(distance=="jaccard"){
  print("jaccard")
  x<-vegdist(ML.baseline.gen.core.binary,method="jaccard")
}


d_vect <- as.dist(x) 

# clustering gerarchico
hcc_bas_maha <- hclust(d_vect, method = clus_method)

patient_clustering_bas_maha <- cutree(hcc_bas_maha, k = clusnum) #divido dendogramma in 3 parti
table(patient_clustering_bas_maha)


## -----------------------------------------------------------------------------
color.divisions <- 100
my_colors <- colorRampPalette(c("navy", "white", "red"))(color.divisions)

# breaks simmetrici attorno a 0
# max_val <- max(abs(ML.baseline.gen.z), na.rm = TRUE)
# breaks <- seq(-max_val, max_val, length.out = color.divisions + 1)
max_val <- 4
breaks <- seq(-max_val, max_val, length.out = color.divisions + 1)

g_gen_bas_maha <- pheatmap(
  ML.baseline.gen.z,
  color = my_colors,
  cluster_rows = FALSE,
  cluster_cols = hcc_bas_maha,   
  clustering_method = clus_method,  
  scale = "none",
  # annotation_col = annotation_col,
  breaks = breaks,
  cutree_cols = 3,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = "black",
  annotation_colors = annotation_colors,
  annotation_legend = TRUE,
  legend = TRUE,
  na_col = "#DDDDDD",
  filename=outputname
)

g_gen_bas_maha #printo figura










