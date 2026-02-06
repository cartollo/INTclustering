# imputation and clr transformation
#è incompleto, perchè metodo GBM è lo stesso se ci scrivo qualsiasi altra cosa, per ora è implementato il bayesian-multiplicative
zero_imputation <- function(count_matrix, method = "CZM") {
  if (method == "skip") { #sostituisce tutti i zeri con dei 1 
    return(count_matrix)
  } else if(method == "pseudocount") {
    count_matrix[count_matrix == 0] <- 1
    return (count_matrix)
  }else{#Bayesian-multiplicative replacement
    return( zCompositions::cmultRepl(X = count_matrix,method = method,z.delete = FALSE)) 
  }
}

#per fare debug, id è numero progressivo, msd è un messaggio che stampa
mybrowser <- function(msg = "", forcequit=FALSE) {
  file <- tryCatch(basename(sys.frame(1)$ofile), error = function(e) "console")
  message(sprintf("=== BROWSER %s | %s | %s ===", isdebug, file, msg))
  browser()
  if(forcequit==TRUE){
    quit(save="no")
  }
  isdebug=isdebug+1
}

dendo_picture<-function(dendo, clusnum=NA, outputname=NA){

  dendo_cut <- cutree(dendo, k = clusnum) # divido dendogramma in 3 parti
  table(dendo_cut)
  color.divisions <- 100
  my_colors <- colorRampPalette(c("navy", "white", "red"))(color.divisions)

  # breaks simmetrici attorno a 0
  # max_val <- max(max(abs(distmatrix, na.rm = TRUE)), min(abs(distmatrix, na.rm = TRUE)))
  max_val <- max(abs(distmatrix), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = color.divisions + 1)
  min_val <- -max_val
  # breaks <- seq(-max_val, max_val, length.out = color.divisions + 1)

  dendo_picture <- pheatmap(
    distinputmatrix,
    color = my_colors,
    cluster_rows = FALSE,
    cluster_cols = dendo,
    clustering_method = clus_method,
    scale = "none",
    # annotation_col = annotation_col,
    # breaks = seq(min_val, max_val, length.out = color.divisions + 1),
    breaks=breaks,
    cutree_cols = clusnum,
    show_rownames = TRUE,
    show_colnames = FALSE,
    border_color = "black",
    # annotation_colors = annotation_colors,
    annotation_legend = TRUE,
    legend = TRUE,
    na_col = "#DDDDDD",
    filename = outputname
  )
}


single_sample_analisys<-function(filename, clusnum){
  cat("start analysis on filename:", filename, "\n")
  database <- readRDS(filename)
  dendo_cut <- cutree(database$results$dendo, k = clusnum)
  sil<-silhouette(dendo_cut, dmatrix=database$results$distmatrix, do.clus.stat=TRUE, do.n.k=TRUE,do.col.sort=TRUE)
  summary(sil)
  avg_sil_width_all <- mean(sil[, "sil_width"])
  cat("Silhouette Score is:", avg_sil_width_all, "\n")
  sil_width_by_cluster <- tapply(sil[, "sil_width"], sil[, "cluster"], mean)
  if(isdebug){
    View(sil)
    View(database$results$distmatrix)
  }
  newfilename<-tools::file_path_sans_ext(filename)
  newfilename<-paste0(newfilename,"clusnum",clusnum,"_silhouette.pdf")
  pdf(file=newfilename)
  cols <- rainbow(clusnum) 
  plot(sil, col=cols)
  dev.off()

  cat(filename,"clusnum=",clusnum,"element_xcclus=",as.numeric(table(dendo_cut))," avg_sil_width_all=",avg_sil_width_all,"  sil_width_xclus=",sil_width_by_cluster,"\n",file="out_compdismatrix.txt",sep=" ",append=TRUE)

  return(sil)
}

# eval_medoids<-function(distmatrix){
#   class(distmatrix)
#   D <- as.matrix(distmatrix)
#   medoids <- sapply(unique(clusters), function(cl) {
#     idx <- which(clusters == cl)
#     idx[ which.min(rowSums(D[idx, idx])) ]

# medoids <- sapply(unique(dendo_cut), function(cl) {
#   samples_cl <- names(dendo_cut)[dendo_cut == cl]
#    subD <- database$results$distmatrix[samples_cl, samples_cl, drop = FALSE]
#    samples_cl[ which.min(rowSums(subD)) ]

# })
# }
