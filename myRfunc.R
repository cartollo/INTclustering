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
  max_val <- max(distmatrix, na.rm = TRUE)
  min_val <- min(distmatrix, na.rm = TRUE)
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

single_sample_analisys<-function(dendo, clusnum, distmatrix, debug=FALSE){
  dendo_cut <- cutree(dendo, k = clusnum)
  sil<-silhouette(dendo_cut, dmatrix=distmatrix, do.clus.stat=TRUE, do.n.k=TRUE,do.col.sort=TRUE)
  summary(sil)
  avg_sil_width <- mean(sil[, "sil_width"])
  if(debug){ 
    cat("Silhouette Score is:", avg_sil_width, "\n")
    View(sil)
    View(corematrix$results$distmatrix)
    plot(sil)
  }
  return(c(sil, avg_sil_width))
}
