#########################################################
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

#########################################################
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


#########################################################
#fare un dendogramma
dendo_picture<-function(dendo_cut, clusnum=NA, distmatrix, distinputmatrix, clus_method, dendo){
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
    na_col = "#DDDDDD"
  )
}


#########################################################
#analisi di un singolo database creato con create_dismatrix
single_sample_analisys<-function(filename, clusnum){
  cat("start analysis on filename:", filename, "\n")
  #carico database
  database <- readRDS(filename)
  dendo_cut <- cutree(database$results$dendo, k = clusnum)
  distmatrix_mat=as.matrix(database$results$distmatrix)
  filename_nopath<-tools::file_path_sans_ext(filename)
  database$config$databasename<-filename_nopath

  #apro file in uscita
  filename_nopath<-paste0(filename_nopath,"_analysis.pdf")
  pdf(file=filename_nopath)
  grid::grid.newpage()
  print_config_on_pdf(database$config)
  # write_report_header(
  # title = filename,
  # lines = database$config
  # )

  #disegno la pheatmap
  mydendo_picture<-dendo_picture(dendo_cut,clusnum, database$results$distmatrix, database$results$distinputmatrix, database$config$clus_method, database$results$dendo)
  print(mydendo_picture)
  
  sil<-silhouette(dendo_cut, dmatrix=distmatrix_mat, do.clus.stat=TRUE, do.n.k=TRUE,do.col.sort=TRUE)
  summary(sil)
  avg_sil_width_all <- mean(sil[, "sil_width"])
  cat("Silhouette Score is:", avg_sil_width_all, "\n")
  sil_width_by_cluster <- tapply(sil[, "sil_width"], sil[, "cluster"], mean)  
  cols <- rainbow(clusnum) 
  plot(sil, col=cols)
  cat(filename,"clusnum=",clusnum,"element_xcclus=",as.numeric(table(dendo_cut))," avg_sil_width_all=",avg_sil_width_all,"  sil_width_xclus=",sil_width_by_cluster,"\n",file="out_compdismatrix.txt",sep=" ",append=TRUE)

  #UMAP analysis 0<kmin<180, 0<min_dist<0.99
  umap_res<-umap_from_distmatrix(distmatrix_mat=distmatrix_mat,dendo_cut=dendo_cut, kvalue=10, min_dist=0.05)
  umap_res<-umap_from_distmatrix(distmatrix_mat=distmatrix_mat,dendo_cut=dendo_cut, kvalue=15, min_dist = 0.1)
  umap_res<-umap_from_distmatrix(distmatrix_mat=distmatrix_mat,dendo_cut=dendo_cut, kvalue=15, min_dist = 0.3)
  umap_res<-umap_from_distmatrix(distmatrix_mat=distmatrix_mat,dendo_cut=dendo_cut, kvalue=30, min_dist = 0.1)

  #tsne analysis perplexity<60
  tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = 5)
  tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = 10)
  tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = 15)
  tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = 30)

  on.exit(dev.off(), add = TRUE)

  return(c(sil, dendo_cut,database$results$distmatrix))
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


#########################################################
# Funzione wrapper per UMAP su distance matrix
umap_from_distmatrix <- function(distmatrix_mat, dendo_cut,
                                 kvalue, min_dist, seed = 123,
                                 plot_result = TRUE) {
  # Controlli a caso
  stopifnot(is.matrix(distmatrix_mat))
  stopifnot(nrow(distmatrix_mat) == ncol(distmatrix_mat))
  stopifnot(all(diag(distmatrix_mat) == 0))
  stopifnot(all(distmatrix_mat == t(distmatrix_mat)))
  stopifnot(length(dendo_cut) == nrow(distmatrix_mat))
  
  # Imposta seed per riproducibilità
  set.seed(seed)
  ### umap###
  nn_list <- uwot:::dist_nn(distmatrix_mat, k = kvalue)
  umap_res <- umap( distmatrix_mat, nn_method=nn_list, n_neighbors=kvalue, min_dist=0.1)

  # Plot
  if(plot_result) {
    plot(
      umap_res[,1], umap_res[,2],
      col = dendo_cut,
      pch = 19,
      xlab = "UMAP1",
      ylab = "UMAP2",
      main = "UMAP from distance matrix"
    )
    legend("topleft",
           legend = sort(unique(dendo_cut)),
           col = sort(unique(dendo_cut)),
           pch = 19)
    text(
      x = min(umap_res[,1]), 
      y = min(umap_res[,2]), 
      labels = paste0("K=", kvalue, "  min_dist=",min_dist),
      pos = 4,       
      cex = 0.9,     # dimensione testo
      col = "black"
    )           
    
  }

  # Restituisci la matrice di coordinate
  return(umap_res)
}

########################################################
#tsne from distmatrix
tsne_from_distmatrix <- function(distmatrix_mat, dendo_cut,
                                 perplexity = 15, seed = 123,
                                 plot_result = TRUE) {

  set.seed(seed)
  ###tsne###
  tsne_res <- Rtsne(as.dist(distmatrix_mat), is_distance = TRUE, perplexity = perplexity, verbose = TRUE)  
  #plot risultati
  if(plot_result) {
    plot(
      tsne_res$Y[,1], tsne_res$Y[,2],
      col = dendo_cut,
      pch = 19,
      xlab = "t-SNE1",
      ylab = "t-SNE2",
      main = "t-SNE from distance matrix"
    )
    text(
      x = min(tsne_res$Y[,1]), 
      y = min(tsne_res$Y[,2]), 
      labels = paste0("perplexity=", perplexity),
      pos = 4,       
      cex = 0.9,     # dimensione testo
      col = "black"
    )    
    legend(
      "topleft",
      legend = sort(unique(dendo_cut)),
      col = sort(unique(dendo_cut)),
      pch = 19
    )
  }
return(tsne_res)
}


#########################################################
#scrivere testo in un pdf in pagina nuova
# write_report_header <- function(lines, title = NULL) {
#   grid.newpage()
  
#   if (!is.null(title)) {
#     grid.text(
#       title,
#       gp = gpar(fontsize = 10)
#     )
#   }
  
#   grid.text(
#     lines,
#     just = "left",
#     gp = gpar(fontsize = 10)
#   )
# }

print_config_on_pdf <- function(config_list, x = 0.05, y = 0.95, fontsize = 10, lineheight = 0.03) {
  # config_list = tobesaved$config
  # x, y = coordinate di partenza nel PDF
  # fontsize = dimensione testo
  # lineheight = distanza verticale tra le righe
  
  # costruisci vettore di stringhe
  lines <- sapply(names(config_list), function(n) {
    val <- config_list[[n]]
    #toglie path intero da inputfile
     if(n == "inputfile") val <- basename(val)
    # converte TRUE/FALSE in stringa
    if(is.logical(val)) val <- ifelse(val, "TRUE", "FALSE")
    paste(n, "=", val)
  })
  
  # crea viewport per il testo
  vp <- viewport(
    x = x, y = y,
    width = 0.95, height = 0.95,
    just = c("left", "top")
  )
  
  pushViewport(vp)
  
  # coordina y per le righe
  nlines <- length(lines)
  # y_positions <- rev(seq(1, 1 - lineheight * (nlines-1), length.out = nlines))
  y_positions <- 1 - seq(0, by = lineheight, length.out = nlines)
  
  for(i in seq_along(lines)) {
    tg <- textGrob(
      lines[i],
      x = 0, y = y_positions[i],
      just = c("left", "top"),
      gp = gpar(fontsize = fontsize)
    )
    grid.draw(tg)
  }
  
  popViewport()
}