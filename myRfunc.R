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
dendo_picture<-function(distinputmatrix, clus_method, dendo, kmeans_cluster = NULL,split_by_kmeans = FALSE){
  # table(dendo_cut)
  color.divisions <- 100
  my_colors <- colorRampPalette(c("navy", "white", "red"))(color.divisions)

  # breaks simmetrici attorno a 0
  #zscore esplicito anche se dovrebbe essere inutile
  distinputmatrix <- t(scale(t(distinputmatrix), center = TRUE, scale = TRUE))
  max_val <- max(abs(distinputmatrix), na.rm = TRUE)
  breaks <- seq(-max_val, max_val, length.out = color.divisions + 1)
  breaks <- unique(breaks)
  gaps_col <- NULL
  annotation_col <- NULL
  #se è kmeans:
  if (split_by_kmeans) {
    missing <- setdiff(colnames(distinputmatrix), names(kmeans_cluster))
    if (length(missing) > 0) {
      stop(paste0("Missing kmeans labels for samples: ", paste(missing, collapse = ", ")))
    }
    cl <- kmeans_cluster[colnames(distinputmatrix)]
    cl <- as.integer(as.factor(cl)) #voglio che i cluster siano numerici
    annotation_col <- data.frame(kmeans = factor(cl))
    rownames(annotation_col) <- colnames(distinputmatrix)
    #ordino i cluster in modo che rimanghino allineati
    ord <- order(cl)
    distinputmatrix <- distinputmatrix[, ord, drop = FALSE]
    annotation_col <- annotation_col[ord, , drop = FALSE]
    cl_ord <- cl[ord]

    tab <- table(cl_ord)
    gaps_col <- cumsum(as.integer(tab))
    gaps_col <- gaps_col[-length(gaps_col)]

  cluster_cols_arg<-FALSE
  }else{
    cluster_cols_arg<-dendo
  }

  mydendo_picture <- pheatmap(
    distinputmatrix,
    color = my_colors,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols_arg,
    clustering_method = clus_method,
    scale = "none",
    # annotation_col = annotation_col,
    # breaks = seq(min_val, max_val, length.out = color.divisions + 1),
    breaks=breaks,
    gaps_col = gaps_col,         
    annotation_col = annotation_col, 
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
single_sample_analisys<-function(filename, clusnum, noprint, txtoutfilename, newfile, folder=".", iscore, iskmeans=TRUE, plot_result=TRUE){
  cat("start analysis on filename:", filename, "\n")
  #carico database
  database <- readRDS(filename)
  dendo_cut <- cutree(database$results$dendo, k = clusnum)
  distmatrix_mat=as.matrix(database$results$distmatrix)
  sil<-silhouette(dendo_cut, dmatrix=distmatrix_mat, do.clus.stat=TRUE, do.n.k=TRUE,do.col.sort=TRUE)
  filename_nopath<-tools::file_path_sans_ext(filename)
  database$config$databasename<-filename_nopath
  if(noprint){
    return(list(sil=sil, dendo_cut=dendo_cut,distmatrix=database$results$distmatrix, filename=filename_nopath))
  }
  #apro file in uscita
  if(newfile){
    txtoutfilename<-paste0(folder,"/",filename_nopath,"_analysis.txt")
    filename_nopath<-paste0(folder,"/",filename_nopath,"clusnum_",clusnum,"_analysis.pdf")
    pdf(file=filename_nopath)
  }
  grid::grid.newpage()
  print_config_on_pdf(database$config)
  # write_report_header(
  # title = filename,
  # lines = database$config
  # )

  #disegno la pheatmap
  if(database$config$ispca==FALSE && database$config$whitemethod=="NOWHITE"){
    if(iscore){
      dendofile <- readRDS("NOWHITE_isclr_isrelab_iscore.rds")
      dendomap<-dendofile$results$distinputmatrix
    }else{
      dendofile <- readRDS("NOWHITE_isclr_isrelab_isfull.rds")
      dendomap<-dendofile$results$distinputmatrix
    }
  }else{
    dendomap<-database$results$distinputmatrix
  }

  mydendo_picture<-dendo_picture(distinputmatrix=dendomap, clus_method=database$config$clus_method, dendo=database$results$dendo)
  print(mydendo_picture)
  avg_sil_width_all <- mean(sil[, "sil_width"])
  cat("Silhouette Score is:", avg_sil_width_all, "\n")
  sil_width_by_cluster <- tapply(sil[, "sil_width"], sil[, "cluster"], mean)  
  cols <- rainbow(clusnum) 
  plot(sil, col=cols)
  cat(filename,"clusnum=",clusnum,"element_xcclus=",as.numeric(table(dendo_cut))," avg_sil_width_all=",avg_sil_width_all,"  sil_width_xclus=",sil_width_by_cluster,"\n",file=txtoutfilename,sep=" ",append=TRUE)

  #UMAP analysis 0<kmin<180, 0<min_dist<0.99
  umap_res<-umap_from_distmatrix(distmatrix_mat=distmatrix_mat,dendo_cut=dendo_cut, dendo=dendo, kvalue=10, min_dist=0.05, iskmeans=iskmeans, plot_result=plot_result, clusnum=clusnum, distinputmatrix=dendomap)
  umap_res<-umap_from_distmatrix(distmatrix_mat=distmatrix_mat,dendo_cut=dendo_cut,dendo=dendo, kvalue=15, min_dist = 0.1, iskmeans=iskmeans,plot_result=plot_result, clusnum=clusnum,distinputmatrix=dendomap)
  umap_res<-umap_from_distmatrix(distmatrix_mat=distmatrix_mat,dendo_cut=dendo_cut, dendo=dendo, kvalue=15, min_dist = 0.3, iskmeans=iskmeans, plot_result=plot_result, clusnum=clusnum,distinputmatrix=dendomap)
  umap_res<-umap_from_distmatrix(distmatrix_mat=distmatrix_mat,dendo_cut=dendo_cut, dendo=dendo, kvalue=30, min_dist = 0.1, iskmeans=iskmeans, plot_result=plot_result, clusnum=clusnum,distinputmatrix=dendomap)

  #tsne analysis perplexity<60
  tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = 5, plot_result=plot_result)
  tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = 10, plot_result=plot_result)
  tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = 15, plot_result=plot_result)
  tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = 30, plot_result=plot_result)

# if(database$config$clus_method=="euclidan" && iskmeans==TRUE){
#     kmeans_end<-kmeans_scan(t(dendofile$results$distinputmatrix), k_range = 2:10, seed = 123, plot_result=plot_result)
# }

  return(list(sil=sil, dendo_cut=dendo_cut,distmatrix=database$results$distmatrix, filename=filename, txtoutfilename=txtoutfilename))
}


compare_twoclustering<-function(fullres, coreres, txtoutfilename, folder="."){
  #adjusted mutual information 0=randomico, 1=perfetto
  ami <- AMI(fullres$dendo_cut, coreres$dendo_cut)
  #adjusted rand index 0=randomico, 1=perfetto
  ari <- ARI(fullres$dendo_cut, coreres$dendo_cut)
  
  #contingency table
  cont_tab <- table(
  full = fullres$dendo_cut,
  core  = coreres$dendo_cut
)

#scrivo cose su pdf
grid::grid.newpage()
grid.text(
  paste("Clustering comparison between:\n",fullres$filename,"\n", coreres$filename),
  x = 0.05, y = 0.95,
  just = c("left", "top"),
  gp = gpar(fontsize = 10)
)
grid.text(
  sprintf("adjusted rand index 0=randomico, 1=perfetto\n ARI = %.3f\n adjusted mutual information 0=randomico, 1=perfetto \nAMI = %.3f", ari, ami),
  x = 0.05, y = 0.80,
  just = c("left", "top"),
  gp = gpar(fontsize = 10)
)
tab_lines <- capture.output(print(cont_tab))
grid.text(
  paste("contingency table\n",paste(tab_lines, collapse = "\n")),
  x = 0.05, y = 0.60,
  just = c("left", "top"),
  gp = gpar(fontsize = 10)
)

pheatmap(
  prop.table(cont_tab, 1),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Contingency table pheatmap"
)

#print results on txt file
outfilename=paste0(folder,"/",txtoutfilename)
cat("AMI=",ami,"ARI=",ari,"\n","\n",sep=" ",file=outfilename,append=TRUE)

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

########################################################
applypca<-function(distinputmatrix, threshold=0.5){
  distinputpca <- princomp(x=t(distinputmatrix), scores=TRUE)
  var_pc <- distinputpca$sdev^2
  prop_var <- var_pc / sum(var_pc) #varianza di ogni nuova variabile
  cum_prop_var <- cumsum(prop_var) #cumulativa di varianza
  k <- which(prop_var >= threshold)
  distinputmatrix <- t(distinputpca$scores[, 1:length(k), drop = FALSE])
  return (distinputmatrix)
}

#########################################################
# Funzione wrapper per UMAP su distance matrix
umap_from_distmatrix <- function(distmatrix_mat, dendo_cut, dendo,
                                 kvalue, min_dist, seed = 123,
                                 plot_result = TRUE, iskmeans=FALSE, clusnum, mydendo_picture=NULL,distinputmatrix=NULL) {
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
  umap_res <- umap( distmatrix_mat, nn_method=nn_list, n_neighbors=kvalue, min_dist=min_dist)
  labels<-paste0("K=", kvalue, "  min_dist=",min_dist)
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
      labels = labels,
      pos = 4,       
      cex = 0.9,     # dimensione testo
      col = "black"
    )           
    
  }

  if(iskmeans){
    # kmeans_end<-kmeans_scan(umap_res, k_range = 2:10, seed = seed, plot_result = plot_result)
    #faccio k means su umap
    kmeans_end<-kmeans_single(umap_res, clusnum=clusnum, seed = seed, plot_result = plot_result, label=labels)  
  
    mydendo_picture<-dendo_picture(distinputmatrix=distinputmatrix, clus_method="euclidean", dendo=dendo, kmeans_cluster=kmeans_end$km$cluster, split_by_kmeans=TRUE)
    print(mydendo_picture)

  }


  # Restituisci la matrice di coordinate
  return(umap_res)
}

###### kmenas scan  ##############
#Y deve essere una matrice campioni x features 
kmeans_scan <- function(Y,
                        k_range = 3:10,
                        seed = 123,
                        nstart = 50,
                        iter.max = 200,
                        plot_result = FALSE) {
  
  # controlli di base su input
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (nrow(Y) < 5) stop("Need at least 5 samples (rows).")
  if (ncol(Y) < 2) stop("Need at least 2 dimensions (columns).")
  if (any(!is.finite(Y))) stop("Y contains NA/NaN/Inf.")
  if (any(k_range >= nrow(Y))) stop("All K must be < number of samples.")

  dY <-stats::dist(Y)

  km_list <- vector("list", length(k_range))
  names(km_list) <- as.character(k_range)

  wss <- rep(NA_real_, length(k_range))
  sil_avg <- rep(NA_real_, length(k_range))
  names(wss) <- names(sil_avg) <- as.character(k_range)

  for (idx in seq_along(k_range)) {
    k <- k_range[idx]
    set.seed(seed)
    km <- stats::kmeans(Y, centers = k, nstart = nstart, iter.max = iter.max)
    km_list[[idx]] <- km
    wss[idx] <- km$tot.withinss
    #silhouette stuff
    sil <- cluster::silhouette(km$cluster, dY)
    sil_avg[idx] <- mean(sil[, "sil_width"])
  }

  best_k <- as.integer(names(which.max(sil_avg)))

  metrics <- data.frame(
    clusnum = k_range,
    WSS = as.numeric(wss),
    Silhouette = as.numeric(sil_avg)
  )

  if (plot_result) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(1, 2))
    plot(k_range, wss, type = "b", xlab = "K", ylab = "WSS",
          main = paste0("kmeans_scan - Elbow"))
    if (!is.na(best_k)) abline(v = best_k, lty = 2)
    plot(k_range, sil_avg, type = "b", xlab = "K", ylab = "Avg silhouette",
          main = paste0("kmeans_scan - Silhouette"))
    if (!is.na(best_k)) abline(v = best_k, lty = 2)
  }

  kmeans_end<-list(kmeans = km_list,metrics = metrics,best_k = best_k)
  return (kmeans_end)
}

########### k means singolo  ###########
kmeans_single <- function(Y,
                          clusnum,
                          seed = 123,
                          nstart = 50,
                          iter.max = 200,
                          plot_result = TRUE,
                          label = NULL) {

  #controlli di base
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (nrow(Y) < 5) stop("Need at least 5 samples (rows).")
  if (ncol(Y) < 2) stop("Need at least 2 dimensions (columns).")
  if (any(!is.finite(Y))) stop("Y contains NA/NaN/Inf.")

  if (!is.numeric(clusnum) || length(clusnum) != 1) stop("clusnum must be a single numeric value.")
  if (clusnum < 2) stop("clusnum must be >= 2.")
  if (clusnum >= nrow(Y)) stop("clusnum must be < number of samples.")

  # distance (for silhouette)
  dY <- stats::dist(Y)

  # kmeans
  set.seed(seed)
  km <- stats::kmeans(Y, centers = clusnum, nstart = nstart, iter.max = iter.max)

  # metrics
  wss <- km$tot.withinss
  sil <- cluster::silhouette(km$cluster, dY)
  sil_avg <- mean(sil[, "sil_width"])

  metrics <- data.frame(
    clusnum = clusnum,
    WSS = as.numeric(wss),
    Silhouette = as.numeric(sil_avg)
  )

if (plot_result) {
  # grid::grid.newpage()
  plot(Y[, 1], Y[, 2],
        col = km$cluster, pch = 19,
        xlab = "Dim1", ylab = "Dim2",
        main = paste0("kmeans_single clusters (clusnum=", clusnum, ")"),
        # sub=label
        )
  legend("topleft",
          legend = sort(unique(km$cluster)),
          col = sort(unique(km$cluster)),
          pch = 19, bty = "n")
  # grid::grid.newpage()
  plot(sil,
        main = paste0("kmeans_single silhouette (avg=", signif(sil_avg, 3), ", clusnum=", clusnum, ")"),
        border = NA)
  }

  kmeans_end <- list(
    kmeans = km,
    metrics = metrics,
    silhouette = sil,
    cluster = km$cluster,
    centers = km$centers
  )

  return(kmeans_end)
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
  lines<-append(lines,paste("clusnum = ",clusnum,sep=" "))
  
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