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
#fare un dendogramma
dendo_picture<-function(distinputmatrix, clus_method, dendo, other_cluster = NULL,split_by_other = FALSE, hasnoise= FALSE, clusnum=0, title=title){
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
  #se è kmeans o snn o altro metodo per clustering e quindi faccio solo disegno senza clustering:
  if (split_by_other) {
    missing <- setdiff(colnames(distinputmatrix), names(other_cluster))
    if (length(missing) > 0) {
      stop(paste0("Missing labels for samples: ", paste(missing, collapse = ", ")))
    }
    cl_raw <- other_cluster[colnames(distinputmatrix)]
    annotation_col <- data.frame(cluster = factor(cl_raw))
    rownames(annotation_col) <- colnames(distinputmatrix)
    
    #ordino i cluster in modo che rimanghino allineati
    ord <- order(cl_raw)
    distinputmatrix <- distinputmatrix[, ord, drop = FALSE]
    annotation_col <- annotation_col[ord, , drop = FALSE]
    cl_ord <- cl_raw[ord]

    tab <- table(cl_ord)
    gaps_col <- cumsum(as.integer(tab))
    gaps_col <- gaps_col[-length(gaps_col)]

    cluster_cols_arg<-FALSE
  }else{
    cluster_cols_arg<-dendo
  }

  mydendo_picture <- pheatmap(
    distinputmatrix,
    main=title,
    color = my_colors,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols_arg,
    clustering_method = clus_method,
    scale = "none",
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

  #dendo picture iniziale
  mydendo_picture<-dendo_picture(distinputmatrix=dendomap, clus_method=database$config$clus_method, dendo=database$results$dendo, clusnum=clusnum, title="Dendrogram clustering")
  print(mydendo_picture)
  avg_sil_width_all <- mean(sil[, "sil_width"])
  cat("Silhouette Score is:", avg_sil_width_all, "\n")
  sil_width_by_cluster <- tapply(sil[, "sil_width"], sil[, "cluster"], mean)  
  cols <- rainbow(clusnum) 
  plot(sil, col=cols)
  cat(filename,"clusnum=",clusnum,"element_xcclus=",as.numeric(table(dendo_cut))," avg_sil_width_all=",avg_sil_width_all,"  sil_width_xclus=",sil_width_by_cluster,"\n",file=txtoutfilename,sep=" ",append=TRUE)


  #clustering con snn
  
  K_values<- c(10,20, 25, 30, 35, 40,50)
  Eps_values<- seq(from = 2, to = 20, by = 2)
  snn_results_single<-list()
  for (k in K_values) {
    for (eps in Eps_values) {
      if (eps >= k) next
    append(snn_results_single,snn_single( k = k, eps = eps, minPts = 5, distinputmatrix = dendomap, distmatrix_mat=distmatrix_mat, plot_result = plot_result, dendo=database$results$dendo, minclusnum=2, maxclusnum=10, maxnoise=15))
    }
  }

  #printo risultati di snn
  if(plot_result){
    plot.new()
    righe <- c(
      "SNN analysis results:",
      for (i in snn_results_single) {
        paste("k:", i$params$k, "eps:", i$params$eps, "minPts:", i$params$minPts, "n_noise:", i$n_noise, "silhouette_avg:", round(i$sil_avg, 3), "isok=", i$isok, "n_clusters:", i$n_clusters, "clustersize:", paste(i$cluster_sizes, sep=" "))
        if(i$isok==TRUE){
          if(i$n_clusters==clusnum)
            compare_twoclustering(first = dendo_cut, second = i$cl, firstfilename="dendrogram clustering", secondfilename=paste0("SNN clustering (k=",snn_results_single$params$k,", eps=",i$params$eps,", minPts=",i$params$minPts,")"))
        }
      }
    )
    text(0.1, seq(0.9, 0.9 - 0.05*(length(righe)-1), by = -0.05),
     labels = righe,
     adj = 0)
  }

  #faccio loop kmeans se è euclidean
  if(iskmeans && database$config$distance=="euclidean")
    kmeans_end<-kmeans_loop_metrics(Y=t(dendomap), k_range = 2:10, plot_result=plot_result, dendo=database$results$dendo, distinputmatrix=dendomap, dendo_cut=dendo_cut)

  #UMAP analysis 0<kmin<180, 0<min_dist<0.99
k_vals  <- c(10, 15, 15, 30)
min_vals <- c(0.05, 0.1, 0.3, 0.1)

for (i in seq_along(k_vals)) {
  umap_res <- umap_from_distmatrix(
    distmatrix_mat = distmatrix_mat,
    dendo_cut = dendo_cut,
    dendo = dendo,
    kvalue = k_vals[i],
    min_dist = min_vals[i],
    iskmeans = iskmeans,
    plot_result = plot_result,
    clusnum = clusnum,
    distinputmatrix = dendomap,
    title=paste0("UMAP from distance matrix (k=", k_vals[i], ", min_dist=", min_vals[i], ")")
  )
}

  #tsne analysis perplexity<60
  perplexity_val  <- c(5, 10, 15, 30)
  for (i in seq_along(perplexity_val)) {
    tsne_res<-tsne_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=dendo_cut, perplexity = perplexity_val[i], plot_result=plot_result, title(paste0("t-SNE from distance matrix (perplexity=", perplexity_val[i], ")"))
    )
  }

  return(list(sil=sil, dendo_cut=dendo_cut,distmatrix=database$results$distmatrix, filename=filename, txtoutfilename=txtoutfilename))
}


compare_twoclustering<-function(first, second, firstfilename, secondfilename, txtoutfilename=NULL, folder="."){
  #adjusted mutual information 0=randomico, 1=perfetto
  ami <- AMI(first, second)
  #adjusted rand index 0=randomico, 1=perfetto
  ari <- ARI(first, second)

  #contingency table
  cont_tab <- table(
    first = first,
    second  = second
  )
  cont_mat <- as.matrix(cont_tab)
  perm <- solve_LSAP(cont_mat, maximum = TRUE)
  # riordina colonne
  cont_mat_align <- cont_mat[, perm]
  cont_mat_align

  #scrivo cose su pdf
  grid::grid.newpage()
  grid.text(
    paste("Clustering comparison between:\n",firstfilename,"\n", secondfilename),
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
  tab_lines <- capture.output(print(cont_mat_align))
  grid.text(
    paste("contingency table\n",paste(tab_lines, collapse = "\n")),
    x = 0.05, y = 0.60,
    just = c("left", "top"),
    gp = gpar(fontsize = 10, fontfamily = "mono")
  )

  pheatmap(
    prop.table(cont_mat_align, 1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Contingency table pheatmap"
  )

  #print results on txt file
  if(length(txtoutfilename)>0){
    outfilename=paste0(folder,"/",txtoutfilename)
    cat("AMI=",ami,"ARI=",ari,"\n","\n",sep=" ",file=outfilename,append=TRUE)
  }
}

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
                                 kvalue, min_dist,
                                 plot_result = TRUE, iskmeans=FALSE, clusnum,distinputmatrix=NULL, title="UMAP from distance matrix") {
  # Controlli a caso
  stopifnot(is.matrix(distmatrix_mat))
  stopifnot(nrow(distmatrix_mat) == ncol(distmatrix_mat))
  stopifnot(all(diag(distmatrix_mat) == 0))
  stopifnot(all(distmatrix_mat == t(distmatrix_mat)))
  stopifnot(length(dendo_cut) == nrow(distmatrix_mat))
  
  ### umap###
  nn_list <- uwot:::dist_nn(distmatrix_mat, k = kvalue)
  umap_res <- umap( distmatrix_mat, nn_method=nn_list, n_neighbors=kvalue, min_dist=min_dist)
  labels<-paste0("K=", kvalue, "  min_dist=",min_dist)
  # Plot
  if(plot_result) {
    plot(
      main=title,
      umap_res[,1], umap_res[,2],
      col = dendo_cut,
      pch = 19,
      xlab = "UMAP1",
      ylab = "UMAP2"
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
    #faccio k means su umap
    kmeans_end<-kmeans_single(umap_res, kclusnum=clusnum, plot_result = plot_result, label=labels, isumap=TRUE, distinputmatrix = distinputmatrix, dendo=dendo, dendo_cut=dendo_cut)  
    compare_twoclustering(first = dendo_cut, second = kmeans_end$km$cluster, firstfilename="dendrogram clustering", secondfilename=paste0("kmeans on UMAP (",labels,")"))
  }


  # Restituisci la matrice di coordinate
  return(umap_res)
}


#######################
#k means loop   

kmeans_loop_metrics <- function(Y,
                                k_range = 3:10,
                                plot_result = FALSE, dendo=NULL, distinputmatrix=NULL, dendo_cut=NULL) {

  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (any(k_range >= nrow(Y))) stop("All K must be < number of samples.")
  if (plot_result==TRUE && is.null(dendo)) stop("dendo must be provided if plot_result is TRUE.")
  if (plot_result==TRUE && is.null(distinputmatrix)) stop("distinputmatrix must be provided if plot_result is TRUE.")

  # salva risultati per ogni K
  res_list <- vector("list", length(k_range))
  names(res_list) <- as.character(k_range)
  #inizializzo vettori
  wss <- rep(NA_real_, length(k_range))
  sil_avg <- rep(NA_real_, length(k_range))
  names(wss) <- names(sil_avg) <- as.character(k_range)

  for (i in seq_along(k_range)) {
    k <- k_range[i]
    res <- kmeans_single(Y=Y, kclusnum=k, plot_result=plot_result, isumap=FALSE, distinputmatrix=distinputmatrix, dendo=dendo, dendo_cut=dendo_cut)
    res_list[[i]] <- res
    wss[i] <- res$metrics$WSS
    sil_avg[i] <- res$metrics$Silhouette
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
         main = "kmeans elbow (WSS)")
    abline(v = best_k, lty = 2)

    plot(k_range, sil_avg, type = "b", xlab = "K", ylab = "Avg silhouette",
         main = "kmeans silhouette vs K")
    abline(v = best_k, lty = 2)
  }

  list(
    results = res_list,   
    metrics = metrics,
    best_k = best_k
  )
}

######################
kmeans_single <- function(Y,
                          kclusnum,
                          nstart = 50,
                          iter.max = 200,
                          plot_result = TRUE,
                          isumap=FALSE,
                          label = NULL,
                          distinputmatrix,
                          dendo,
                          dendo_cut) {

  #controlli di base
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (nrow(Y) < 5) stop("Need at least 5 samples (rows).")
  if (ncol(Y) < 2) stop("Need at least 2 dimensions (columns).")
  if (any(!is.finite(Y))) stop("Y contains NA/NaN/Inf.")

  if (!is.numeric(kclusnum) || length(kclusnum) != 1) stop("kclusnum must be a single numeric value.")
  if (kclusnum < 2) stop("kclusnum must be >= 2.")
  if (kclusnum >= nrow(Y)) stop("kclusnum must be < number of samples.")

  # distance (for silhouette)
  dY <- stats::dist(Y)

  # kmeans
  km <- stats::kmeans(Y, centers = kclusnum, nstart = nstart, iter.max = iter.max)

  # metrics
  wss <- km$tot.withinss
  sil <- cluster::silhouette(km$cluster, dY)
  sil_avg <- mean(sil[, "sil_width"])

  metrics <- data.frame(
    clusnum = kclusnum,
    WSS = as.numeric(wss),
    Silhouette = as.numeric(sil_avg)
  )


if (plot_result) {
    mydendo_picture<-dendo_picture(distinputmatrix=distinputmatrix, clus_method="euclidean", dendo=dendo, other_cluster=km$cluster, split_by_other=TRUE, hasnoise=FALSE, title=paste0("kmeans with K=", kclusnum))
    print(mydendo_picture)  
    if(kclusnum==clusnum){
      compare_twoclustering(first = dendo_cut, second = km$cluster, firstfilename="dendrogram clustering", secondfilename="kmeans with same number of clusters")
    }
  if(isumap){
    plot(Y[, 1], Y[, 2],
          col = km$cluster, pch = 19,
          xlab = "Dim1", ylab = "Dim2",
          main = paste0("kmeans_single clusters (clusnum=", kclusnum, ")"),
          # sub=label
          )
    legend("topleft",
            legend = sort(unique(km$cluster)),
            col = sort(unique(km$cluster)),
            pch = 19, bty = "n")
  }
  # grid::grid.newpage()
  plot(sil,
        main = paste0("kmeans_single silhouette (avg=", signif(sil_avg, 3), ", clusnum=", kclusnum, ")"),
        border = NA)
  }

  kmeans_end <- list(
    kmeans = km,
    metrics = metrics,
    silhouette = sil
  )

  return(kmeans_end)
}

##################################################

snn_single <- function(distmatrix_mat,
                      k, eps, minPts,
                       distinputmatrix,         # matrice feature×campioni per dendo/pheatmap
                       plot_result = TRUE,
                       dendo=dendo,
                       minclusnum=2, 
                       maxclusnum=10,
                       maxnoise_percent=15 #percentuale massima di noise accettabile per disegnare risultati
                       ) {

  # --- clustering ---
  X_snn<-as.dist(distmatrix_mat)
  snn <- dbscan::sNNclust(X_snn, k = k, eps = eps, minPts = minPts)
  cl <- snn$cluster
  cl <- ifelse(cl == 0, 104, cl) #sostituisco 15 a noise per colore plot
  names(cl) <- colnames(distmatrix_mat) #assegno nomi   
  # check cluster e noise
  n_noise <- sum(cl == 104)
  n_clusters <- length(setdiff(unique(cl), 104)) # toglie elementi di noise che hanno cluster==104
  maxnoise_num=maxnoise_percent*length(cl)/100

  # silhouette su spazio originale (escludendo noise=104) 
  if (n_clusters >= minclusnum && n_clusters <= maxclusnum && n_noise <= maxnoise_num ) {
    sil <- NULL
    sil_avg <- NA_real_
    keep <- cl !=104
    d_keep<-as.dist(distmatrix_mat[keep, keep, drop = FALSE])
    sil <- cluster::silhouette(cl[keep], d_keep)
    sil_avg <- mean(sil[, "sil_width"])

  # --- plotting ---
    if(plot_result) {
      #dendo/pheatmap split
      distinputmatrix <- as.matrix(distinputmatrix)  # feature×campioni
      # check: colonne devono essere campioni
      missing2 <- setdiff(colnames(distinputmatrix), names(cl))
      if (length(missing2) > 0) {
        stop(paste0("distinputmatrix columns not found in cluster labels: ", paste(missing2, collapse = ", ")))
      }

      mydendo_picture<-dendo_picture(distinputmatrix = distinputmatrix, clus_method="euclidean", dendo=dendo, other_cluster=cl, split_by_other=TRUE, hasnoise=TRUE, title=paste0("SNN clustering (k=", k, ", eps=", eps, ", minPts=", minPts, ")"))

      umap_res <- umap_from_distmatrix(distmatrix_mat=distmatrix_mat, dendo_cut=cl, dendo=dendo, kvalue=10, min_dist=0.05, plot_result=TRUE, iskmeans=FALSE, distinputmatrix=distinputmatrix, title=paste0("SNN clustering (k=", k, ", eps=", eps, ", minPts=", minPts, ")"))

      #silhouette
      if(!is.null(sil)) {
        plot(sil, main = paste0("SNN silhouette on ORIGINAL space (avg=", signif(sil_avg, 3), ")"),   border = NA)
      }
    }

    tobereturned<-list(
      cluster = cl,
      n_clusters = n_clusters,
      n_noise = n_noise,
      sil = sil,
      sil_avg = sil_avg,
      params = list(k = k, eps = eps, minPts = minPts),
      isok=TRUE
    )
  }    
  tobereturned<-list(
    cluster = cl,
    n_clusters = n_clusters,
    cluster_sizes = table(cl),
    n_noise = n_noise,
    sil = NULL,
    sil_avg = NA_real_,
    params = list(k = k, eps = eps, minPts = minPts),
    isok=FALSE
  )
    return(tobereturned)
}

########################################################
#tsne from distmatrix
tsne_from_distmatrix <- function(distmatrix_mat, dendo_cut,
                                 perplexity = 15,
                                 plot_result = TRUE, title="t-SNE from distance matrix") {

  ###tsne###
  tsne_res <- Rtsne(as.dist(distmatrix_mat), is_distance = TRUE, perplexity = perplexity, verbose = TRUE)  
  #plot risultati
  if(plot_result) {
    plot(
      tsne_res$Y[,1], tsne_res$Y[,2],
      main=title,
      col = dendo_cut,
      pch = 19,
      xlab = "t-SNE1",
      ylab = "t-SNE2"
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



######################################## old not used anymore functions ########################################  


#########################
#Y deve essere una matrice campioni x features 
kmeans_scan <- function(Y,
                        k_range = 3:10,
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



#########################################################
#scrivere testo in un pdf in pagina nuova
write_report_header <- function(lines, title = NULL) {
  grid.newpage()
  
  if (!is.null(title)) {
    grid.text(
      title,
      gp = gpar(fontsize = 10)
    )
  }
  
  grid.text(
    lines,
    just = "left",
    gp = gpar(fontsize = 10)
  )
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
