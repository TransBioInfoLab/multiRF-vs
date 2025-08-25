#' Multiple plot functions
#' @param dat A data frame or matrix
#' @param group Class group
#' @param position Legend position. The default is "bottomright"
#' @inheritParams graphics::plot
#' @return A plot object
#'
#' @export plot_tSNE

# plot tSNE
plot_tSNE <- function(dat, group = NULL, label_group = T, position = "right", perplexity = 30, pca = F, ncomp = 70, main = "tSNE", size = .5,...){


  t <- Rtsne::Rtsne(dat, perplexity = perplexity, pca = pca, initial_dims = ncomp, ...)
  embed <- t$Y
  
  if(is.null(group)) group <- rep("black", nrow(embed))

  df <- data.frame(tSNE1 = embed[,1], tSNE2 = embed[,2], group = group)
  df <- na.omit(df)
  p1 <- ggplot(df, aes(tSNE1, tSNE2)) +
    geom_point(aes(color = group), size = size) +
    xlab("tSNE dim 1") +
    ylab("tSNE dim 2") +
    ggtitle(main) +
    theme_classic()

  if(length(unique(group)) == 1) {
    p1 <- p1 + guides(color = "none")
  } else {
    p1 <- p1 +
      ggsci::scale_color_igv() +
      guides(color = guide_legend(override.aes = list(size=3))) +
      theme(legend.position = position,
            legend.title = element_blank(),
            legend.text = element_text(size = 6))
  }

  if(label_group & length(unique(group)) != 1) {
    data2 <- df %>% group_by(group) %>% dplyr::select(tSNE1, tSNE2) %>% summarize_all(mean)

    p1 <- p1 +
      ggrepel::geom_text_repel(data = data2, aes(label = group),
                               size = 2.5, color = "grey20",
                               max.overlaps = 15)
  }


  p1


}

#' @param dat A data frame or matrix
#' @param group Class group
#' @param label A logical parameter that determine whether to show labels or not.
#' @param cutoff Cutoff of edgeweights. The default is mean value of dat.
#' @param position Legend position. The default is "bottomright"
#' @inheritParams igraph::plot.igraph
#'
#' @export
#' @rdname plot_tSNE

# plot sample network
plot_network <- function(dat, group = NULL, label = F, cutoff = mean(dat),
                         layout = layout_with_fr, vertex.size = 5, label.dist =1,
                         edge.width = 0.3, vertex.frame.color = "white", vertex.label.color = "black",
                         vertex.label.cex = 0.5, vertex.label.dist = 0,
                         edge.curved = 0.5,vertex.label.degree =  -pi/2,
                         position = "bottomright", ...){

  diag(dat) <- 0
  network <- igraph::graph_from_adjacency_matrix(dat,
                                                 mode = "undirect",
                                                 weighted = T)

  if(!is.null(group)){
    n_class <- length(unique(group))
    get_p <- ggpubr::get_palette("jco", k = n_class)
    names(get_p) <- sort(unique(group))
    col <- get_p[group]
    V(network)$class <- col
  }

  if(!label){
    label <- NA
  } else {label <- V(network)$ids}

  network <- delete_edges(network, E(network)[which(E(network)$weight < cutoff)])
  network <- delete.vertices(network, igraph::degree(network) == 0)

  plot(network,
       vertex.color = V(network)$class,
       vertex.size = vertex.size,
       vertex.label = label,
       label.dist = label.dist,
       edge.width = edge.width,
       vertex.frame.color = vertex.frame.color,
       vertex.label.color = vertex.label.color,
       vertex.label.cex = vertex.label.cex,
       vertex.label.dist = vertex.label.dist,
       edge.curved = edge.curved,
       vertex.label.degree = -pi/2,
       layout = layout,
       ...)
  if(!is.null(group)){
    legend(position, fill = get_p, legend = unique(group))
  }

}


#' @param dat A data frame or matrix
#' @param group Class group
#' @param top The number of top weighted variables to show in the plot. The default is 20. Can be chosen from numerical values or all for setting the parameter = NULL.
#' @import ggplot2
#'
#' @export
#' @rdname plot_tSNE

plot_weights <- function(weights, plot.which = "all", top = 20, labels = NULL){

  if(toupper(plot.which) == "ALL"){
    imp <- weights
  } else {
    imp <- weights[plot.which]
  }

  df <-  purrr::map(imp,
                    ~data.frame(Weights = sort(., decreasing = T),
                                Var_names = names(.)[order(., decreasing = T)],
                                Ranks = c(1:length(.))))

  if(!is.null(top)){

    df <- purrr::map(df, ~.[1:top,])

    p <- plyr::llply(
      df,
      .fun = function(g){
        p <- ggplot(
          data = g,
          aes(x = reorder(Var_names, -order(Ranks, decreasing = F)), y = Weights, fill = Weights)
        ) +
          geom_bar(stat = 'identity') +
          ggpubr::theme_pubr(base_size = 9) +
          theme(axis.text.x = element_text(size = .5)) +
          xlab("Variables") +
          coord_flip() +
          scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(8, "Oranges"))(10))

      }
    )

  } else {

    p <- plyr::llply(
      df,
      .fun = function(g){
        p <- ggplot(
          data = g,
          aes(x = Var_names, y = Weights)
        ) +
          geom_point() +
          ggpubr::theme_pubr(base_size = 9) +
          theme(axis.text.x = element_blank()) +
          xlab("Variables")

      }
    )

  }
  if(is.null(labels)) labels <- names(df)
  ggpubr::ggarrange(plotlist = p, nrow = 1, labels = labels)

}


#' @param weights A list of weights get from mrf3
#' @param group Class group
#' @param order_by_group A logical parameter that determine whether the plot is ordered by group or by hierarchical clustering. The default is FALSE.
#' @param symm A logical parameter that determine whether the data is symmetric.
#' @inheritParams pheatmap::pheatmap
#'
#' @export
#' @rdname plot_tSNE

plot_heatmap <- function(dat, group = NULL, cluster_cols = T,
                         fontsize_row = 5, border_color = NA,
                         ...){

  if(!is.null(group)){

    annotation <- data.frame(Cluster = as.factor(group))
    rownames(annotation) <- rownames(dat)
    mat <- dat

  } else {
    mat <- dat
    annotation <- NA
  }

  pheatmap(t(mat), annotation_col = annotation,
           show_colnames = F,
           cluster_cols = cluster_cols,
           fontsize_row = fontsize_row,
           border_color = border_color,
           ...)
}

#' @export
#' @rdname plot_tSNE
plot_circos <- function(mat, names.list, group = NULL, cut.off = NULL, highlight = NULL, ...){

  vimp <- mat

  group <- names(names.list)

  if(!is.null(cut.off)){
    vimp[vimp <= cut.off] <- 0
    # vimp <- vimp[ rowSums(vimp) != 0 ,]
    # vimp <- vimp[, colnames(vimp) %in% rownames(vimp) ]
    # names.list <- lapply(names.list, function(l) l[l %in% colnames(vimp)])
  }

  # Circular Network Diagram Plot

  group.list <- lapply(1:length(group), function(g){
    structure(rep(group[g], length(names.list[[g]])), names =  names.list[[g]])
  }
  )

  col_mat <- vimp
  col_mat[col_mat > 1e-05] <- 2
  col_mat[col_mat == 1e-05] <- "#FFFFFF00"
  col_mat[col_mat < 1e-05] <- 3

  circos <- chordDiagram(vimp + 1e-05, annotationTrack = "grid",
                         group = unlist(group.list),
                         preAllocateTracks = list(
                           track.height = mm_h(4),
                           track.margin = c(mm_h(4), 0),
                           symmetric = T,
                           col = col_mat,
                           ...)
  )


  circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index.all = get.cell.meta.data("sector.index")
    for(i in c(1:length(names.list))){

      sector.index = sector.index.all[sector.index.all %in% names.list[[i]]]
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")

      circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.4 , facing = "clockwise", niceFacing = TRUE)

    }}, bg.border = NA)



  for(i in c(1:length(names.list))){
    highlight.sector(names.list[[i]], track.index = 1, col = i,
                     text = group[i], cex = 0.8, text.col = "white", niceFacing = TRUE)
  }

  if(!is.null(highlight)){
    highlight <- intersect(highlight,names.list$X)
    print(length(names.list$X))
    print(length(highlight))
    highlight.sector(highlight, col = "#00FF0040")
  }


  circos.clear()


}



#' @export
#' @rdname plot_tSNE
plot_embed <- function(dat, group = NULL, position = "bottom", pch = 20, ...){


  if(!is.null(group)){
    l <- sort(na.omit(unique(group)))
    get_p <- ggpubr::get_palette("jco", k = length(l))
    names(get_p) <- l
    col <- get_p[group]
  } else {col = 1}

  pairs(dat, col = col, pch = pch,  oma=c(10,4,6,4), ...)
  par(xpd=TRUE)
  if(!is.null(group)){
    legend(position, fill = get_p, legend = l, horiz=TRUE,
           xpd=TRUE, bty="n", cex = .5)
  }


}


#' @export
#' @rdname plot_tSNE
plot_umap <- function(dat, group = NULL, group2 = NULL, main = "UMAP", label_group = T,
                      position = "right", seed = 1080,...){

  set.seed(seed)
  t <- uwot::umap(dat, ...)

  if(is.null(group)) group <- rep("black", nrow(dat))
  
  df <- data.frame(umap1 = t[,1], umap2 = t[,2], group = group)
  if(!is.null(group2)) df$group2 <- group2
  
  df <- na.omit(df)

  p1 <- ggplot(df, aes(umap1, umap2)) +
    geom_point(aes(color = group), size = 1) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    ggtitle(main) +
    theme_classic()

  if(length(unique(group)) == 1) {
    p1 <- p1
  } else {
    p1 <- p1 +
      ggsci::scale_color_igv() +
      guides(color = guide_legend("PAN-Cancer Types", override.aes = list(size=3))) +
      theme(legend.position = position,
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            title = element_text(size = 12, face = "bold"))
  }
  
  if(!is.null(group2)) {
    p1 <- p1 + 
      geom_point(aes(shape = group2)) +
      scale_shape_manual("NMF Clusters", values=1:nlevels(df$group2)) 
  }
  
  if(label_group & length(unique(group)) != 1) {
    data2 <- df %>%
      group_by(group) %>%
      dplyr::select(umap1, umap2) %>%
      summarize_all(mean)

    p1 <- p1 +
      ggrepel::geom_text_repel(data = data2, aes(label = group),
                               size = 3.5, color = "red",
                               max.overlaps = 20)
  }


  p1



}

# KM plot
KM_plot <- function(test_var, time_var, event_var, pheno_mat, cut = "median", ...){
  
  if(is.numeric(test_var)){
    
    if(cut == "median"){
      m <- median(test_var)
    }
    
    if(cut == "mean"){
      m <- mean(test_var)
    }
    
    if(cut == "maxstat"){
      
      df <- data.frame(cluster = test_var, time = pheno_mat[[time_var]], death = pheno_mat[[event_var]])
      m <- maxstat::maxstat.test(Surv(time, death) ~ cluster, data = df, smethod = "LogRank")$estimate
      
    }
    gene_cut <- ifelse(test_var < m, "low", "high")
    
  } else {
    
    gene_cut = test_var
    
  }
  
  df <- data.frame(cluster = gene_cut, time = pheno_mat[[time_var]], death = pheno_mat[[event_var]])
  
  fo <- as.formula(paste0("Surv(time, death) ~ cluster"))
  fit <- surv_fit(fo, data = df)
  p <- round(-log10(survdiff(fo, data = df)$p),4)
  
  survminer::ggsurvplot(
    fit,
    data = df,
    size = 1,
    palette = "npg",
    conf.int = F,
    pval = T,
    risk.table = F,
    xlab = "Time (Days)",
    ...
  )
}
