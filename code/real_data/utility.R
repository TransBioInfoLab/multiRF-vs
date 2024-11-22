library(msigdbr)
library(clusterProfiler)

ORA <- function(gene_list, pathway = "KEGG", method = "ORA", ...){
  
  if (pathway == "KEGG"){
    c2.cp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="KEGG")
    pathway.df <- c2.cp %>% dplyr::select(gs_name, human_gene_symbol)
  } else if (pathway == "REACTOME") {
    c2.cp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory ="REACTOME")
    pathway.df <- c2.cp %>% dplyr::select(gs_name, human_gene_symbol)
  } else if (pathway == "c5.bp"){
    c5.bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
    pathway.df <- c5.bp %>% dplyr::select(gs_name, human_gene_symbol)
  } else if (pathway == "c7.IMMUNESIGDB") {
    c7.IMMUNESIGDB <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB")
    pathway.df <- c7.IMMUNESIGDB %>% dplyr::select(gs_name, human_gene_symbol)
  } else if (pathway == "H"){
    h <- msigdbr(species = "Homo sapiens", category = "H")
    pathway.df <- h %>% dplyr::select(gs_name, human_gene_symbol)
  } else if (pathway == "c2.cp") {
    p1 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
    p2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
    p3 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = toupper("BioCarta"))
    p4 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "PID")
    p5 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = toupper("WikiPathways"))
    pathway.df <- rbind(p1,p2,p3,p4,p5)  %>% dplyr::select(gs_name, human_gene_symbol)
  } else {
    stop("Please specify the right pathway list or enter your own pathway list")
  }
  
  if(method == "ORA") {
    res <- enricher(
      gene = gene_list,
      TERM2GENE = pathway.df,
      ...
    )
  }
 
  if(method == "GSEA") {
    res <- GSEA(
      geneList = gene_list,
      TERM2GENE = pathway.df,
      ...
    )
  }
  
  return(res)
  
}

umap_fn <- function(dat.list, ...) {
  
  plyr::llply(
    dat.list,
    .fun = function(l) {
      umap::umap(l)$layout
    }
  )
 
}

plot_pathway <- function(results, topn = 10, ylim = 1e-5, padj = 0.05) {
  
  results <- results %>%
    group_by(category) %>% 
    slice_min(pvalue, n = topn) %>%
    filter(p.adjust < padj) %>%
    mutate(logp = -log10(pvalue),
           Description = gsub("_", " ", Description),
           Description = gsub("GOBP |REACTOME ", "", Description))
  
  scaleFUN <- function(x) sprintf("%.0f", x)
  results <- results %>%
    mutate(ordering = -as.numeric(as.factor(category)) + pvalue,
           Description = fct_reorder(Description, ordering, .desc = T))
  cbf_1 <- c( "#E69F00", "#56B4E9", "#009E73", 
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  ggplot(results, aes(x = Description, y = logp)) + 
    geom_col(aes(fill = category), orientation = "x") +
    scale_fill_discrete(type = cbf_1[1:length(unique(results$category))]) + 
    # guides(fill = "none") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0), limits = c(0,-log10(ylim)), labels=scaleFUN) +
    facet_grid(rows = vars(category), scales = "free_y", switch = "x", space = "free_y") + 
    theme_classic() +
    ylab(expression("-log"["10"]~"(P-value)")) + 
    xlab("") + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey20" ) + 
    theme(
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
      plot.title = element_text(size = 15, face = "bold"),
      strip.text.y = element_blank(),
      strip.placement = "outside",
      axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      legend.title = element_blank()
    )
  
  
  
}

plot_umap <- function(dat, group = NULL, main = "UMAP", label_group = T,
                      pca = T, ncomp = 70, position = "right", pch = 20, config = umap::umap.defaults,  method = "umap-learn",...){
  
  if(pca) {
    x <- prcomp(dat)$x[,1:ncomp]
  } else {
    x <- dat
  }
  
  t <- umap::umap(x, config = config, method = method)
  
  if(is.null(group)) group <- rep("black", nrow(dat))
  
  df <- data.frame(umap1 = t$layout[,1], umap2 = t$layout[,2], group = group)
  df <- na.omit(df)
  
  p1 <- ggplot(df, aes(umap1, umap2)) +
    geom_point(aes(color = group), size = .5) +
    xlab("UMAP1") +
    ylab("UMAP2") +
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
    data2 <- df %>%
      group_by(group) %>%
      dplyr::select(umap1, umap2) %>%
      summarize_all(mean)
    
    p1 <- p1 +
      ggrepel::geom_text_repel(data = data2, aes(label = group),
                               size = 2.5, color = "grey20",
                               max.overlaps = 15)
  }
  
  
  p1
  
  
  
}
