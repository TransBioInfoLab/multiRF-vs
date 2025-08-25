library(igraph)
#' Find pairwise IMD between variables across datasets
#' @param x A mrf3 object
#' @param all_var A logical parameter that determines whether to compute the connections of all variables or selected variables.
#' The default is FALSE for memory saving.
#' @return A pairwise adjacency matrix between variables and a pairwise connection matrix between datasets
#' @export pairwise_imd

pairwise_imd <- function(mod, all_var = F, normalized = F){
  
 
  net <- mod$net
  connection <- mod$connection
  dat_names <- names(mod$weights)
  
  nt <- mod$ntree
  
  var_names <- purrr::map(mod$weights, ~names(.))
  
  num_dat <- length(mod$weights)
  am <- matrix(0, nrow = num_dat, ncol = num_dat,
               dimnames = list(dat_names, dat_names))
  
  d <- do.call(rbind, connection) 
  length_count <- table(d)[dat_names]
  
  for (i in 1:nrow(d)) {
    am[d[i, 1], d[i, 2]] <- 1
    am[d[i, 2], d[i, 1]] <- 1
  }
  diag(am) <- length_count
  
  large_mat <- matrix(0, nrow = length(unlist(var_names)), ncol = length(unlist(var_names)),
                      dimnames = list(unlist(var_names), unlist(var_names)))
  
  l_ply(
    1:length(net),
    .fun = function(i){
      ne <- net[[i]]
      vn <- var_names[connection[[i]]]
      names(vn) <- c("yvar", "xvar")

      mat <- pairwise_imd_mod(ne, vn)

      mat <- mat[rownames(mat) %in% rownames(large_mat), colnames(mat) %in% colnames(large_mat)]
      large_mat[rownames(mat),colnames(mat)] <<- large_mat[rownames(mat),colnames(mat)] + mat
    }
  )
  
  for (i in 1:ncol(am)) {
    for(j in (i:ncol(am))){
      m1 <- colnames(am)[i]
      m2 <- colnames(am)[j]
      if(am[i,j] != 0){
        large_mat[var_names[[m1]],var_names[[m2]]] <- large_mat[var_names[[m1]],var_names[[m2]]]/am[i,j]
        if(i != j){
          large_mat[var_names[[m2]],var_names[[m1]]] <- large_mat[var_names[[m2]],var_names[[m1]]]/am[i,j]
        }
      }
    }
  }
  
  if(!all_var) {
    large_mat <- large_mat[names(Reduce(c,mod$weights))[Reduce(c,mod$weights) != 0],
                           names(Reduce(c,mod$weights))[Reduce(c,mod$weights) != 0]]
    var_names <- purrr::map(mod$weights, ~names(.)[. > 0])
    
  }

  large_mat <- large_mat/nt
  if(normalized) {
    d <- rowSums(large_mat)
    l <- diag(d^(-1/2)) %*% large_mat %*% diag(d^(-1/2))
    dimnames(l) <- dimnames(large_mat)
    large_mat <- l
  }
  
  
  out <- list(
    adj_var_mat = large_mat,
    adj_dat_mat = am,
    var_use = var_names
  )
  
  return(out)
  
}


pairwise_imd_mod <- function(net_list, var_name) {
  doParallel::registerDoParallel(20)
  mat_list <- plyr::llply(
    net_list,
    .fun = function(net) {
      mat <<- pairwise_imd_tree(net, var_name)
    }, .parallel = T
  )
  Reduce("+",mat_list)
}

pairwise_imd_tree <- function(net, var_name) {
  
  subtrees_list <- find_subtrees(net)
  
  # 4) init output matrix over your x/y vars
  all_vars <- unique(c(var_name$xvar, var_name$yvar))
  mat <- matrix(
    0,
    nrow = length(all_vars),
    ncol = length(all_vars),
    dimnames = list(all_vars, all_vars)
  )
  depth <- 1/net$inv_d
  names(depth) <- net$from

  # 5) IMD weighting function
  imd_w <- function(x,y) {
    dx <- depth[x]; dy <- depth[y]
    if (is.na(dx)||is.na(dy)) return(0)
    1 / (abs(dx - dy) * pmax(dx, dy))
  }
  
  # 6) for each first‐level subtree, accumulate all pairwise IMD
  plyr::l_ply(subtrees_list, .fun = function(sub) {
    root <- sub[1]
    root_node <- sub("^(.*)_.*", "\\1", root)  
    root_y <- na.omit(net$Y_id[which(net$from == root)])
    mat[root_node, root_node] <<- max(mat[root_node, root_node], 1/(depth[root]))
    mat[root_y, root_y] <<- max(mat[root_y, root_y], 1/(depth[root]))
    for (child in sub[-1]) {
      if(grepl("<leaf>", child)) next
      
      child_y <- na.omit(net$Y_id[which(net$from == child)])
      w <- imd_w(root, child)
      
      child_node <- sub("^(.*)_.*", "\\1", child)
      mat[root_node, child_node] <<- max(mat[root_node, child_node], w)
      mat[child_node, root_node] <<- max(mat[child_node, root_node], w)
      mat[root_y,child_y] <<- max(mat[root_y, child_y], w)
      mat[child_y,root_y] <<- max(mat[child_y, root_y], w)
      mat[root_node, child_y] <<- max(mat[root_node, child_y], w)
      mat[child_y, root_node] <<- max(mat[child_y, root_node], w)
      
    }
  })
  
  # mat[lower.tri(mat)] <- mat[upper.tri(mat)]
  
  mat
}


build_tree <- function(net) {
  # # strip off any "_…" suffix to get the true node IDs
  # net$parent <- sub("^(.*)_.*$", "\\1", net$from)
  # net$child  <- sub("^(.*)_.*$", "\\1", net$to)
  
  # unique edge list

  edges <- unique(net[, c("from","to")])
  graph_from_data_frame(edges, directed = TRUE)
}

get_subtree <- function(graph, root) {
  if (! root %in% V(graph)$name) {
    stop("'", root, "' is not a node in the tree")
  }
  subcomponent(graph, v = root, mode = "out")$name
}

find_subtrees <- function(net) {
  g <- build_tree(net)
  # internal nodes = those that appear as a parent
  # internal_nodes <- unique(sub("^(.*)_.*$", "\\1", net$from))
  nodes <- net$from
  # for each, get its subtree
  setNames(
    lapply(nodes, function(r) get_subtree(g, r)),
    nodes
  )
}
