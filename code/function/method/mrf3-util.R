## Utility functions for mrf3.R

# Get tree net from random forest
get_tree_net <- function(mod, tree.id){

  xvar.names <- mod$xvar.names
  xvar.factor <- mod$xvar.factor
  native.array <- mod$forest$nativeArray
  native.f.array <- mod$forest$nativeFactorArray

  node.stat <- mod$node.stats
  native.array <- cbind(native.array, node.stat) %>% data.frame


  tree.df <- native.array %>% dplyr::filter(treeID == tree.id)

  #converted.tree <- display.tree
  vars.id <- data.frame(var = c("<leaf>", xvar.names), parmID = 0:length(xvar.names), stringsAsFactors = FALSE)
  tree.df$var <- vars.id$var[match(tree.df$parmID, vars.id$parmID)]

  num.node <- tree.df$nodeID %>% unique %>% length

  var.count <- 1:nrow(tree.df)
  lapply(unique(tree.df$var), function(v) {
    pt <- tree.df$var == v
    var.count[which(pt)] <<- 1:sum(pt)
  })

  tree.df$var_count <- var.count
  tree.df$var_conc <- paste0(tree.df$var, "_", tree.df$var_count)

  var.id <- c((num.node+1): nrow(tree.df))
  tree.df$var.tip.id <- 1:nrow(tree.df)
  tree.df$var.tip.id[tree.df$var != "<leaf>"] <- var.id
  tree.df$var.tip.id[tree.df$var == "<leaf>"] <- tree.df$var_count[tree.df$var == "<leaf>"]

  from_node <- ""
  network <- data.frame()
  num.children <- data.frame(tree.df, children = 0)
  num.children <- num.children[num.children$var != "<leaf>",, drop = FALSE]
  num.children <- num.children[!duplicated(num.children$var_conc),, drop = FALSE]
  num_children <- as.list(rep(0, nrow(num.children)))
  names(num_children) <- num.children$var_conc

  lapply(1:nrow(tree.df), function(i) {
    rowi <- tree.df[i, ]
    if (i == 1){
      from_node <<- rowi$var_conc
      from_id <<- rowi$var.tip.id
      ns_all <<- rowi$nodeSZ
      dpthST <- 1/(rowi$dpthST)
    } else{
      to_node <- rowi$var_conc
      to_id <- rowi$var.tip.id
      ns_part <- rowi$nodeSZ
      dpthST <- 1/(rowi$dpthST)
      new_node <- list(from = from_node, to = to_node,
                       from_id = from_id, to_id = to_id,
                       inv_d = dpthST,
                       edge = ns_all/ns_part,
                       nodesize = ns_part)
      network <<- data.frame(rbind(network, new_node, stringsAsFactors = FALSE))
      num_children[[from_node]] <<- num_children[[from_node]] + 1
      if(rowi$var != "<leaf>"){
        from_node <<- to_node
        from_id <<- to_id
        ns_all <<- ns_part
      }else{
        if(i != nrow(tree.df)){
          while(num_children[[from_node]] == 2){
            from_node <<- network$from[network$to == from_node]
            from_id <<- network$from_id[network$to_id == from_id]
            ns_all <<- ns_part
          }
        }
      }
    }
  })


  return(network)

}

# Get leaf node information
get_leaf_ds <- function(mod, tree.membership, net){
  
  # Find upper level node
  prev_leaf <- dplyr::filter(.data = net, is_leaf == 1)
  
  # Update mem_id and membership
  net$mem_id_old <- net$mem_id
  mem_update <- slice_max(.data = prev_leaf, order_by = mem_id, by = from, with_ties = F)
  mem_old <- slice_min(.data = prev_leaf, order_by = mem_id, by = from, with_ties = F)
  mem_new <- setNames(mem_update$mem_id, mem_update$from)
  mem_org_id <- mem_old$mem_id
  names(mem_org_id) <- mem_new
  sapply(1:nrow(mem_old), function(id) net$mem_id[net$mem_id == mem_org_id[id]] <<- as.numeric(names(mem_org_id)[id]))
  
  find_node <- unique(prev_leaf$from)
  
  leaf_select <- filter(net, to %in% find_node)
  map_id <- unique(leaf_select$from)
  
  leaf_select <- filter(net, from %in% map_id)
  leaf_select$mem_id[match(names(mem_new), leaf_select$to)] <- mem_new
  drop_leaf <- leaf_select$from[leaf_select$mem_id == 0]
  leaf_select <- leaf_select[!leaf_select$from %in% drop_leaf,]
  
  # Update leaf information
  
  net$is_leaf[net$to %in% prev_leaf$to] <- 0
  net$is_leaf[net$from %in% leaf_select$from] <- 1
  
  net$mem_id[match(names(mem_new), net$to)] <- mem_new
  
  tree_mem <- tree.membership
  
  sapply(1:length(mem_org_id),
         function(id){
           tree_mem[tree_mem == mem_org_id[id]] <<-
             as.numeric(names(mem_org_id)[id])
         })
  
  
  if(any(leaf_select$mem_id %in% 0)) {
    from_id <- leaf_select[leaf_select$mem_id %in% 0,"from"]
    leaf_select[leaf_select$from %in% from_id,"is_leaf"] <- 0
    leaf_drop <- leaf_select[leaf_select$is_leaf == 0,]
    net[net$to %in% leaf_drop$to,"is_leaf"] <- 0
    net[net$to %in% leaf_drop$to,"leaf_stack"] <- 1
  }
  
  return(
    list(net = net,
         tree.mem = tree_mem)
  )
}

# Get response outcome splitting scores for splits in a tree
get_Y_imp <- function(net, tree.membership, dat, w = NULL, yprob = 1, seed = -5){

  node_id <- unique(net$from)
  var_imp <- llply(
    node_id,
    .fun = function(id){
      mem_id <- filter(.data = net, from %in% id)
      mem_id <- mem_id$mem_id
      mem_selected <- tree.membership[tree.membership %in% mem_id]
      dat0 <- scale(dat[tree.membership %in% mem_id,])
      split_stat <- llply(
        mem_id,
        .fun = function(m_id){
          if(length(nrow(dat0[mem_selected %in% m_id,])) > 0) {
            dd <- dat0[mem_selected %in% m_id,]
            colSums(dd)^2/nrow(dd)
          } else {
            dat0[mem_selected %in% m_id,]^2
          }
          
          
        })

      split_stat <- Reduce("+", split_stat)
      samp <- rep(1, length(split_stat))
      if(!is.null(w)) {
        ns <- min(ceiling(length(split_stat)*yprob), length(w[w != 0]))
        w <- w/sum(w)
        set.seed(seed)
        samp0 <- sample.int(length(split_stat), ns, prob = w)
        samp0 <- (1:length(split_stat))[-samp0]
        samp[samp0] <- 0
      } else {
        set.seed(seed)
        samp0 <- sample.int(length(split_stat), ceiling(length(split_stat)/3))
        samp[samp0] <- 0
      }
     
      idx <- which.max(split_stat * samp)
      varY <- mean(split_stat)
      names(varY) <- colnames(dat)[idx]

      return(varY)
      }, .parallel = F)

  var_imp <- unlist(var_imp)

  return(var_imp)

}

# Get importance for leaves in a tree
get_tree_imp <- function(mod, dat = NULL, tree.membership, net, calc = "Both", w = NULL, yprob = 1, weighted = F, seed = -5){

  if(is.null(dat)){
    dat <- mod$xvar
  }

  if(calc %in% c("Y", "Both")){
    datY <- mod$yvar[rownames(dat),]
  }

  updated_net <- net

  top_node_info <- filter(.data = updated_net, is_leaf == 1)

  node_ds <- setNames(top_node_info$inv_d, top_node_info$from)

  old_net <- updated_net[updated_net$from %in% top_node_info$from,]
  scores_imp <- NULL

  old_net_corr <- group_by(.data = old_net, from)
  old_net_corr <- summarise_at(old_net_corr, .var = "inv_d", .funs = mean)

  node_ds <- node_ds[old_net_corr$from]

  use_sample <- old_net_corr$from

  match_old_net <- old_net[old_net$from %in% use_sample,]

  # Get important variable
  # Get X
  if(calc %in% c("Both", "X")){
    imp_var <- top_node_info[match(use_sample, top_node_info$from),"from"]
    scores_imp <- (old_net_corr$inv_d)
    if(weighted) {
      scores_imp <- scores_imp * old_net_corr$edge
    }
    # scores_imp <- drop_case[match(imp_var,drop_case$to),] %>% pull(corr)
    names(scores_imp) <- gsub("^(.*)_.*", "\\1",imp_var)
  }

  if(calc %in% c("Both", "Y")){
    # Get Y
    impY <- get_Y_imp(net = match_old_net, tree.membership = tree.membership, dat = datY, w = w, yprob = yprob, seed = seed)
    updated_net$Y_id[unique(match(match_old_net$from, updated_net$from))] <- names(impY)
    scores_impY <- updated_net$inv_d[match(unique(match_old_net$from), updated_net$from)]
    if(weighted) {
      scores_impY <- scores_impY * updated_net$edge[match(unique(match_old_net$from), updated_net$from)]
    }
    names(scores_impY) <- names(impY)

  }


  if(calc == "Both"){
    scores_imp <- list(X = scores_imp, Y = scores_impY)
  } else if(calc == "X"){
    scores_imp <- list(X = scores_imp)
  } else {
    scores_imp <- list(Y = scores_impY)
  }



  return(
    list(
      net = updated_net,
      dat = dat,
      imp_var = scores_imp
    )
  )
}

# Update tree importance from bottom to top
update_iter_imp <- function(mod, tree.id, calc = "Both", lambda, w = NULL, yprob = 1, weighted = F, seed = -5) {

  # Get the tree structure for the specified tree.id
  net <- get_tree_net(mod = mod, tree.id = tree.id)

  # Get the membership information for the specified tree.id
  mem <- mod$membership[, tree.id]

  # Initialize an empty data frame 'dat'
  dat <- NULL

  # Initialize an index 'i' and k
  i <- 1
  k <- 1

  # Initialize an empty list 'imp'
  imp <- list()

  net$Y_id <- NA

  # Iterate until all nodes in the tree are leaves or there are only 2 or fewer nodes left
  while (any(is.null(net$is_leaf), sum(net$used_leaf == 0) > 0, is.null(net$used_leaf))) {

    if(k > 1){
      if(is.null(net$used_leaf)) {
        net$used_leaf <- ifelse(net$is_leaf == 1, 1, 0)
      } else{
        net$used_leaf[net$is_leaf == 1] <- 1
      }
      # Get initial leaf information
      update_net <- get_leaf_ds(mod = mod,
                        tree.membership = mem,
                        net = net)
      mem <- update_net$tree.mem
      net <- update_net$net
    } else {
      # Create new columns
      if(is.null(net$is_leaf)) {
        net$is_leaf <- ifelse(grepl("leaf" ,net$to), 1, 0)
      }
      if(is.null(net$mem_id)) {
        net$mem_id <- ifelse(net$is_leaf == 1, net$to_id, 0)
      }
      find_node <- net[net$is_leaf == 1,]
      find_node <- dplyr::summarise(.data = find_node, n = n(), .by = from)
      find_node <- filter(.data = find_node, n == 1)
      net[net$is_leaf == 1,"is_leaf"][net[net$is_leaf == 1,]$from %in% find_node$from] <- 0
    }

    # Get the updated tree importance statistics
    update_ls <- get_tree_imp(
      mod = mod,
      dat = dat,
      tree.membership = mem,
      net = net,
      calc = calc,
      w = w,
      yprob = yprob,
      weighted = weighted,
      seed = seed
    )

    # Update 'net', 'mem', and 'dat' with the results from 'update_ls'
    net <- update_ls$net
    dat <- update_ls$dat

    # If importance for a variable is available, add it to the 'imp' list
    if (!is.null(update_ls$imp_var)) {
      imp[[i]] <- update_ls$imp_var
      # if(i == 1) {
      #   varls <- update_ls$imp_var
      #   varls <- purrr::map(varls, ~.-.)
      #   imp[[i]] <- varls
      # }
      i <- i + 1
    }

    k <- k + 1

  }


  # Filter out the root node(s)
  root_net <- filter(net, is_leaf == 1)

  # Initialize variables for variable importance for X and Y
  # if (calc %in% c("X", "Both")) {
  #
  #   imp_root <- 1/2
  #   names(imp_root) <- gsub("^(.*)_.*", "\\1", unique(root_net$from))
    impX <- unlist(purrr::map(imp, "X"))
  #   impX <- c(impX, imp_root)
  # }
  #
  if (calc %in% c("Y", "Both")) {
    impY <- unlist( purrr::map(imp, "Y") )
  }

  # Create lists of variable importance and variable names
  if (calc == "Both") {
    imp_ls <- list(impX, impY)
    var_name <- list(colnames(mod$xvar), colnames(mod$yvar))
  } else if (calc == "X") {
    imp_ls <- list(impX)
    var_name <- list(colnames(mod$xvar))
  } else {
    imp_ls <- list(impY)
    var_name <- list(colnames(mod$yvar))
  }

  # Calculate variable importance scores for each variable
  imp_col <- plyr::llply(
    1:length(imp_ls),
    .fun = function(l) {
      get_iv(var_name[[l]], imp_ls[[l]])
    }
  )

  # Set appropriate names for the variable importance scores
  if (calc == "Both") {
    names(imp_col) <- c("X", "Y")
  } else if (calc == "X") {
    names(imp_col) <- "X"
  } else {
    names(imp_col) <- "Y"
  }

  # add penalty
  x_freq <- mod$var.used[tree.id,]
  x_freq <- x_freq[x_freq != 0]
  imp_col <- add_lambda(imp_col, net, x_freq, lambda)

  out <- list(imp_ls = imp_col, net = net)
  return(out)
}

# Add penalty to weights
add_lambda <- function(imp_ls, net, x_freq, lambda){

  if(all(is.na(net$Y_id))){
    freq <- list(X = x_freq)
  } else {
    freq <- list(X = x_freq, Y = table(net$Y_id))
  }

  new_ww <- plyr::llply(
    names(imp_ls),
    .fun = function(i){

      var <- freq[[i]]
      var_imp_all <- imp_ls[[i]]
      var_imp <- var_imp_all[var_imp_all != 0]
      var <- var[names(var_imp)]

      if(any(var == 1)){
        var_imp[var == 1] <- var_imp[var == 1] * lambda
        var_imp_all[names(var_imp)] <- var_imp
      }

      var_imp_all
    }
  )

  names(new_ww) <- names(imp_ls)

  new_ww

}

#' Get forest importance
#' @rdname get_imp_forest
#' @export

get_imp_forest <- function(mod, parallel = T, calc = "Both", weighted = F, use_depth = F, normalized = F, lambda = 1, w = NULL, yprob = 1, cores = detectCores() - 2, seed = -5){

  nt <- mod$ntree

  if(parallel){
    if(Sys.info()["sysname"] == "Windows"){
      cluster <- parallel::makeCluster(cores)
    } else {cluster <- cores}
    doParallel::registerDoParallel(cluster)
  }

    if(calc == "Both") {
      cc <- "Both"
      idx <- c("X", "Y")
    }
    if(calc == "X") {cc <- "X";idx <- c("X")}
    if(calc == "Y") {cc <- "Y";idx <- c("Y")}

    results <- plyr::llply(1:nt,
                       .fun = function(t){
                         update_iter_imp(mod,
                                         tree.id = t,
                                         calc = cc,
                                         lambda = lambda,
                                         yprob = yprob,
                                         w = w,
                                         weighted = weighted,
                                         seed = seed)
                       }, .parallel = parallel)

    imp <- purrr::map(results, "imp_ls")
    net <- purrr::map(results, "net")

    imp_ls <- plyr::llply(idx,
                          .fun = function(im){
                            im_df <- Reduce(cbind, purrr::map(imp, im))
                            if(!is.null(ncol(im_df))){
                              if(use_depth) {
                                rs <- rowSums(im_df != 0)
                                rs[rs == 0] <- -999
                                imp <- rowSums(im_df)/rs
                              } else {
                                imp <- rowMeans(im_df)
                              }
                            } else imp <- im_df
                            if(normalized){
                              imp <- imp/sqrt(sum(imp^2))
                            }
                            return(imp)
                          })
    names(imp_ls) <- idx


  out <- list(
    imp_ls = imp_ls,
    imp_ls_init = imp,
    net = net
  )
  return(out)
}

# Match importance name to column name of the data frame
get_iv <- function(var_name, imp){

  imp_df <- data.frame(var_name = names(imp), importance = imp)
  imp_df <- dplyr::group_by(.data = imp_df,var_name)
  imp_df <-  dplyr::summarise_at(imp_df,"importance", max)

  var_v <- rep(0, length(var_name))
  names(var_v) <- var_name
  var_v[imp_df$var_name] <- imp_df$importance

  return(var_v)

}

#' Get multi-omics weights
#' @rdname get_multi_weights
#' @export

get_multi_weights <- function(mod_list, dat.list, y = NULL, weighted = F, lambda = 1, use_depth = F,
                              parallel = T, normalized = T, calc = "Both", yprob = 1,
                              w = NULL, cores = max(detectCores() - 2,20), seed = -5,  ...){

  mod_names <- names(mod_list)
  if(length(lambda) == 1) {
    lambda <- rep(lambda, length(mod_list))
    names(lambda) <- mod_names
  }

  results <- get_results(mod_list = mod_list, parallel = parallel, weighted = weighted, normalized = F, use_depth = use_depth,
                         calc = calc, lambda = lambda, w = w, cores = cores, yprob = yprob, seed = seed)

  net <- purrr::map(results, "net")
  weight_l <- purrr::map(results, "wl")
  weight_l_init <- purrr::map(results, "wl_init")
  
  if(length(weight_l) > 1){
    weight_list <- plyr::llply(
      names(dat.list),
      .fun = function(i){
        w <- purrr::map(weight_l, i)
        w <- purrr::compact(w)
        w <- (Reduce("+", w))/length(w)
        if(normalized) {
          w <- w/sqrt(sum(w^2))
        }
        w
      }
    )
    names(weight_list) <- names(dat.list)
  } else {
    weight_list <- weight_l[[1]]
    wl <- plyr::llply(
      names(weight_list),
      .fun = function(i){
        w <- weight_list[[i]]
        if(normalized) {
          w <- w/sqrt(sum(w^2))
        }
        w
      }
    )
    names(wl) <- names(weight_list)
    weight_list <- wl
    weight_list <- weight_list[names(dat.list)]
  }

  out <- list(weight_list = weight_list,
              weight_list_init = weight_l_init,
              net = net)

  return(out)
}


cal_freq <- function(mod, net){

  x_freq <- colMeans(mod$var.used != 0)
  freq <- list(X = x_freq)

  if(is.null(mod$yvar) | class(mod)[3] == "class+"){
    return(freq)
  } else {

    f <- rep(0, ncol(mod$yvar))
    names(f) <- colnames(mod$yvar)
    y_freq <- purrr::map(
      net,
      ~{
        id <- .[["Y_id"]]
        tb <- table(id)
        names(tb)
      }
    )
    y_freq <- table(unlist(y_freq))/mod$ntree
    f[names(y_freq)] <- y_freq

    freq <- c(freq,list(Y = f))

    return(freq)
  }


}

step_two_weight <- function(two_step, rm_noise, normalized, weight_l, dat.list, freq_ls, s){

  if(two_step){
    if(length(freq_ls) > 1){
      freq <- plyr::llply(
        names(dat.list),
        .fun = function(i){

          fr <- purrr::map(freq_ls, i)
          fr <- purrr::compact(fr)
          Reduce("+", fr)/length(fr)
        }
      )
      names(freq) <- names(dat.list)
    } else {
      freq <- freq_ls[[1]]
    }
  }

  if(length(weight_l) > 1){
    weight_list <- plyr::llply(
      names(dat.list),
      .fun = function(i){
        w <- purrr::map(weight_l, i)
        w <- purrr::compact(w)
        w <- (Reduce("+", w))/length(w)
        if(two_step){
          x <- freq[[i]]
          # o <- sort(x, decreasing = T)[1:ceiling(a[i] * length(x))]
          # x[x < min(o)] <- 0
          w <- w * x
        }
        if(rm_noise){
          x <- s * sd(w)
          w[w < x] <- 0
          w[w > x] <- w[w > x] - x
          
        } 
        if(normalized) {
          w <- w/sqrt(sum(w^2))
        }
        
        w
      }
    )
    names(weight_list) <- names(dat.list)
  } else {
    weight_list <- weight_l[[1]]
    wl <- plyr::llply(
      names(weight_list),
      .fun = function(i){
        w <- weight_list[[i]]
        if(two_step){

          x <- freq[[i]]
          # o <- sort(x, decreasing = T)[1:ceiling(a[i] * length(x))]
          # x[x < min(o)] <- 0
          w <- w * x
        }
        if(rm_noise){
          x <- s * sd(w)
          w[w < x] <- 0
          w[w > x] <- w[w > x] - x

        } 
        if(normalized) {
          w <- w/sqrt(sum(w^2))
        }
        
        w
      }
    )
    names(wl) <- names(weight_list)
    weight_list <- wl
    weight_list <- weight_list[names(dat.list)]
  }

  weight_list
}

get_results <- function(mod_list, parallel, normalized = F, weighted = F, use_depth = F, calc, lambda, w = NULL, yprob = 1, cores = detectCores() - 2, seed = -5){

  mod_names <- names(mod_list)
  plyr::llply(
    mod_names,
    .fun = function(m_name){

      mod <- mod_list[[m_name]]
      l <- lambda[m_name]
      if(!is.null(w)) {
        w0 <- w[[gsub("_.*", "", m_name)]]
      } else {w0 <- NULL}
      results <- get_imp_forest(mod, parallel = parallel, normalized = normalized, 
                                weighted = weighted, calc = calc, lambda = l, w = w0, 
                                yprob = yprob, cores = cores, use_depth = use_depth, seed = seed)
      wl <- results$imp_ls
      wl_init <- results$imp_ls_init
      net <- results$net
      m_name_sep <- unlist(stringr::str_split(m_name, "_"))

      names(wl) <- rev(m_name_sep)
      # freq <- cal_freq(mod, net)
      # names(freq) <- rev(m_name_sep)
      return(list(
        wl = wl,
        wl_init = wl_init,
        net = net
        # freq = freq
      ))
    },.parallel = F
  )


}
