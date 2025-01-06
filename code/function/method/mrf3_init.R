#' Fit an initial MRF model
#'
#' @param dat.list A list containing multi-omics datasets with samples in columns and features in rows. Samples should be matched across datasets.
#' @param ntree Number of trees for fitting MRF model. Default is 300.
#' @param scale Whether to z-standardize each feature. Default is TRUE.
#' @param yprob Probability of response features being selected in each node split. Default is 0.5.
#' @param connect_list Pre-defined connection list between datasets. If provided, variable selection uses this list. If not, the algorithm finds optimal connections between datasets.
#' @param var_prop Proportion of variance explained by PC datasets when finding optimal connections. Default is 0.6.
#' @param direct Logical; determines whether to keep both directions in the connection list for optimal connections.
#' @param lambda Penalizes variables selected only once in a tree. Experimental parameter. Default is 1.
#' @param normalized Logical; determines whether to use normalized variable weights. Default is FALSE.
#' @param use_depth Logical; determines whether to compute the average IMD selected in a tree. Default is FALSE.
#' @param calc Select which weights to calculate: "X", "Y", or "Both". Use when fewer than two datasets are in the model. Default is "Both".
#' @param parallel Logical; determines whether to use parallel computation for weight calculation.
#' @param return_data Whether to return the data list. Default is FALSE.
#' @param cores Number of cores to use for computation.
#'
#' @inheritParams randomForestSRC::rfsrc
#' @return mrf3 object
#' @export
#' @import randomForestSRC
#'
#'
# ---------------------------------------------------------------------------------
# Wrapper function for mrf direct variable selection
# ---------------------------------------------------------------------------------
mrf3_init <- function(dat.list,
                 # Number of trees
                 ntree = 300,
                 scale = T,
                 yprob = .5,
                 # Find connections
                 connect_list = NULL,
                 var_prop = .6,
                 direct = T,
                 keep_prop = NULL,
                 lambda = 1,
                 normalized = F,
                 use_depth = F,
                 calc = "Both",
                 parallel = T,
                 return_data = F,
                 cores = detectCores() - 2,
                 ...){

  if(length(dat.list) == 1) type = "unsupervised"
  if(length(dat.list) > 1) type = "regression"

  # Find connection
  if(is.null(connect_list) & length(dat.list) > 1){
    message("Finding maximum connections..")
    connection <- findConnection(dat.list = dat.list, var_prop = var_prop, direct = direct, keep_prop = keep_prop)
    connect_list <- stringr::str_split(connection, "_")
  }
  if(length(dat.list) == 1){
    connect_list <- list(c(names(dat.list)))
    calc <- "X"
  }
  if(scale) {
    new_dat <- purrr::map(dat.list, ~scale(.))
  } else {
    new_dat <- dat.list
  }


  message("Fitting models..")

  mod_list <- fit_multi_rfsrc(new_dat, connect_list = connect_list, ntree = ntree, type = type, yprob = yprob,...)
  oob_err <- purrr::map(mod_list, ~get_r_sq(.))
  oob_err <- Reduce("+", oob_err)

  w <- NULL

  message("Calculating weights..")
  multi_weights_mod <- get_multi_weights(mod_list = mod_list,
                                         dat.list = new_dat,
                                         y = y,
                                         parallel = parallel,
                                         ntree = ntree1,
                                         normalized = normalized,
                                         calc = calc,
                                         type = type,
                                         lambda = lambda,
                                         w = w,
                                         yprob = yprob,
                                         use_depth = use_depth,
                                         cores = cores,
                                         ...)
  multi_weights <- multi_weights_mod$weight_list

  net <- multi_weights_mod$net
  names(net) <- names(mod_list)

  if(!return_data) dat.list <- NULL
  out <- list(
    weights = multi_weights,
    mod = mod_list,
    oob_err = oob_err,
    type = type,
    connection = connect_list,
    net = net,
    ntree = ntree,
    yprob = yprob,
    dat.list = dat.list,
    weights_ls = multi_weights_mod$weight_list_init
  )

  attr(out, "class") <- "mrf3"

  message("Done!")

  return(out)
}
