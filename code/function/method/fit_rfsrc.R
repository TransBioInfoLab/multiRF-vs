#' Fit random forest model and create model lists
#'
#' @param X A data frame that consider as predictor set
#' @param Y A data frame that consider as response set. If Y = NULL, an unsupervised RF is conducted
#' @param type Select the type of RF model. The default is regression. Can select from "regression", "classification", and "unsupervised"
#' @param nodedepth rfsrc parameter. Maximum depth to which a tree should be grown. Parameter is ignored by default.
#' @param ntree rfsrc parameter. Number of trees.
#' @param forest.wt rfsrc parameter. Creates an nxn matrix which can be used for prediction and constructing customized estimators. The default is "all".
#' Can select from "all", "inbag", "oob", TRUE, or FALSE. Setting forest.wt = TRUE is equivalent to forest.wt = "inbag".
#' @param proximity rfsrc parameter. Proximity of cases as measured by the frequency of sharing the same terminal node. This is an nxn matrix, which can be large.
#' Choices are inbag, oob, all, TRUE, or FALSE. Setting proximity = TRUE is equivalent to proximity = "inbag". The default is "all".
#' @inheritParams randomForestSRC::rfsrc
#' @return A model list
#'
#' @export fit_rfsrc
#'
#'
fit_rfsrc <-  function(X, Y = NULL, type = "regression", nodedepth = NULL,
                       ntree = 200, forest.wt = "all", proximity = "all", ytry = NULL, ...){

  X <- data.frame(X)

  if(type == "classification"){

    mrf <- rfsrc(
      Y ~ .,
      data = data.frame(Y = Y, X),
      membership = T,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      na.action = "na.impute",
      statistics = T,
      proximity = proximity,
      ...
    )
  }

  if(type == "regression") {

    Y <- data.frame(Y)

    mrf <- rfsrc(
      get.mv.formula(colnames(Y)),
      data = data.frame(X,Y),
      membership = T,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      na.action = "na.impute",
      statistics = T,
      proximity = proximity,
      ytry = ytry,
      ...
    )

  }
  if(type == "unsupervised"){

    if(is.null(ytry)) ytry <- 10

    mrf <- rfsrc(
      data = X,
      membership = T,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      statistics = T,
      proximity = proximity,
      ytry = ytry,
      ...
    )
  }

  return(mrf)
}


#' @export
#' @param dat.list A list that contains multi-omics datasets with samples in columns and features in rows. Samples should be matched in each dataset.
#' @param connect_list A pre-defined connection list between datasets. If provided, variable selection will be conducted based on this connection list.
#' If not provided, the algorithm will find the optimal connection between datasets.
#' @rdname fit_rfsrc
#'
fit_multi_rfsrc <- function(dat.list, connect_list = NULL, var.wt = NULL, yprob = 1, ...){

  if(is.null(connect_list)){

    mod_l <- fullConnect(dat.list, ...)

  } else {

    mod_l <- plyr::llply(
      connect_list,
      .fun = function(d){

        dat_fit <- dat.list[d]
        if(!is.null(var.wt)) {
          varwt <- var.wt[d]
        } else {
          varwt <- NULL
        }

        if(length(dat_fit) == 1){

          ytry <- min(ceiling(ncol(dat_fit[[1]]) * yprob), 20)
          mod <- fit_rfsrc(dat_fit[[1]], xvar.wt = varwt[[1]], ytry = ytry, ...)
        } else {

          if(is.null(varwt)) {
            ytry <- min(ceiling(ncol(dat_fit[[1]]) * yprob), 5000)
          } else {
            ytry <- min(ceiling(ncol(dat_fit[[1]]) * yprob), length(varwt[[1]][varwt[[1]] != 0]))
          }

          if(ytry == ncol(dat_fit[[1]])) ytry <- NULL

          mod <- fit_rfsrc(dat_fit[[2]], dat_fit[[1]], xvar.wt = varwt[[2]], yvar.wt = varwt[[1]], ytry = ytry, ...)
        }


        return(mod)
      }
    )

    names(mod_l) <- purrr::map(connect_list, ~paste0(., collapse = "_"))
  }

  return(mod_l)

}
