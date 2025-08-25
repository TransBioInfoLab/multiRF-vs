#' Find optimal connections among multi-omics data
#'
#' @param dat.list A list containing multi-omics datasets with samples in columns and features in rows.
#' Samples should be matched across all datasets.
#' @param var_prop Proportion of variance explained by PC datasets when finding optimal connections. Default is 0.1.
#' @param keep_prop Pre-defined proportion of connections to retain. Default is median.
#' @param return_score Logical value determining whether to return the connections score (scaled RCV) matrix.
#' @param direct Logical parameter specifying whether to keep both directions in the connection list for optimal connections.
#'
#' @inheritParams randomForestSRC::rfsrc
#' @return model connection vector indicates the direction of connect between datasets. The connection is indicated as: responses_predictors
#' @export findConnection
#' @import randomForestSRC
# ---------------------------------------------------------------------------------
# Find optimal connections among multi-omics data
# ---------------------------------------------------------------------------------
findConnection <- function(dat.list, var_prop = .1, return_score = F, top = 2, ...){
  
  # # Create dimension reduction matrices
  reduced_list <- purrr::map(dat.list, ~getEmbed(., var_prop = var_prop))

  
  # Fit full connection model
  mod <- fullConnect(reduced_list, ...)
  
  # Get oob error
  oob_err <- plyr::laply(
    names(mod),
    .fun = function(i){
      get_r_sq(mod[[i]])
    }
  )
  
  names(oob_err) <- names(mod)
  
  dat_combn <- combn(names(dat.list), m = 2)
  
  oob_err_selected <- sapply(1:ncol(dat_combn), function(d) {
    mod_name1 <- paste0(dat_combn[,d][1], "_", dat_combn[,d][2])
    mod_name2 <- paste0(dat_combn[,d][2], "_", dat_combn[,d][1])
    v <- c(oob_err[mod_name1], oob_err[mod_name2])
    v[which.min(v)]
  })
  model_connection <- names(sort(oob_err_selected)[1:top])
  mod_list <- stringr::str_split(model_connection, "_")
  # # Note that the first var is response and the second var is predictor
  # if(is.null(keep_prop)){
  #   m_oob <- median(oob_err)
  #   model_connection <- names(oob_err)[oob_err < m_oob]
  # 
  # } else {
  #   oob_sort <- sort(oob_err)
  #   m_oob <- oob_sort[floor(length(oob_err) * keep_prop)]
  #   model_connection <- names(oob_sort)[1:floor(length(oob_err) * keep_prop)]
  # }
  
  # if(!direct){
  #   mod_list <- stringr::str_split(model_connection, "_")
  #   oob_err_selected <- oob_err[oob_err < m_oob]
  #   idx <- c()
  #   for(i in 1:(length(mod_list) - 1)){
  #     for(j in (i+1):length(mod_list)){
  #       if(setequal(mod_list[[i]], mod_list[[j]])){
  #         idx <- c(idx, which(oob_err_selected == max(oob_err_selected[c(i,j)])))
  #       }
  #     }
  #   }
  #   model_connection <- model_connection[-idx]
  # }
  
  if(return_score){
    
    d <- matrix(0, nrow = length(dat.list), ncol = length(dat.list))
    dimnames(d) <- list(names(dat.list), names(dat.list))
    d[row(d) != col(d)] <- oob_err
    model_connection <- list(model_connection = model_connection,
                             score = d)
    
  }
  return(model_connection)
  
}

#' @export
#' @rdname findConnection
#' @inheritParams randomForestSRC::rfsrc

fullConnect <- function(dat.list, seed = seed, ...){

  dat_names <- names(dat.list)
  mod_l <- plyr::llply(
    dat_names,
    .fun = function(d){
      response_d <- dat.list[[d]]
      predict_d <- dat.list[!names(dat.list) %in% d]

      mod <- plyr::llply(predict_d, .fun = function(pred){
        fit_rfsrc(X = pred, Y = response_d, seed = seed, ...)
      })

      mod_names <- paste0(d, "_", names(predict_d))
      names(mod) <- mod_names

      return(mod)
    }
  )

  mod_l <- Reduce(c, mod_l)

  return(mod_l)

}


# Get embedding of one data
getEmbed <- function(dat, var_prop = 0.6){

  dat <- scale(dat, scale = F)
  s <- kpca(~.,data.frame(dat))
  vars <- cumsum((s@eig^2) / sum(s@eig^2))
  dim <- which(vars > var_prop)[1]
  s@pcv[,1:dim]

}

get_r_sq <- function(mod){
  
  fw <- mod$forest.wt
  
  ex <- mean(colMeans((fw %*% as.matrix(mod$xvar) - as.matrix(mod$xvar))^2)/matrixStats::colVars(as.matrix(mod$xvar)))
  if(is.null(mod$yvar)){
    return(ex)
  } else {
    if(class(mod)[3] == "class+"){
      ey <- na.omit(mod$err.rate)
    } else {
      ey <- mean(colMeans((fw %*% as.matrix(mod$yvar) - as.matrix(mod$yvar))^2)/matrixStats::colVars(as.matrix(mod$yvar)))
    }
  }
  
  ey #+ ex
  
}
