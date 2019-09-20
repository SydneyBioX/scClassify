# functions get ensemble results



# function to transfer the error
alpha <- function(e) {
  e[e < 0.001] <- 0.001
  e[e > 0.999] <- 0.999
  log((1 - e)/e)
}



# function to get the ensemble res
getEnsembleRes <- function(testRes, trainRes, exclude = NULL, weighted_ensemble = TRUE) {
#
#    trainRes <- testRes$train
  # Get the weighted
  if (weighted_ensemble) {
    if (!is.null(exclude)) {
      keep_method <- !grepl(paste(exclude, collapse = "|"), names(trainRes))
      trainRes <- trainRes[keep_method]
      weight <- weight[keep_method]
    }

    errClass_train <- do.call(rbind, lapply(trainRes, function(x) table(x$classifyRes)/length(x$classifyRes)))
    weighted_train <- alpha(1 - (errClass_train[, 1] + errClass_train[, 2]))

  } else {
    weighted_train <- NULL
  }

  # ensembleRes <- lapply(testRes, function(x) {
  #   calEnsembleRes(x,
  #                  exclude = exclude,
  #                  weight = weighted_train)
  # })
  ensembleRes <-   calEnsembleRes(testRes,
                                  exclude = exclude,
                                  weight = weighted_train)
  return(ensembleRes)

}

# Function to combine the ensemble res

calEnsembleRes <- function(res, exclude = NULL, weight = NULL){


  if (!is.null(exclude)) {
    keep_method <- !grepl(paste(exclude,collapse = "|"), names(res))
    res <- res[keep_method]
    weight <- weight[keep_method]
  }

  if (is.null(weight)) {
    weight <- rep(1, length(res))
  }

  ensembleResMat <- do.call(cbind, lapply(res, function(x) x$predRes))

  ensembleRes <- apply(ensembleResMat, 1, function(x) {
    names(x) <- colnames(ensembleResMat)
    # keep <- x!="unassigned"
    keep <- rep(TRUE, length(x))
    if (sum(keep) == 0) {
      data.frame(cellTypes = "unassigned", scores = 0)
    }else{
      getResByWeights(x[keep], weight[keep])
    }

  })

  ensembleRes <- do.call(rbind, ensembleRes)
  return(ensembleRes)
}




# function to calculate the weight for each base classifier

getResByWeights <- function(res, weight) {
  resType <- unique(res)
  if (length(resType) == 1) {
    final <- resType
    scores <- 1
  } else {


    mat <- sapply(1:length(resType), function(i){
      ifelse(res %in% resType[i], 1, 0) * weight
    })
    colnames(mat) <- resType
    rownames(mat) <- names(res)
    mat_colMeans <- Matrix::colMeans(mat)
    final <- names(mat_colMeans)[which.max(mat_colMeans)]
    scores <- max(mat_colMeans)/sum(mat_colMeans)
  }

  return(data.frame(cellTypes = final, scores = scores))
}
