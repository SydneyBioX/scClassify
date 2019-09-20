#' @title Run sample size calculation for pilot data for reference dataset
#'
#' @param exprsMat A matrix of expression matrix of pilot dataset (log-transformed, or normalised)
#' @param cellTypes A vector of cell types of pilot dataset
#' @param n_list A vector of integer indicates the sample size to run.
#' @param num_repeat An integer indicates the number of run for each sample size will be repeated.
#' @param ncores An integer indicates the number of cores that can be used to run the algorithm.
#' @param subset_test A ogical input indicates whether we used a subset of data (fixed number for each sample size)
#' to test instead of all remaining data. By default, it is FALSE.
#' @param num_test An integer indicates the size of the test data.
#' @param cellType_tree A list indicates the cell type tree (optional), if it is NULL, the accuracy rate is calculate
#' based on the provided cellTypes.
#' @param level An integer indicates the accuracy rate is calculate based on the n-th level from top of cell type tree.
#' If it is NULL (by default), it will be the bottom of the cell type tree.
#' It can not be larger than the total number of levels of the tree.
#' @param ... other parameter from scClassify
#'
#' @importFrom pbmcapply pbmclapply
#' @export


runSampleCal <- function(exprsMat,
                          cellTypes,
                          n_list = c(20, 40, 60, 80, 100, seq(200, 500, 100)),
                          num_repeat = 20,
                          level = NULL,
                          cellType_tree = NULL,
                          ncores = 1,
                          subset_test = FALSE,
                          num_test = NULL,

                          ...) {




  if (length(n_list) < 5) {
    stop("length of n_list provided is too short... wont be enough point to fit the learning curve")
  }




  # If there is an input of cell type tree, relabelled the cell types based on input level

  if (!is.null(cellType_tree)) {
    if (is.null(level)) {
      level <- length(cellType_tree)
    } else {
      if (level <= 0 | level > length(cellType_tree)) {
        stop("The input of level is invalid")
      }

      level <- level + 1
    }

    cellTypes_ind <- cellType_tree[[level]][as.character(cellTypes)]
    cellTypes_relabelled <- sapply(1:length(cellTypes_ind), function(i) {
      paste(names(cellType_tree[[level]])[cellType_tree[[level]] %in% cellTypes_ind[i]],
            collapse = "_")
    })

    cellTypes <- cellTypes_relabelled
    cellType_tree_relabelled <- reLevelCellTypeTree(cellType_tree, level = level)

  } else {
    cellType_tree_relabelled <- NULL
  }


  exprsMat <- as(exprsMat, "dgCMatrix")
  # n_list <- c(20, 40, 60, 80, 100, seq(200, 500, 100))

  res_sub <- list()

  for (i in 1:length(n_list)) {
    print(paste("n=",n_list[i]))
    system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:num_repeat, function(x) {
      tryCatch({
        l <- runSubSampling(exprsMat, cellTypes, n = n_list[[i]],
                            subset_test = subset_test, num_test = num_test,
                            cellType_tree = cellType_tree_relabelled,
                            cutoff_method = "dynamic",
                            prob_threshold = 0.6
                            )

      table(l)/length(l)},
      error = function(e){NULL})
    }, mc.cores = ncores))
    print(do.call(cbind, res_sub[[i]]))
    gc(reset = TRUE)
  }

  names(res_sub) <- n_list

  res_sub <- lapply(res_sub, function(x) unlist(lapply(x, "[[", "correct")))
  res_sub <- res_sub[!is.null(res_sub)]
  accuracy_mat <- do.call(cbind, res_sub)
  return(accuracy_mat)
}




# function to run subsampling


runSubSampling <- function(exprsMat,
                           cellTypes,
                           n = 50,
                           subset_test = F,
                           num_test = 2000,
                           ...){


  # make sure the smallest training type has at least 3 cells.
  n_subset <- round(table(cellTypes)/length(cellTypes)*n)
  n_subset <- ifelse(n_subset < 3, 3, n_subset)

  trainIdx <- unlist(lapply(1:length(n_subset), function(x) sample(which(cellTypes == names(n_subset)[x]), n_subset[x])))

  exprsMat_train <- exprsMat[,trainIdx]
  cellTypes_train <- cellTypes[trainIdx]

  print(table(cellTypes_train))

  if (subset_test) {
    testIdx <- sample(c(1:ncol(exprsMat))[-trainIdx], num_test)
    exprsMat_test <- exprsMat[, testIdx]
    cellTypes_test <- cellTypes[testIdx]
  }else{
    exprsMat_test <- exprsMat[, -trainIdx]
    cellTypes_test <- cellTypes[-trainIdx]
  }


  trainRes <- scClassify(exprsMat_train = exprsMat_train,
                         cellTypes_train = cellTypes_train,
                         exprsMat_test = list(
                           test = exprsMat_test),
                         cellTypes_test = list(
                           test = cellTypes_test),
                         ...)

  acc_cls <- trainRes$testRes$test$pearson_WKNN_limma$classifyRes
  return(acc_cls)
}



reLevelCellTypeTree <- function(cellType_tree, level) {

  cellTypes_newTree <- list()

  newCellTypeNames <- sapply(seq_len(max(cellType_tree[[level]])), function(i) {
    paste(names(cellType_tree[[level]])[cellType_tree[[level]] == i], collapse = "_")
  })

  cellTypes_newTree[[level]] <- 1:max(cellType_tree[[level]])
  names(cellTypes_newTree[[level]]) <- newCellTypeNames

  for (i in 1:(level - 1)) {

    # cellTypes_newTree[[i]] <- cellType_tree[[i]]
    #
    cellTypes_newTree[[i]] <- sapply(seq_len(length(newCellTypeNames)),
                                     function(j)
                                       unique(cellType_tree[[i]][names(cellType_tree[[i]]) %in% unlist(strsplit(newCellTypeNames[j], "_"))])
    )
    names(cellTypes_newTree[[i]]) <- newCellTypeNames

  }

  return(cellTypes_newTree)
}
