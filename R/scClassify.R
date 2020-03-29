#' @title Train and test scClassify model
#'
#' @param exprsMat_train A matrix of log-transformed expression matrix of reference dataset
#' @param cellTypes_train A vector of cell types of reference dataset
#' @param exprsMat_test A list or a matrix indicates the expression matrices of the query datasets
#' @param cellTypes_test A list or a vector indicates cell types of the query datasets (Optional).
#' @param tree A vector indicates the method to build hierarchical tree, set as "HOPACH" by default.
#' This should be one of "HOPACH" and "HC" (using hclust).
#' @param selectFeatures A vector indicates the gene selection method, set as "limma" by default.
#' This should be one or more of "limma", "DV", "DD", "chisq", "BI".
#' @param algorithm A vector indicates the KNN method that are used, set as
#' "WKNN" by default. Thisshould be one or more of "WKNN", "KNN", "DWKNN".
#' @param similarity A vector indicates the similarity measure that are used, set as "pearson" by default.
#' This should be one or more of "pearson",  "spearman", "cosine", "jaccard",
#' "kendall", "binomial", "weighted_rank","manhattan"
#' @param cutoff_method A vector indicates the method to cutoff the correlation distribution.
#' Set as "dynamic" by default.
#' @param weighted_ensemble A logical input indicates in ensemble learning,
#' whether the results is combined by a
#' weighted score for each base classifier.
#' @param weights A vector indicates the weights for ensemble
#' @param weighted_jointClassification A logical input indicates in joint classification
#' using multiple training datasets,
#' whether the results is combined by a weighted score for each training model.
#' @param cellType_tree A list indicates the cell type tree provided by user.
#' (By default, it is NULL) (Only for one training data input)
#' @param k An integer indicates the number of neighbour
#' @param topN An integer indicates the top number of features that are selected
#' @param hopach_kmax An integer between 1 and 9 specifying the maximum number of
#' children at each node in the HOPACH tree.
#' @param pSig A numeric indicates the cutoff of pvalue for features
#' @param prob_threshold A numeric indicates the probability threshold for KNN/WKNN/DWKNN.
#' @param cor_threshold_static A numeric indicates the static correlation threshold.
#' @param cor_threshold_high A numeric indicates the highest correlation threshold
#' @param returnList A logical input indicates whether the output will be class of list
#' @param parallel A logical input indicates whether running in paralllel or not
#' @param BPPARAM  A \code{BiocParallelParam} class object
#' from the \code{BiocParallel} package is used. Default is SerialParam().
#' @param verbose A logical input indicates whether the intermediate steps will be printed
#'
#' @return A list of the results, including testRes storing the results of the testing information,
#' and trainRes storing the training model inforamtion.
#'
#' @author Yingxin Lin
#'
#' @examples
#'
#' data("scClassify_example")
#' xin_cellTypes <- scClassify_example$xin_cellTypes
#' exprsMat_xin_subset <- scClassify_example$exprsMat_xin_subset
#' wang_cellTypes <- scClassify_example$wang_cellTypes
#' exprsMat_wang_subset <- scClassify_example$exprsMat_wang_subset
#'
#' scClassify_res <- scClassify(exprsMat_train = exprsMat_xin_subset,
#' cellTypes_train = xin_cellTypes,
#' exprsMat_test = list(wang = exprsMat_wang_subset),
#' cellTypes_test = list(wang = wang_cellTypes),
#' tree = "HOPACH",
#' algorithm = "WKNN",
#' selectFeatures = c("limma"),
#' similarity = c("pearson"),
#' returnList = FALSE,
#' verbose = FALSE)
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam
#' @export


scClassify <- function(exprsMat_train = NULL,
                       cellTypes_train = NULL,
                       exprsMat_test = NULL,
                       cellTypes_test = NULL,
                       tree = "HOPACH",
                       algorithm = "WKNN",
                       selectFeatures = "limma",
                       similarity = "pearson",
                       cutoff_method = c("dynamic", "static"),
                       weighted_ensemble = FALSE,
                       weights = NULL,
                       weighted_jointClassification = TRUE,
                       cellType_tree = NULL,
                       k = 10,
                       topN = 50,
                       hopach_kmax = 5,
                       pSig = 0.01,
                       prob_threshold = 0.7,
                       cor_threshold_static = 0.5,
                       cor_threshold_high = 0.7,
                       returnList = TRUE,
                       parallel = FALSE,
                       BPPARAM = BiocParallel::SerialParam(),
                       verbose = FALSE) {


  # check input
  if (is.null(exprsMat_train) | is.null(cellTypes_train) | is.null(exprsMat_test)) {
    stop("exprsMat_train or cellTypes_train or exprsMat_test is NULL!")
  }

  if (!is.null(cellTypes_test)) {
    if ("character" %in% is(cellTypes_test)) {
      if (length(cellTypes_test) != ncol(exprsMat_test)) {
        stop("Length of testing cell types does not match
             with number of column of testing expression matrix")
      }
    }

    if ("list" %in% is(cellTypes_test)) {
      if (sum(unlist(lapply(cellTypes_test, length)) !=
              unlist(lapply(exprsMat_test, ncol))) != 0) {
        stop("Length of testing cell types does not match
             with number of column of testing expression matrix")
      }
    }

  }

  if ("list" %in% is(exprsMat_train) ) {
    if (sum(unlist(lapply(cellTypes_train, length)) !=
            unlist(lapply(exprsMat_train, ncol))) != 0) {
      stop("Length of training cell types does not match with
           number of column of training expression matrix")
    }
  }else {
    if (length(cellTypes_train) != ncol(exprsMat_train)) {
      stop("Length of training cell types does not match with
           number of column of training expression matrix")
    }
  }

  if ("list" %in% is(exprsMat_train)) {
    if (any(lapply(cellTypes_train, function(x) any(table(x) == 1)))) {
      stop("There is cell type with only one cell,
           please check cellTypes_train")
    }
  }else {
    if (any(table(cellTypes_train) == 1)) {
      stop("There is cell type with only one cell,
           please check cellTypes_train")
    }
  }




  tree <- match.arg(tree, c("HOPACH", "HC"), several.ok = FALSE)
  selectFeatures <- match.arg(selectFeatures,
                              c("limma", "DV", "DD", "chisq", "BI"),
                              several.ok = TRUE)

  algorithm <- match.arg(algorithm,
                         c("WKNN", "KNN", "DWKNN"),
                         several.ok = TRUE)

  similarity <- match.arg(similarity, c("pearson",  "spearman",
                                        "cosine", "jaccard", "kendall",
                                        "weighted_rank","manhattan"), several.ok = TRUE)
  cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"))



  # To check if need to run weighted ensemble learning

  if ((length(selectFeatures) > 1 | length(algorithm) > 1 | length(similarity) > 1) ) {
    if (weighted_ensemble) {
      weighted_ensemble <- TRUE
      ensemble <- TRUE

      if (verbose) {
        cat("Performing weighted ensemble learning... \n")
      }

    } else {
      weighted_ensemble <- FALSE
      ensemble <- TRUE

      if (verbose) {
        cat("Performing unweighted ensemble learning... \n")
      }

    }
  } else {
    weighted_ensemble <- FALSE
    ensemble <- FALSE

    if (verbose) {
      cat("Ensemble learning is disabled... \n")
    }

  }

  # To check if need to run weighted joint classification

  if ("list" %in% is(exprsMat_train) &
      length(exprsMat_train) > 1 &
      weighted_jointClassification) {
    cat("Performing weighted joint classification \n")
    weighted_jointClassification <- TRUE
  } else {
    weighted_jointClassification <- FALSE
  }



  ensemble_methods <- as.matrix(expand.grid(similarity = similarity,
                                            algorithm = algorithm,
                                            features = selectFeatures))

  if (!is.null(weights)) {
    if (length(weights) != nrow(ensemble_methods)) {
      stop("The length of weights is not equal to
           the number of combination of ensemble methods")
    }
  }



  # calculate the weights for train model
  if (weighted_jointClassification | (weighted_ensemble & is.null(weights))) {
    weightsCal = TRUE
  } else {
    weightsCal = FALSE
  }



  ### train_scClassify
  trainRes <- train_scClassify(exprsMat_train,
                               cellTypes_train,
                               tree = tree,
                               selectFeatures = selectFeatures,
                               topN = topN,
                               hopach_kmax = hopach_kmax,
                               pSig = pSig,
                               cellType_tree = cellType_tree,
                               weightsCal = weightsCal,
                               parallel = parallel,
                               BPPARAM = BPPARAM,
                               verbose = verbose,
                               k = k,
                               prob_threshold = prob_threshold,
                               cor_threshold_static = cor_threshold_static,
                               cor_threshold_high = cor_threshold_high,
                               algorithm = algorithm,
                               similarity = similarity,
                               cutoff_method = cutoff_method
                               )





  if (verbose) {
    cat("Predicting using followings parameter combinations: \n")
    print(ensemble_methods)
  }

  ### if there are multiple testing datasets
  if (verbose) {
    cat("=====================  Start classifying on test dataset  ========================== \n")
  }

  if ("list" %in% is(exprsMat_test)) {
    testRes <- list()

    for (testDataset_idx in seq_len(length(exprsMat_test))) {

      if (verbose) {
        cat("Predicting: ")
        print(names(exprsMat_test)[testDataset_idx])
      }

      if ("list" %in% is(exprsMat_train)) {
        # for the case there are multiple training datasets
        #
        predictRes <- list()


        for (train_list_idx in seq_len(length(exprsMat_train))) {

          if (verbose) {
            cat("Training using ")
            cat(names(exprsMat_train)[train_list_idx], "\n")
          }

          predictRes[[train_list_idx]] <- predict_scClassify(exprsMat_test = exprsMat_test[[testDataset_idx]],
                                                             trainRes  = trainRes[[train_list_idx]],
                                                             cellTypes_test = cellTypes_test[[testDataset_idx]],
                                                             k = k,
                                                             prob_threshold = prob_threshold,
                                                             cor_threshold_static = cor_threshold_static,
                                                             cor_threshold_high = cor_threshold_high,
                                                             algorithm = algorithm,
                                                             features = selectFeatures,
                                                             similarity = similarity,
                                                             cutoff_method = cutoff_method,
                                                             weighted_ensemble = weighted_ensemble,
                                                             weights = weights,
                                                             parallel = parallel,
                                                             BPPARAM = BPPARAM,
                                                             verbose = verbose)

        }
        names(predictRes) <- paste("Trained_by", names(trainRes), sep = "_")
      }else {
        predictRes <- predict_scClassify(exprsMat_test = exprsMat_test[[testDataset_idx]],
                                         trainRes  = trainRes,
                                         cellTypes_test = cellTypes_test[[testDataset_idx]],
                                         k = k,
                                         prob_threshold = prob_threshold,
                                         cor_threshold_static = cor_threshold_static,
                                         cor_threshold_high = cor_threshold_high,
                                         algorithm = algorithm,
                                         features = selectFeatures,
                                         similarity = similarity,
                                         cutoff_method = cutoff_method,
                                         weighted_ensemble = weighted_ensemble,
                                         weights = weights,
                                         parallel = parallel,
                                         BPPARAM = BPPARAM,
                                         verbose = verbose)


      }

      testRes[[testDataset_idx]] <- predictRes
    }

    names(testRes) <- names(exprsMat_test)
  }else{
    # else only one dataset as a matrix in the test
    if ("list" %in% is(exprsMat_train)) {
      # for the case there are multiple training datasets
      #
      testRes <- list()
      for (train_list_idx in seq_len(length(exprsMat_train))) {

        if (verbose) {
          cat("Training using ")
          cat(names(exprsMat_train)[train_list_idx], "\n")
        }
        testRes[[train_list_idx]] <- predict_scClassify(exprsMat_test = exprsMat_test,
                                                        trainRes  = trainRes[[train_list_idx]],
                                                        cellTypes_test = cellTypes_test,
                                                        k = k,
                                                        prob_threshold = prob_threshold,
                                                        cor_threshold_static = cor_threshold_static,
                                                        cor_threshold_high = cor_threshold_high,
                                                        algorithm = algorithm,
                                                        features = selectFeatures,
                                                        similarity = similarity,
                                                        cutoff_method = cutoff_method,
                                                        weighted_ensemble = weighted_ensemble,
                                                        weights = weights,
                                                        parallel = parallel,
                                                        BPPARAM = BPPARAM,
                                                        verbose = verbose)

      }
      names(testRes) <- paste("Trained_by", names(trainRes), sep = "_")
    }else {
      predictRes <- predict_scClassify(exprsMat_test = exprsMat_test,
                                       trainRes  = trainRes,
                                       cellTypes_test = cellTypes_test,
                                       k = k,
                                       prob_threshold = prob_threshold,
                                       cor_threshold_static = cor_threshold_static,
                                       cor_threshold_high = cor_threshold_high,
                                       algorithm = algorithm,
                                       features = selectFeatures,
                                       similarity = similarity,
                                       cutoff_method = cutoff_method,
                                       weighted_ensemble = weighted_ensemble,
                                       weights = weights,
                                       parallel = parallel,
                                       BPPARAM = BPPARAM,
                                       verbose = verbose)
      testRes <- list(test = predictRes)
    }




  }


  if (returnList) {
    return(list(testRes = testRes, trainRes = trainRes))
  } else {

    if ("list" %in% is(exprsMat_train)) {
      trainClassList <- list()
      for (train_list_idx in seq_len(length(trainRes))) {
        trainClassList[[train_list_idx]] <- scClassifyTrainModel(
          name = names(trainRes)[train_list_idx],
          cellTypeTree = trainRes[[train_list_idx]]$cutree_list,
          cellTypeTrain = as.character(trainRes[[train_list_idx]]$cellTypes_train),
          features = names(trainRes[[train_list_idx]]$hierarchyKNNRes),
          model = trainRes[[train_list_idx]]$hierarchyKNNRes,
          modelweights = trainRes[[train_list_idx]]$modelweights,
          metaData = S4Vectors::DataFrame())

      }
      trainClassList <- scClassifyTrainModelList(trainClassList)
    } else {
      trainClassList <- scClassifyTrainModel(
        name = "training",
        cellTypeTree = trainRes$cutree_list,
        cellTypeTrain = as.character(trainRes$cellTypes_train),
        features = names(trainRes$hierarchyKNNRes),
        model = trainRes$hierarchyKNNRes,
        modelweights = trainRes$modelweights,
        metaData = S4Vectors::DataFrame())
    }


    return(list(testRes = testRes, trainRes = trainClassList))


  }


}
