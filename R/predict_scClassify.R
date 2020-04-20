#' Testing scClassify model
#'
#' @param exprsMat_test A list or a matrix indicates the log-transformed
#' expression matrices of the query datasets
#' @param trainRes A `scClassifyTrainModel` or a `list` indicates
#' scClassify trained model
#' @param cellTypes_test A list or a vector indicates cell types
#' of the qurey datasets (Optional).
#' @param k An integer indicates the number of neighbour
#' @param prob_threshold A numeric indicates the probability threshold
#' for KNN/WKNN/DWKNN.
#' @param cor_threshold_static A numeric indicates the static
#' correlation threshold.
#' @param cor_threshold_high A numeric indicates the highest
#' correlation threshold
#' @param features A vector indicates the gene selection method,
#' set as "limma" by default.
#' This should be one or more of "limma", "DV", "DD", "chisq", "BI".
#' @param algorithm A vector indicates the KNN method that are used,
#' set as "WKNN" by default.
#' This should be one or more of "WKNN", "KNN", "DWKNN".
#' @param similarity A vector indicates the similarity measure that are used,
#' set as "pearson" by default.
#' This should be one or more of "pearson",  "spearman", "cosine",
#' "jaccard", "kendall", "binomial", "weighted_rank","manhattan"
#' @param cutoff_method A vector indicates the method to cutoff the
#' correlation distribution. Set as "dynamic" by default.
#' @param weighted_ensemble A logical input indicates in ensemble learning,
#' whether the results is combined by a
#' weighted score for each base classifier.
#' @param weights A vector indicates the weights for ensemble
#' @param parallel A logical input indicates whether running in paralllel or not
#' @param BPPARAM  A \code{BiocParallelParam} class object
#' from the \code{BiocParallel} package is used. Default is SerialParam().
#' @param verbose A logical input indicates whether the intermediate steps
#' will be printed
#' @return list of results
#'
#' @examples
#' data("scClassify_example")
#' wang_cellTypes <- scClassify_example$wang_cellTypes
#' exprsMat_wang_subset <- scClassify_example$exprsMat_wang_subset
#' data("trainClassExample_xin")
#'
#' pred_res <- predict_scClassify(exprsMat_test = exprsMat_wang_subset,
#' trainRes = trainClassExample_xin,
#' cellTypes_test = wang_cellTypes,
#' algorithm = "WKNN",
#' features = c("limma"),
#' similarity = c("pearson"),
#' prob_threshold = 0.7,
#' verbose = TRUE)
#'
#' @author Yingxin Lin
#'
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam bplapply
#'
#' @export

# function to predict multiple feature + distance combinations

predict_scClassify <- function(exprsMat_test,
                               trainRes,
                               cellTypes_test = NULL,
                               k = 10,
                               prob_threshold = 0.7,
                               cor_threshold_static = 0.5,
                               cor_threshold_high = 0.7,
                               features = "limma",
                               algorithm = "WKNN",
                               similarity = "pearson",
                               cutoff_method = c("dynamic", "static"),
                               weighted_ensemble = FALSE,
                               weights = NULL,
                               parallel = FALSE,
                               BPPARAM = BiocParallel::SerialParam(),
                               verbose = FALSE){

    # checking input

    algorithm <- match.arg(algorithm, c("WKNN", "KNN", "DWKNN"),
                           several.ok = TRUE)
    similarity <- match.arg(similarity, c("pearson",  "spearman",
                                          "cosine", "jaccard", "kendall",
                                          "weighted_rank","manhattan"),
                            several.ok = TRUE)
    cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"),
                               several.ok = TRUE)


    if ("scClassifyTrainModelList" %in% is(trainRes)) {
        stop("For a list of training model,
         please use predict_scClassifyJoint()
             instead to get joint training results.")
    }

    if (!any(c("scClassifyTrainModel", "list") %in% is(trainRes))) {
        stop("Wrong trainRes input.
             Need to be either scClassifyTrainModel or list")
    }

    if ((length(features) > 1 | length(algorithm) > 1 |
         length(similarity) > 1) ) {
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



    ensemble_methods <- as.matrix(expand.grid(similarity = similarity,
                                              algorithm = algorithm,
                                              features = features))

    if (!is.null(weights)) {
        if (length(weights) != nrow(ensemble_methods)) {
            stop("The length of weights is not equal to
           the number of combination of ensemble methods")
        }
    }



    if (parallel) {

        predictRes <- BiocParallel::bplapply(seq_len(nrow(ensemble_methods)),
                                             function(em){
                                                 predict_scClassifySingle(exprsMat_test =  exprsMat_test,
                                                                          trainRes = trainRes,
                                                                          cellTypes_test = cellTypes_test,
                                                                          k = k,
                                                                          prob_threshold = prob_threshold,
                                                                          cor_threshold_static = cor_threshold_static,
                                                                          cor_threshold_high = cor_threshold_high,
                                                                          features = ensemble_methods[em, 3],
                                                                          algorithm = ensemble_methods[em, 2],
                                                                          similarity = ensemble_methods[em, 1],
                                                                          cutoff_method = cutoff_method,
                                                                          verbose = verbose)
                                             }, BPPARAM = BPPARAM)


    }else{
        predictRes <- list()
        for (em in seq(nrow(ensemble_methods))) {

            if (verbose) {
                cat("Using parameters: \n")
                print(ensemble_methods[em, ])
            }

            predictRes[[em]] <- predict_scClassifySingle(exprsMat_test =  exprsMat_test,
                                                         trainRes = trainRes,
                                                         cellTypes_test = cellTypes_test,
                                                         k = k,
                                                         prob_threshold = prob_threshold,
                                                         cor_threshold_static = cor_threshold_static,
                                                         cor_threshold_high = cor_threshold_high,
                                                         features = ensemble_methods[em, 3],
                                                         algorithm = ensemble_methods[em, 2],
                                                         similarity = ensemble_methods[em, 1],
                                                         cutoff_method = cutoff_method,
                                                         verbose = verbose)
        }


    }
    names(predictRes) <- paste(ensemble_methods[, 1],
                               ensemble_methods[, 2],
                               ensemble_methods[, 3],
                               sep = "_")

    if (is.null(weights)) {
        if ("list" %in% is(trainRes)) {
            weights <- trainRes$modelweights
        }

        if ("scClassifyTrainModel" %in% is(trainRes)) {
            weights <- modelweights(trainRes)
        }

    }

    if (verbose) {
        cat("weights for each base method: \n")
        print(weights)
    }



    if (ensemble) {
        ensembleRes <- getEnsembleRes(predictRes,
                                      weights,
                                      exclude = NULL,
                                      weighted_ensemble = weighted_ensemble)


        predictRes$ensembleRes <- ensembleRes
    }


    return(predictRes)
}


#' Testing scClassify model (joint training)
#'
#' @param exprsMat_test A list or a matrix indicates the expression matrices of the testing datasets
#' @param trainRes A `scClassifyTrainModel` or a `list` indicates scClassify training model
#' @param cellTypes_test A list or a vector indicates cell types of the testing datasets (Optional).
#' @param k An integer indicates the number of neighbour
#' @param prob_threshold A numeric indicates the probability threshold for KNN/WKNN/DWKNN.
#' @param cor_threshold_static A numeric indicates the static correlation threshold.
#' @param cor_threshold_high A numeric indicates the highest correlation threshold
#' @param features A vector indicates the method to select features, set as "limma" by default.
#' This should be one or more of "limma", "DV", "DD", "chisq", "BI".
#' @param algorithm A vector indicates the KNN method that are used, set as "WKNN" by default.
#' This should be one or more of "WKNN", "KNN", "DWKNN".
#' @param similarity A vector indicates the similarity measure that are used,
#' set as "pearson" by default.
#' This should be one or more of "pearson",  "spearman", "cosine", "jaccard", "kendall",
#' "binomial", "weighted_rank","manhattan"
#' @param cutoff_method A vector indicates the method to cutoff the correlation distribution.
#' Set as "dynamic" by default.
#' @param parallel A logical input indicates whether running in paralllel or not
#' @param BPPARAM  A \code{BiocParallelParam} class object
#' from the \code{BiocParallel} package is used. Default is SerialParam().
#' @param verbose A logical input indicates whether the intermediate steps will be printed
#' @return list of results
#' @author Yingxin Lin
#'
#' @examples
#' data("scClassify_example")
#' wang_cellTypes <- scClassify_example$wang_cellTypes
#' exprsMat_wang_subset <- scClassify_example$exprsMat_wang_subset
#' data("trainClassExample_xin")
#' data("trainClassExample_wang")
#'
#' trainClassExampleJoint <- scClassifyTrainModelList(trainClassExample_wang,
#' trainClassExample_xin)
#'
#' pred_res <- predict_scClassifyJoint(exprsMat_test = exprsMat_wang_subset,
#' trainRes = trainClassExampleJoint,
#' cellTypes_test = wang_cellTypes,
#' algorithm = "WKNN",
#' features = c("limma"),
#' similarity = c("pearson"),
#' prob_threshold = 0.7,
#' verbose = FALSE)
#'
#' table(pred_res_joint$jointRes$cellTypes, wang_cellTypes)
#'
#' @importFrom methods is slot
#' @importFrom BiocParallel SerialParam
#' @export

predict_scClassifyJoint <- function(exprsMat_test,
                                    trainRes,
                                    cellTypes_test,
                                    k = 10,
                                    prob_threshold = 0.7,
                                    cor_threshold_static = 0.5,
                                    cor_threshold_high = 0.7,
                                    features = "limma",
                                    algorithm = "WKNN",
                                    similarity = "pearson",
                                    cutoff_method = c("dynamic", "static"),
                                    parallel = FALSE,
                                    BPPARAM = BiocParallel::SerialParam(),
                                    verbose = FALSE){

    # checking input

    algorithm <- match.arg(algorithm, c("WKNN", "KNN", "DWKNN"),
                           several.ok = TRUE)
    similarity <- match.arg(similarity, c("pearson",  "spearman",
                                          "cosine", "jaccard", "kendall",
                                          "weighted_rank","manhattan"),
                            several.ok = TRUE)
    cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"),
                               several.ok = TRUE)


    if (!"scClassifyTrainModelList" %in% is(trainRes)) {
        stop("trainRes needs to be a scClassifyTrainModelList object")
    }



    predictResList <- list()


    for (trainListIdx in seq_len(length(methods::slot(trainRes, "listData")))) {

        if (verbose) {
            cat("Training using ")
            cat(names(methods::slot(trainRes, "listData"))[trainListIdx], "\n")
        }

        predictResList[[trainListIdx]] <- predict_scClassify(exprsMat_test = exprsMat_test,
                                                             trainRes = trainRes[[trainListIdx]],
                                                             cellTypes_test = cellTypes_test,
                                                             k = k,
                                                             prob_threshold = prob_threshold,
                                                             cor_threshold_static = cor_threshold_static,
                                                             cor_threshold_high = cor_threshold_high,
                                                             features = features,
                                                             algorithm = algorithm,
                                                             similarity = similarity,
                                                             cutoff_method = cutoff_method,
                                                             parallel = parallel,
                                                             BPPARAM = BPPARAM,
                                                             verbose = verbose)
    }


    names(predictResList) <- paste("Trained_by",
                                   names(methods::slot(trainRes, "listData")),
                                   sep = "_")

    jointRes <- getJointRes(predictResList, trainRes)

    rownames(jointRes) <- colnames(exprsMat_test)


    predictResList$jointRes <- jointRes


    return(predictResList)
}


predict_scClassifySingle <- function(exprsMat_test,
                                     trainRes,
                                     cellTypes_test,
                                     k = 10,
                                     prob_threshold = 0.7,
                                     cor_threshold_static = 0.5,
                                     cor_threshold_high = 0.7,
                                     features = "limma",
                                     algorithm = c("WKNN", "KNN", "DWKNN"),
                                     similarity = c("pearson",  "spearman",
                                                    "cosine", "jaccard",
                                                    "kendall", "weighted_rank",
                                                    "manhattan"),
                                     cutoff_method = c("dynamic", "static"),
                                     verbose = TRUE){


    if (!is.null(cellTypes_test)) {
        if (length(cellTypes_test) != ncol(exprsMat_test)) {
            stop("Length of testing cell types does not match
           with number of column of testing expression matrix")
        }
    }

    if (is.null(rownames(exprsMat_test))) {
        stop("rownames of exprsMat_test is NULL")
    }

    if (is.null(colnames(exprsMat_test)) |
        sum(duplicated(colnames(exprsMat_test))) != 0) {
        stop("colnames of exprsMat_test is NULL or not unique")
    }

    if (all(exprsMat_test %% 1 == 0)) {
        warning("exprsMat_test looks like a count matrix
            (scClassify requires a log-transformed normalised input)")
    }

    # check input
    algorithm <- match.arg(algorithm, c("WKNN", "KNN", "DWKNN"))
    similarity <- match.arg(similarity, c("pearson",  "spearman",
                                          "cosine", "jaccard", "kendall",
                                          "weighted_rank","manhattan"))
    cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"))


    if (ncol(exprsMat_test) < 100 & cutoff_method == "dynamic") {
        warning("Number of cells in test data is small.
                The cutoff method is changed to static.")
        cutoff_method <- "static"
    }


    if ("scClassifyTrainModelList" %in% is(trainRes)) {
        stop("For a list of training model,
         please use predict_scClassifyJoint instead to
             get joint training results.")
    }

    if (!any(c("scClassifyTrainModel", "list") %in% is(trainRes))) {
        stop("wrong trainRes input.
             Need to be either scClassifyTrainModel or list")
    }

    if ("scClassifyTrainModel" %in% is(trainRes)) {

        if (!features %in% features(trainRes)) {
            stop("The selected features are not trained in the provided model!")
        }


        # get the input from train model

        levelModel <- model(trainRes)[[features]]$model
        levelHVG <- model(trainRes)[[features]]$hvg
        cutree_list <- cellTypeTree(trainRes)
        cellTypes_train <- cellTypeTrain(trainRes)
    } else {

        if (!features %in% names(trainRes$hierarchyKNNRes)) {
            stop("The selected features are not trained in the provided model!")
        }


        # get the input from train model

        levelModel <- trainRes$hierarchyKNNRes[[features]]$model
        levelHVG <- trainRes$hierarchyKNNRes[[features]]$hvg
        cutree_list <- trainRes$cutree_list
        cellTypes_train <- trainRes$cellTypes_train

    }

    pred <- list()

    # For each level
    for (i in seq_len(length(levelModel))) {
        # If this level is not NULL (Not Level1)
        if (!is.null(levelModel[[i]])) {


            #If not all the cells are not unassigned in the parent nodes
            if (sum(pred[[i - 1]] != 0) != 0) {
                pred_level <- list()

                for (j in seq_len(length(levelModel[[i]]))) {

                    # If the model of level i-1, cells labeled as j is NOT "no model"
                    if (!"character" %in% is(levelModel[[i]][[j]])) {

                        # Select the cells that are going to classified
                        # (according to what they are classified in last level)
                        predIdx <- which(pred[[i - 1]] == j)

                        if (length(predIdx) != 0) {

                            # features that are in the test dataset
                            common_HVG <- intersect(rownames(exprsMat_test), levelHVG[[i]][[j]])
                            exprsMat_toTest <- exprsMat_test[common_HVG, predIdx, drop = FALSE]
                            exprsMat_toTest <- exprsMat_toTest[Matrix::rowSums(exprsMat_toTest) != 0, ,
                                                               drop = FALSE]

                            if (nrow(exprsMat_toTest) < 5) {

                                message(paste("There are only", nrow(exprsMat_toTest),
                                              "selected genes in reference data expressed in query data"))
                                pred_level[[j]] <- rep(0, length(predIdx))
                            }else{

                                # Calculate the similarity
                                corMat <- calculateSimilarity(exprsMat_train = Matrix::t(levelModel[[i]][[j]]$train[ ,rownames(exprsMat_toTest)]),
                                                              exprsMat_test = exprsMat_toTest,
                                                              similarity = similarity)

                                # Different algorithm
                                if (algorithm == "KNN") {
                                    predRes <- KNNcor(corMat = corMat,
                                                      subLevelModel = levelModel[[i]][[j]],
                                                      cutoff_method = cutoff_method,
                                                      k = k,
                                                      prob_threshold = prob_threshold,
                                                      cor_threshold_static = cor_threshold_static,
                                                      cor_threshold_high = cor_threshold_high,
                                                      topLevel = FALSE,
                                                      verbose = verbose)
                                }
                                if (algorithm == "WKNN") {
                                    predRes <- WKNNcor(corMat = corMat,
                                                       subLevelModel = levelModel[[i]][[j]],
                                                       cutoff_method = cutoff_method,
                                                       k = k,
                                                       prob_threshold = prob_threshold,
                                                       cor_threshold_static = cor_threshold_static,
                                                       cor_threshold_high = cor_threshold_high,
                                                       topLevel = FALSE,
                                                       verbose = verbose)
                                }

                                if (algorithm == "DWKNN") {
                                    predRes <- DWKNNcor(corMat = corMat,
                                                        subLevelModel = levelModel[[i]][[j]],
                                                        cutoff_method = cutoff_method,
                                                        k = k,
                                                        prob_threshold = prob_threshold,
                                                        cor_threshold_static = cor_threshold_static,
                                                        cor_threshold_high = cor_threshold_high,
                                                        topLevel = FALSE,
                                                        verbose = verbose)
                                }

                                pred_level[[j]] <- predRes$predRes
                            }

                        }

                    }

                    # Else, the model of (level i-1, cells labeled as j) IS "no model"
                    # maintain the same class.
                    else {
                        predIdx <- which(pred[[i - 1]] == j)
                        # check the label of current level based on the label of last level
                        pred_level[[j]] <- as.factor(rep(unique(cutree_list[[i]][cutree_list[[i - 1]] == j]),
                                                         length(predIdx)))
                        names(pred_level[[j]]) <- colnames(exprsMat_test)[predIdx]
                    }
                }



                # Get the predict results for level i, change it to the numeric

                pred[[i]] <- unlist(lapply(pred_level, function(x){
                    cellNames <- names(x)
                    x <- as.numeric(as.character(x))
                    names(x) <- cellNames
                    x
                }))

                # reorder the prediction results to consistent
                #with the exprsMat_test
                pred[[i]] <- pred[[i]][colnames(exprsMat_test)]

                names(pred[[i]]) <- colnames(exprsMat_test)
                pred[[i]] <- as.numeric(as.character(pred[[i]]))
                # there will be NA since they are unassigned
                #from the last level,
                #and therefore are not predicted in this level
                pred[[i]][is.na(pred[[i]])] <- 0
                names(pred[[i]]) <- colnames(exprsMat_test)


            }
            else{

                # else, if all the cells in paraent nodes are unassigned
                pred[[i]] <- as.factor(rep(0, ncol(exprsMat_test)))
                names(pred[[i]]) <- colnames(exprsMat_test)
            }

        }else{
            # If this level is NULL (Level 1)
            pred[[i]] <- as.factor(rep(1, ncol(exprsMat_test)))
            names(pred[[i]]) <- colnames(exprsMat_test)
        }


    }

    predMat <- do.call(cbind, pred)
    predMat <- vapply(seq_len(ncol(predMat)),
                      function(x) getPredRes(predMat, cutree_list, x),
                      character(nrow(predMat)))

    predRes <- apply(predMat, 1, function(x){
        unAssIdx <- which(x == "unassigned")
        if (length(unAssIdx) != 0) {
            if (min(unAssIdx) >= 3) {
                x[min(unAssIdx) - 1]
            }else{
                "unassigned"
            }
        }else{
            x[length(x)]
        }
    })

    if (is.null(cellTypes_test)) {
        return(list(predRes = predRes, predLabelMat = predMat))
    }else{

        classify_res <-  ClassifyError(cellTypes_pred = predRes,
                                       cellTypes_test = cellTypes_test,
                                       cellTypes_train = cellTypes_train)

        if (verbose) {
            print(table(classify_res)/length(classify_res))
        }

        return(list(predRes = predRes, predLabelMat = predMat,
                    classifyRes = classify_res))

    }

}



# Function to calculate similarity
#' @importFrom proxyC simil
#' @importFrom proxy as.dist
#' @importFrom proxy dist
#' @importFrom methods as


calculateSimilarity <- function(exprsMat_train,
                                exprsMat_test,
                                similarity = c("pearson",  "spearman",
                                               "cosine", "jaccard", "kendall",
                                               "weighted_rank","manhattan")) {


    similarity <- match.arg(similarity, c("pearson",  "spearman",
                                          "cosine", "jaccard", "kendall",
                                          "weighted_rank","manhattan"))

    if ("dgCMatrix"  %in% is(exprsMat_test) &
        !"dgCMatrix" %in% is(exprsMat_train)) {
        exprsMat_train <- methods::as(exprsMat_train, "dgCMatrix")
    }

    if (!"dgCMatrix" %in% is(exprsMat_test) &
        "dgCMatrix" %in% is(exprsMat_train)) {
        exprsMat_test <- methods::as(exprsMat_test, "dgCMatrix")
    }


    if (similarity == "cosine") {

        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- proxyC::simil(Matrix::t(exprsMat_train),
                                    Matrix::t(exprsMat_test),
                                    method = "cosine")
        } else {
            corMat <- 1 - as.matrix(proxy::dist(t(as.matrix(exprsMat_train)),
                                                t(as.matrix(exprsMat_test)),
                                                method = "cosine"))
        }

        corMat[is.na(corMat) | is.infinite(corMat)] <- min(corMat)

    } else if (similarity == "kendall") {



        corMat <- stats::cor(as.matrix(exprsMat_train),
                             as.matrix(exprsMat_test), method = "kendall")
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1

    } else if (similarity == "jaccard") {

        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- proxyC::simil(Matrix::t(exprsMat_train),
                                    Matrix::t(exprsMat_test),
                                    method = "jaccard")
        } else {
            corMat <- 1 -
                as.matrix(proxy::dist(t(as.matrix(exprsMat_train > 0)),
                                      t(as.matrix(exprsMat_test > 0)),
                                      method = "Jaccard"))

        }
        corMat[is.na(corMat) | is.infinite(corMat)] <- min(corMat)
    }else if (similarity == "weighted_rank") {

        corMat <- wtd_rank2(as.matrix(exprsMat_train),
                            as.matrix(exprsMat_test), method = "pearson")
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1

    }else if (similarity == "manhattan") {
        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- 1 - as.matrix(proxy::dist(t(as.matrix(exprsMat_train)),
                                                t(as.matrix(exprsMat_test)),
                                                method = "Manhattan"))
        } else {
            corMat <- 1 - as.matrix(proxy::dist(t(as.matrix(exprsMat_train)),
                                                t(as.matrix(exprsMat_test)),
                                                method = "Manhattan"))
        }
        corMat[is.na(corMat) | is.infinite(corMat)] <- min(corMat)

    }else if (similarity == "spearman") {

        corMat <- stats::cor(as.matrix(exprsMat_train),
                             as.matrix(exprsMat_test), method = "spearman")
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1

    }else if (similarity == "pearson") {

        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- proxyC::simil(Matrix::t(exprsMat_train),
                                    Matrix::t(exprsMat_test),
                                    method = "correlation")
        } else {
            corMat <- stats::cor(as.matrix(exprsMat_train),
                                 as.matrix(exprsMat_test), method = "pearson")
        }
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1

    }else{

        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- proxyC::simil(Matrix::t(exprsMat_train),
                                    Matrix::t(exprsMat_test),
                                    method = "correlation")
        } else {
            corMat <- stats::cor(as.matrix(exprsMat_train),
                                 as.matrix(exprsMat_test), method = "pearson")
        }
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1
    }

    return(corMat)
}


# Function to caculate weighted rank correlation
# the codes are modified from the wtd_rank function in dismay package

wtd_rank2 <- function(mat1, mat2 = NULL, method = "pearson") {

    method <- match.arg(method, c("pearson", "spearman"))

    ranks1 <- apply(mat1, 2, function(x) rank(-x, ties.method = "average"))
    # weight the ranks
    # calculate the savage scores
    n1 <- nrow(mat1)
    reciprocals1 <- 1 / seq_len(n1)
    savage1 <- vapply(seq_len(n1),
                      function(i) sum(reciprocals1[i:n1]),
                      numeric(1L))
    # replace each rank with the savage score
    savages1 <- ranks1
    savages1[] <- savage1[ranks1]

    if (!is.null(mat2)) {
        ranks2 <- apply(mat2, 2, function(x) rank(-x, ties.method = "average"))
        # weight the ranks
        # calculate the savage scores
        n2 <- nrow(mat2)
        reciprocals2 <- 1 / seq_len(n2)
        savage2 <- vapply(seq_len(n2),
                          function(i) sum(reciprocals2[i:n2]),
                          numeric(1L))
        # replace each rank with the savage score
        savages2 <- ranks2
        savages2[] <- savage2[ranks2]

        cor <- stats::cor(savages1, savages2, method = method)
    } else {
        cor <- stats::cor(savages1, method = method)
    }


    # calculate pearson correlation

    return(cor)
}

# Function to get the prediction labels according to the tree
getPredRes <- function(predMat, cutree_list, level){
    res <- predMat[,level]
    for (i in seq_len(length(predMat[,level]))) {
        if (predMat[i,level] == 0) {
            res[i] <- "unassigned"
        }else{
            res[i] <- paste(names(cutree_list[[level]])[cutree_list[[level]]
                                                        %in% predMat[i,level]],
                            collapse = "_")
        }

    }
    return(res)
}


# Functions to classify the error for the predicted cell types

ClassifyError <- function(cellTypes_pred, cellTypes_test, cellTypes_train){

    errClass <- c("correct", "correctly unassigned",
                  "intermediate", "incorrectly unassigned",
                  "error assigned", "misclassified")


    if (length(cellTypes_pred) != length(cellTypes_test)) {
        stop("wrong input")
    }
    train_ref <- unique(cellTypes_train)
    res <- vapply(seq_len(length(cellTypes_pred)), function(i){
        if (cellTypes_test[i] %in% train_ref) {
            if (cellTypes_pred[i] %in% c("unassigned", "Unassigned")) {
                "incorrectly unassigned"
            } else if (cellTypes_pred[i] == "intermediate") {
                "intermediate"
            }else{
                if (cellTypes_test[i] == cellTypes_pred[i]) {
                    "correct"
                }else if (grepl(cellTypes_test[i], cellTypes_pred[i])) {
                    "intermediate"
                }
                else{
                    "misclassified"
                }
            }
        }else{
            if (cellTypes_pred[i] %in% c("unassigned","Unassigned")) {
                "correctly unassigned"
            }else{
                "error assigned"
            }
        }
    }, character(1L))

    res <- factor(res, levels = errClass)
    return(res)
}



