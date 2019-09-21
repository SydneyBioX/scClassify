#' Training scClassify model
#'
#' @param exprsMat_train A matrix of expression matrix of reference dataset
#' @param cellTypes_train A vector of cell types of reference dataset
#' @param tree A vector indicates the method to build hierarchical tree, set as "HOPACH" by default.
#' This should be one of "HOPACH" and "HC" (using stats::hclust).
#' @param selectFeatures A vector indicates the gene selection method, set as "limma" by default.
#' This should be one or more of "limma", "DV", "DD", "chisq", "BI".
#' @param topN An integer indicates the top number of features that are selected
#' @param hopach_kmax An integer between 1 and 9 specifying the maximum number of
#' children at each node in the HOPACH tree.
#' @param pSig A numeric indicates the cutoff of pvalue for features
#' @param cellType_tree A list indicates the cell type tree provided by user. (By default, it is NULL)
#' @param parallel A logical input indicates whether the algorihms will run in parallel
#' @param ncores An integer indicates the number of cores that are used
#' @param verbose A logical input indicates whether the intermediate steps will be printed
#' @return list of results
#' @author Yingxin Lin
#'
#' @importFrom stats na.omit
#' @export




train_scClassify <- function(exprsMat_train,
                             cellTypes_train,
                             tree = "HOPACH",
                             selectFeatures = "limma",
                             topN = 50,
                             hopach_kmax = 5,
                             pSig = 0.05,
                             cellType_tree = NULL,
                             parallel = F,
                             ncores = 1,
                             verbose = T){

  if (is.null(rownames(exprsMat_train))) {
    stop("rownames of the exprssion matrix is NULL!")
  }

  if (length(cellTypes_train) != ncol(exprsMat_train)) {
    stop("Length of training cell types does not match with number of column of training expression matrix")
  }

  # Matching the argument of the tree construction method
  # TODO: Allow user to provide their own tree
  tree <- match.arg(tree, c("HOPACH", "HC"), several.ok = FALSE)

  # Matching the argument of feature selection method
  # TODO: 1. Allow user to provide their own markers
  selectFeatures <- match.arg(selectFeatures, c("limma", "DV", "DD", "chisq", "BI"), several.ok = TRUE)


  if (verbose) {
    print("Feature Selection...")
  }


  # Select the features to construct tree
  tt <- doLimma(exprsMat_train, cellTypes_train)
  de <- Reduce(union, lapply(tt, function(t) rownames(t)[1:max(min(50, sum(t$adj.P.Val < 0.001)), 30)]))

  if (verbose) {
    print(paste("Number of genes selected to construct HOPACH tree", length(de)))
  }


  if (is.null(cellType_tree)) {
    # Calculate the centroid matrix for tree construction using selected features
    centroidMat <- do.call(cbind, lapply(unique(cellTypes_train),
                                         function(x) Matrix::rowMeans(as.matrix(exprsMat_train[de, cellTypes_train == x]))))

    colnames(centroidMat) <- unique(cellTypes_train)

    # Constructing the tree using selected tree method
    if (verbose) {
      print("Constructing tree ...")
    }

    cutree_list <- constructTree(centroidMat,
                                 tree = tree,
                                 hopach_kmax = hopach_kmax,
                                 plot = verbose)


  } else {
    cutree_list <- cellType_tree
  }

  if (verbose) {
    print("Training....")
  }



  # if (featureHierachy) {
  #   hierachyKNNRes <-  hierchyKNNcor_prob2(exprsMat_train,
  #                                          cellTypes_train,
  #                                          cutree_list,
  #                                          topN = topN,
  #                                          feature = selectFeatures)
  # }
  #

  if (parallel) {
    hierarchyKNNRes <- pbmcapply::pbmclapply(1:length(selectFeatures), function(ft) hierarchyKNNcor(exprsMat_train,
                                                                             cellTypes_train,
                                                                             cutree_list,
                                                                             feature = selectFeatures[ft],
                                                                             topN = topN,
                                                                             pSig = pSig,
                                                                             verbose = verbose), mc.cores = ncores)
    names(hierarchyKNNRes) <- selectFeatures
  }else{
    hierarchyKNNRes <- list()
    for (ft in 1:length(selectFeatures)) {

      if (verbose) {
        print(paste("============ selecting features by:", selectFeatures[ft], "============"))
      }

      hierarchyKNNRes[[ft]] <- hierarchyKNNcor(exprsMat_train,
                                               cellTypes_train,
                                               cutree_list,
                                               feature = selectFeatures[ft],
                                               topN = topN,
                                               pSig = pSig,
                                               verbose = verbose)
    }

    names(hierarchyKNNRes) <- selectFeatures
  }


  return(list(hierarchyKNNRes = hierarchyKNNRes,
              cutree_list = cutree_list,
              cellTypes_train = cellTypes_train))
}



# Function to construct tree using HOPACH or HC

constructTree <- function(centroidMat,
                          tree = c("HOPACH", "HC"),
                          hopach_kmax = hopach_kmax,
                          plot = T){


  tree <- match.arg(tree, c("HOPACH", "HC"), several.ok = FALSE)



  if (tree == "HOPACH") {

    res <- runHOPACH(data = t(centroidMat),
                     plot = plot,
                     kmax = hopach_kmax)
    cutree_list <- res$cutree_list

  }else{
    distMat <- as.dist(1 - stats::cor(centroidMat))
    hc <- stats::hclust(distMat, method = "complete")
    # plot(hc)
    cutree_list <- cutree_iterative(hc, depth = 1)
  }

  return(cutree_list)
}


# Cut the Hierarchical tree for each level
# depth here indicates the deep the tree is cut

cutree_iterative <- function(hc, depth = 1){

  height_sort <- sort(hc$height, decreasing = T)
  cutree_list <- list()

  for (i in 1:sum(height_sort >= height_sort[round(length(height_sort)*depth)])) {
    cutree_list[[i]] <- cutree(hc, h = height_sort[i])
  }

  # Last level is distinct number
  cutree_list[[length(cutree_list) + 1]] <- 1:length(hc$labels)
  names(cutree_list[[length(cutree_list)]]) <- hc$labels

  return(cutree_list)
}

currentClass <- function(cellTypes, cutree_res){
  cellTypes <- as.character(cellTypes)
  res <- cutree_res[cellTypes]
  return(res)
}

# Function to perform hierarchical feature selection

hierarchyKNNcor <- function(exprsMat,
                            cellTypes,
                            cutree_list,
                            feature = c("limma", "DV", "DD", "chisq", "BI"),
                            topN = 50,
                            pSig = 0.001,
                            verbose = T){
  feature <- match.arg(feature, c("limma", "DV", "DD", "chisq", "BI"))
  numHierchy <- length(cutree_list)
  levelModel <- list()
  levelHVG <- list()
  for (i in 1:numHierchy) {
    # print(paste ("Level", i))


    #Make sure not all same cluster
    if (length(unique(cutree_list[[i]])) != 1) {

      model <- list()
      hvg <- list()
      class_tmp <- currentClass(cellTypes, cutree_list[[i]])
      names(class_tmp) <- colnames(exprsMat)
      for (j in unique(cutree_list[[i - 1]])) {
        # print(paste("Group", j))
        trainIdx <- which(cellTypes %in% names(cutree_list[[i - 1]])[cutree_list[[i - 1]] == j])
        trainClass <- class_tmp[trainIdx]
        if (length(unique(trainClass)) != 1) {

          hvg[[j]] <- featureSelection(exprsMat[,trainIdx],
                                        trainClass,
                                        feature = feature,
                                        topN = topN,
                                        pSig = pSig
                                        )
          # if (verbose) {
          #   print(hvg[[j]])
          # }

          model[[j]] <- list(train = Matrix::t(exprsMat[na.omit(hvg[[j]]), trainIdx, drop = FALSE]),
                             y = as.factor(trainClass))

        }else{
          model[[j]] <- "noModel"
        }
      }
      levelHVG[[i]] <- hvg
      levelModel[[i]] <- model


    }
  }
  res <- list(model = levelModel, hvg = levelHVG)
  return(res)
}
