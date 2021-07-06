# functions get ensemble results



# function to transfer the error
alpha <- function(e) {
    e[e < 0.001] <- 0.001
    e[e > 0.999] <- 0.999
    log((1 - e)/e)
}



getTrainWeights <- function(trainRes) {
    errClass_train <- do.call(rbind, lapply(trainRes, function(x)
        table(x$classifyRes)/length(x$classifyRes)))
    weighted_train <- alpha(1 - (errClass_train[, 1] + errClass_train[, 2]))
    return(weighted_train)
}


# function to get the ensemble res
getEnsembleRes <- function(testRes, weighted_train,
                           exclude = NULL, weighted_ensemble = TRUE) {
    #
    #    trainRes <- testRes$train
    # Get the weighted
    if (weighted_ensemble) {
        if (!is.null(exclude)) {
            keep_method <- !grepl(paste(exclude, collapse = "|"),
                                  names(weighted_train))
            weighted_train <- weighted_train[keep_method]
        }
        
    } else {
        weighted_train <- NULL
    }
    
    ensembleRes <- calEnsembleRes(testRes,
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
    

    ensembleResMat <- do.call(cbind, lapply(res, function(x) x$predRes))
    
    if (is.null(weight)) {
        weight <- rep(1, ncol(ensembleResMat))
    }
    
    ensembleRes <- apply(ensembleResMat, 1, function(x) {
        names(x) <- colnames(ensembleResMat)
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
        
        
        mat <- vapply(seq_len(length(resType)), function(i){
            ifelse(res %in% resType[i], 1, 0) * weight
        }, numeric(length(res)))
        
        colnames(mat) <- resType
        rownames(mat) <- names(res)
        mat_colMeans <- Matrix::colMeans(mat)
        final <- names(mat_colMeans)[which.max(mat_colMeans)]
        scores <- max(mat_colMeans)/sum(mat_colMeans)
    }
    
    return(data.frame(cellTypes = final, scores = scores))
}




# Function to summarise joint classification results

getJointRes <- function(predictResList, trainRes) {
    
    predictJointResList <- lapply(seq_along(predictResList), function(l) {
        
        if ("scClassifyTrainModel" %in% is(trainRes[[l]])) {
            weights_train <- modelweights(trainRes[[l]])
        } else {
            weights_train <- trainRes[[l]]$modelweights
        }
        
        getEnsembleRes(predictResList[[l]],
                       weights_train,
                       exclude = NULL,
                       weighted_ensemble = TRUE)
    })
    
    names(predictJointResList) <- names(predictResList)
    
    
    jointCellType <- do.call(cbind, lapply(predictJointResList, function(x)
        x$cellTypes))
    colnames(jointCellType) <- names(predictJointResList)
    
    
    jointCellType <- apply(jointCellType, 2, function(x) {
        len <- unlist(lapply(strsplit(x, "_"), length))
        x[len > 2] <- "intermediate"
        x
    })
    
    jointWeights <- do.call(cbind, lapply(predictJointResList, function(x)
        x$scores))
    colnames(jointWeights) <- names(predictJointResList)
    
    
    
    jointRes <- lapply(seq_len(nrow(jointCellType)), function(idx) {
        cell_res <- jointCellType[idx, ]
        names(cell_res) <- colnames(jointCellType)
        if (length(cell_res) == 2 & sum(cell_res == "unassigned") == 1) {
            data.frame(cellTypes = cell_res[cell_res != "unassigned"],
                       scores = 0.5)
        } else {
            getResByWeights(cell_res, jointWeights[idx, ])
        }
    })
    
    
    jointRes <- do.call(rbind, jointRes)
    
    jointRes <- data.frame(jointRes)
    
    return(jointRes)
}

