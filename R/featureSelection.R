#' @importFrom methods new
#' @import limma


featureSelection <- function(exprsMat,
                             trainClass,
                             feature = c("limma", "DV", "DD", "chisq", "BI"),
                             topN = 50,
                             pSig = 0.001
                             ){

  feature <- match.arg(feature, c("limma", "DV", "DD", "chisq", "BI"), several.ok = FALSE)


 if (feature == "DV") {
    tt <- doDV(exprsMat, trainClass)
    tt <- lapply(tt, function(x)sort(x))
    res <- Reduce(union, lapply(tt, function(t) names(t)[1:max(min(topN, sum(t < pSig)), 30)]))
  } else if (feature == "DD") {
    tt <- doDD(exprsMat, trainClass)
    tt <- lapply(tt, function(x)sort(x))
    res <- Reduce(union, lapply(tt, function(t) names(t)[1:max(min(topN, sum(t < pSig)), 30)]))
  } else if (feature == "chisq") {
    tt <- doChisSquared(exprsMat, trainClass)
    tt <- lapply(tt, function(x)sort(x))
    res <- Reduce(union, lapply(tt, function(t) names(t)[1:max(min(topN, sum(t < pSig)), 30)]))
    #
  } else if (feature == "BI") {
    tt <- doBI(exprsMat, trainClass)
    tt <- lapply(tt, function(x)x)
    res <- Reduce(union, lapply(tt, function(t) names(t)[1:topN]))
  }


  else{
    tt <- doLimma(exprsMat, trainClass)
    res <- Reduce(union, lapply(tt, function(t)
      rownames(t[t$logFC > 0 & (t$meanPct.2 - t$meanPct.1) > 0.05 & t$adj.P.Val < pSig,])[1:topN]))

  }

  return(res)
}


#' @importFrom limma eBayes lmFit
#' @importFrom methods new

doLimma <- function(exprsMat, cellTypes, exprs_pct = 0.05){

  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))
    design <- stats::model.matrix(~tmp_celltype)


    meanExprs <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
    }))

    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowSums(exprsMat[, tmp_celltype == i, drop = FALSE] > 0)/sum(tmp_celltype == i)
    }))

    keep <- meanPct[,2] > exprs_pct

    y <- methods::new("EList")
    y$E <- exprsMat[keep, ]
    fit <- limma::lmFit(y, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)



    if (!is.null(tt[[i]]$ID)) {
      tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
      rownames(tt[[i]]) <- tt[[i]]$ID
    }

    tt[[i]]$meanExprs.1 <- meanExprs[rownames(tt[[i]]), 1]
    tt[[i]]$meanExprs.2 <- meanExprs[rownames(tt[[i]]), 2]
    tt[[i]]$meanPct.1 <- meanPct[rownames(tt[[i]]), 1]
    tt[[i]]$meanPct.2 <- meanPct[rownames(tt[[i]]), 2]
  }



  return(tt)


}




doDV <- function(exprsMat, cellTypes){


  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))

    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowSums(exprsMat[, tmp_celltype == i, drop = FALSE] > 0)/sum(tmp_celltype == i)
    }))


    posNeg <- (meanPct[,2] - meanPct[,1]) > 0.05
    print(sum(posNeg))
    exprsMat_filt <- exprsMat[posNeg,]
    tt[[i]] <- apply(exprsMat_filt, 1, function(x) {
      df <- data.frame(gene = x, cellTypes = as.factor(tmp_celltype))
      stats::bartlett.test(gene~cellTypes, df)$p.value
    })

    tt[[i]] <- stats::p.adjust(tt[[i]], method = "BH")
  }



  return(tt)


}

doDD <- function(exprsMat, cellTypes){

  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- ifelse(cellTypes == levels(cellTypes)[i], 1, 0)


    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowSums(exprsMat[, tmp_celltype == i, drop = FALSE] > 0)/sum(tmp_celltype == i)
    }))

    posNeg <- (meanPct[,2] - meanPct[,1]) > 0.05
    print(sum(posNeg))
    exprsMat_filt <- exprsMat[posNeg,]
    tt[[i]] <- apply(exprsMat_filt, 1, function(x) {
      x1 <- x[tmp_celltype == 0]
      x2 <- x[tmp_celltype == 1]
      stats::ks.test(x1, x2, alternative = "greater")$p.value
    })



    tt[[i]] <- stats::p.adjust(tt[[i]], method = "BH")
  }



  return(tt)


}



doChisSquared <- function(exprsMat, cellTypes, threshold = 1){


  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))


    zerosMat <- ifelse(exprsMat > threshold, 1, 0)

    tt[[i]] <- apply(zerosMat,1,  function(x){
      tab <- c()
      for (i in c(0,1)) {
        tmp <- factor(x[tmp_celltype == i], levels = c(0, 1))
        tab <- rbind(tab, table(tmp))
      }


      suppressWarnings(stats::chisq.test(tab)$p.value)


    })




    tt[[i]] <- stats::p.adjust(tt[[i]], method = "BH")
  }



  return(tt)


}






doBI <- function(exprsMat, cellTypes){
  # Select genes by bimodal index

  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))

    pi <- table(tmp_celltype)/length(tmp_celltype)

    agg_mean <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
    }))

    agg_sd2 <- do.call(cbind, lapply(c(0,1), function(i){
      apply(exprsMat[, tmp_celltype == i, drop = FALSE], 1, stats::var)
    }))

    bi <- abs(agg_mean[,2] - agg_mean[,1])/sqrt(pi[1]*agg_sd2[,1] + pi[2]*agg_sd2[,2])

    bi <- unlist(bi)
    names(bi) <- rownames(exprsMat)
    bi <- bi[order(bi, decreasing = T)]
    tt[[i]] <- bi
  }

  return(tt)


}
