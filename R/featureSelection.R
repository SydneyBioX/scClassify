#' @importFrom methods new
#' @importFrom Cepo Cepo topGenes
#' @importFrom scrapper scoreMarkers

featureSelection <- function(exprsMat,
                             trainClass,
                             sampleID = NULL,
                             feature = c("DM", "DV", "DD", "chisq", "BI",
                                         "Cepo"),
                             topN = 50,
                             pSig = 0.001
){

    feature <- match.arg(feature, c("DM", "DV", "DD", "chisq", "BI", "Cepo"),
                         several.ok = FALSE)


    if (feature == "DV") {
        tt <- doDV(exprsMat, trainClass)
        tt <- lapply(tt, function(x)sort(x))
        res <- Reduce(union, lapply(tt, function(t)
            names(t)[seq_len(min(topN, sum(t < pSig)))]))
    } else if (feature == "DD") {
        tt <- doDD(exprsMat, trainClass)
        tt <- lapply(tt, function(x)sort(x))
        res <- Reduce(union, lapply(tt, function(t)
            names(t)[seq_len(min(topN, sum(t < pSig)))]))
    } else if (feature == "chisq") {
        tt <- doChisSquared(exprsMat, trainClass)
        tt <- lapply(tt, function(x)sort(x))
        res <- Reduce(union, lapply(tt, function(t)
            names(t)[seq_len(min(topN, sum(t < pSig)))]))
        #
    } else if (feature == "BI") {
        tt <- doBI(exprsMat, trainClass)
        tt <- lapply(tt, function(x)x)
        res <- Reduce(union, lapply(tt, function(t)
            names(t)[seq_len(topN)]))
    } else if (feature == "Cepo") {
        tt <- Cepo::Cepo(as.matrix(exprsMat), trainClass, exprsPct = 0.05)
        res <- Reduce(union, Cepo::topGenes(tt, n = topN))
    } else{
        effectSizes <- scrapper::scoreMarkers(exprsMat, trainClass, sampleID, block.weight.policy = "equal")
        res <- Reduce(union, mapply(function(typeD, propDetect)
        {
          bestNames <- rownames(exprsMat)[order(typeD$mean, decreasing = TRUE)]
          bestNames <- bestNames[propDetect$mean > 0.05]
          bestNames[seq_len(topN)]
        }, effectSizes[["cohens.d"]], effectSizes[["delta.detected"]], SIMPLIFY = FALSE))

    }

    return(res)
}

doDV <- function(exprsMat, cellTypes){


    cellTypes <- droplevels(as.factor(cellTypes))
    tt <- list()
    for (i in seq_len(nlevels(cellTypes))) {
        tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))

        meanPct <- do.call(cbind, lapply(c(0,1), function(i){
            Matrix::rowSums(exprsMat[,
                                     tmp_celltype == i,
                                     drop = FALSE] > 0)/sum(tmp_celltype == i)
        }))


        posNeg <- (meanPct[,2] - meanPct[,1]) > 0.05
        # print(sum(posNeg))
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
    for (i in seq_len(nlevels(cellTypes))) {
        tmp_celltype <- ifelse(cellTypes == levels(cellTypes)[i], 1, 0)


        meanPct <- do.call(cbind, lapply(c(0,1), function(i){
            Matrix::rowSums(exprsMat[,
                                     tmp_celltype == i,
                                     drop = FALSE] > 0)/sum(tmp_celltype == i)
        }))

        posNeg <- (meanPct[,2] - meanPct[,1]) > 0.05
        # print(sum(posNeg))
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
    for (i in seq_len(nlevels(cellTypes))) {
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
    for (i in seq_len(nlevels(cellTypes))) {
        tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))

        pi <- table(tmp_celltype)/length(tmp_celltype)

        agg_mean <- do.call(cbind, lapply(c(0,1), function(i){
            Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
        }))

        agg_sd2 <- do.call(cbind, lapply(c(0,1), function(i){
            apply(exprsMat[, tmp_celltype == i, drop = FALSE], 1, stats::var)
        }))

        bi <- abs(agg_mean[,2] - agg_mean[,1])/sqrt(pi[1]*agg_sd2[,1] +
                                                        pi[2]*agg_sd2[,2])

        bi <- unlist(bi)
        names(bi) <- rownames(exprsMat)
        bi <- bi[order(bi, decreasing = TRUE)]
        tt[[i]] <- bi
    }

    return(tt)


}
