#' Create HOPACH tree
#'
#' A function generating HOPACH tree using the average expression matrix for
#' each cell type.
#'
#' @param data A matrix of average expression matrix
#'  (each row indicates the gene, each column indicates the cell type)
#' @param plot Indicate whether plot or not
#' @param kmax Integer between 1 and 9 specifying the maximum number of children
#' at each node in the tree.
#' @return Return a \code{list} where
#' \itemize{
#' \item{cutree_list}: A list indicates the hierarchical cell type tree
#' \item{plot}: A \code{ggplot} visualise the cell type tree
#' }
#' @author Yingxin Lin
#' @importFrom stats cutree median lm
#' @importFrom hopach vectmatrix hdist is.hdist improveordering distancevector
#' as.hdist
#' @importFrom methods as is
#' @importFrom cluster pam
#' @importFrom graphics plot text
#' @import statmod
#'
#' @examples
#'
#' data("scClassify_example")
#' wang_cellTypes <- factor(scClassify_example$wang_cellTypes)
#' exprsMat_wang_subset <- scClassify_example$exprsMat_wang_subset
#' avgMat_wang <- apply(exprsMat_wang_subset, 1, function(x)
#' aggregate(x, list(wang_cellTypes), mean)$x)
#' rownames(avgMat_wang) <- levels(wang_cellTypes)
#' res_hopach <- runHOPACH(avgMat_wang)
#' res_hopach$plot
#'
#' @references van der Laan, M. J. and Pollard, K. S. (2003)
#' ‘A new algorithm for hybrid hierarchical clustering with
#' visualization and the bootstrap’,
#' Journal of Statistical Planning and Inference.
#' doi: 10.1016/S0378-3758(02)00388-9.
#' @export

runHOPACH <- function(data, plot= TRUE, kmax = 5) {

    # Compute the distance matrix using correlation
    dist <- distancematrix(data, d = "cor")

    # Run HOPACH function
    clustresult <- hopach(data, dmat = dist, kmax = kmax)


    # Rearrange HOPACH results
    final_labels <- strsplit(as.character(format(clustresult$final$labels,
                                                 scientific = FALSE)), "")

    # Get each level seperate
    level_list <- list()
    for (i in seq_len(nchar(clustresult$final$labels[1]))) {
        if (i != 1) {
            level_list[[i]] <- paste(level_list[[i - 1]],
                                     unlist(lapply(final_labels, "[[", i)),
                                     sep = "")
        }else{
            level_list[[i]] <- unlist(lapply(final_labels, "[[", i))
        }
    }



    # Generate cutree results
    cutree_list <- list()

    # First level (All ones)
    cutree_list[[1]] <- rep(1, length(rownames(data)))
    names(cutree_list[[1]]) <- rownames(data)

    cutree_list[2:(length(level_list) + 1)] <- lapply(level_list, function(x) {
        x <- as.numeric(as.factor(x))
        names(x) <- rownames(data)
        x
    })

    # Last levels will be all distinct numbers
    cutree_list[[length(cutree_list) + 1]] <-
        seq_len(length(clustresult$final$labels))
    names(cutree_list[[length(cutree_list)]]) <- rownames(data)

    # To plot using ggraph
    if (plot) {

        treePlot <- plotCellTypeTree(cutree_list)
        return(list(cutree_list = cutree_list, plot = treePlot))

    }else{
        return(list(cutree_list = cutree_list))
    }

}

#' To plot cell type tree
#'
#'
#' @param cutree_list A list indicates the hierarchical cell type tree
#' @param group_level Indicate whether plot or not
#'
#' @return A ggplot object visualising the HOPACH tree
#'
#' @importFrom ggraph geom_edge_diagonal geom_node_text geom_node_point ggraph
#' @import igraph
#'
#' @examples
#'
#' data("trainClassExample_xin")
#'
#' plotCellTypeTree(cellTypeTree(trainClassExample_xin))
#'
#' @export


plotCellTypeTree <- function(cutree_list, group_level = NULL) {

    cutree_list <- cutree_list

    if (is.null(group_level)) {
        group_level <- 3
    }

    E_list <- list()
    for (i in seq_len(length(cutree_list))) {
        if (i != 1) {
            parent <- cutree_list[[i - 1]]
            E_list[[i]] <- paste(paste(i - 1, parent, sep = ""),
                                 paste(i, cutree_list[[i]], sep = ""),
                                 sep = "_")
        }
    }

    V_list <- list()
    for (i in seq_len(length(cutree_list))) {
        V_list[[i]] <- paste(i, cutree_list[[i]], sep = "")
        names(V_list[[i]]) <- names(cutree_list[[i]])
    }


    g <- igraph::make_graph(unlist(strsplit(unique(unlist(lapply(E_list,
                                                                 sort))),
                                            "_")))



    ### Which level is used to color the node
    group_level <- group_level
    group_col <- (V_list[[group_level]])
    names(group_col) <- V_list[[length(cutree_list)]]


    V_group <- rep(NA, length(igraph::V(g)))
    names(V_group) <- names(V(g))
    V_group[V_list[[length(V_list)]]] <- group_col[V_list[[length(V_list)]]]
    V(g)$group <- V_group  # note only the bottom level are coloured


    ### Next, get the label for the bottom level

    V_name <- rep(NA, length(igraph::V(g)))
    names(V_name) <- names(V(g))
    V_name[V_list[[length(V_list)]]] <- names(V_list[[length(V_list)]])
    V(g)$name <- V_name




    treePlot <- ggraph::ggraph(g, layout = 'dendrogram') +
        ggraph::geom_edge_diagonal(colour = "black") +
        ggraph::geom_node_text(aes(vjust = 1, hjust = -0.2, label = name),
                               size = 2.7, alpha = 1) +
        ggraph::geom_node_point(aes(filter = leaf, alpha = 0.2, size = 3,
                                    colour = as.factor(V(g)$group))) +
        ggplot2::theme_void() +
        ggplot2::scale_colour_hue() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::coord_flip() +
        NULL
    treePlot
}



#Hierarchical Ordered Partitioning and Collapsing Hybrid (HOPACH)#
#The codes here are originated from HOPACH pacakge with small edits
#
#Reference: ‘A new algorithm for hybrid hierarchical clustering with
# visualization and the bootstrap’,
# Journal of Statistical Planning and Inference.
# doi: 10.1016/S0378-3758(02)00388-9.

labelstomss <- function(labels, dist, khigh = 9, within = "med",
                        between = "med", hierarchical = TRUE) {
    if (!is.vector(labels))
        stop("First arg to labelstomss() must be a vector")

    p = methods::slot(dist, "Size")
    if (length(labels) != p)
        stop("Distance matrix and labels dimensions
             do not agree in labelstomss()")

    unlabels <- sort(unique(labels))
    k <- length(unlabels)
    ss <- NULL
    for (i in seq_len(k)) {
        labs <- (seq_len(p))[labels == unlabels[i]]
        pp <- length(labs)
        if (pp < 3)
            ss[i] <- NA
        else{
            dissvec <- methods::slot(dist[labs,labs], "Data")
            bestk <- silcheck(dissvec,min(khigh,(length(labs) - 1),
                                          na.rm = TRUE), diss = TRUE)[1]
            if (within == "med")
                ss[i] <- median(pam(dissvec, bestk,
                                    diss = TRUE)$silinfo$widths[, 3])
            if (within == "mean")
                ss[i] <- pam(dissvec,bestk,diss = TRUE)$silinfo$avg.width
        }
    }
    if (sum(is.na(ss)) == k)
        out <- NA
    else{
        if (hierarchical == TRUE) {
            parentlab <- trunc(labels/10)
            unparentlab <- sort(unique(parentlab))
            pk <- length(unparentlab)
            pss <- NULL
            for (i in seq_len(pk)) {
                if (between == "med")
                    pss[i] <- median(ss[trunc(unlabels/10) == unparentlab[i]],
                                     na.rm = TRUE)
                if (between == "mean")
                    pss[i] <- mean(ss[trunc(unlabels/10) == unparentlab[i]],
                                   na.rm = TRUE)
            }
            if (between == "med")
                out <- median(pss,na.rm = TRUE)
            if (between == "mean")
                out <- mean(pss,na.rm = TRUE)
        }
        else{
            if (between == "med")
                out <- median(ss,na.rm = TRUE)
            if (between == "mean")
                out <- mean(ss,na.rm = TRUE)
        }
    }
    return(out)
}


silcheck <- function(data, kmax=9, diss=FALSE, echo=FALSE, graph=FALSE)
{
    if ( !diss ) {
        if ( inherits(data, "dist"))
            stop("data argument is a dist object, but diss is FALSE")
        if ( is.matrix(data) && (nrow(data)  ==  ncol(data) ) )
            warning("data argument is square, could be a dissimilarity")
    }
    if ( diss && is.matrix(data) && nrow(data) != ncol(data) )
        stop("should be a dissimilarity matrix - but is not square")

    sil <- NULL
    m <- min(kmax, max((!diss)*(dim(data)[1] - 1),
                       (diss)*(0.5*(sqrt(1 + 8*length(data)) - 1)),
                       na.rm = TRUE))
    if (m < 2)
        out <- c(1, NA)
    else{
        for (i in seq_len((m - 1)))
            sil[i] <- pam(data, k = (i + 1), diss = diss)$silinfo$avg.width
        if (echo)
            cat("best k = ", order(sil)[length(sil)] + 1, ", sil(k) = ",
                round(max(sil),4), "\n")
        if (graph) {
            plot(2:m, sil, type = "n", xlab = "Number of Clusters",
                 ylab = "Average Silhouette")
            text(2:m, sil, 2:m)
        }
        out <- c(order(sil)[length(sil)] + 1,max(sil))
    }
    return(out)
}

#msscheck (mss)
msscheck <- function(dist, kmax=9, khigh=9, within="med", between="med",
                     force=FALSE, echo=FALSE, graph=FALSE) {
    p = methods::slot(dist, "Size")
    if (p < 3)
        out <- c(1,NA)
    else{
        dvec <- methods::slot(dist, "Data")

        kmax <- min(kmax,p - 1,na.rm = TRUE)
        if (force)
            mss <- 0
        else
            mss <- labelstomss(rep(1,p),dist,khigh,within,between)
        for (k in 2:kmax)
            mss[k] <- labelstomss(pam(dvec, k, diss = TRUE)$clust, dist,
                                  khigh, within, between)
        shift <- 0
        if (force) {
            mss <- mss[-1]
            shift <- 1
        }
        if (echo)
            cat("best k = ", order(mss)[1] + shift, ", mss(k) = ",
                round(min(mss),4), "\n")
        if (graph) {
            kmin <- ifelse(force,2,1)
            plot(kmin:kmax, mss, type = "n" ,xlab = "Number of Clusters",
                 ylab = "MSS")
            text(kmin:kmax,mss,kmin:kmax)
        }
        out <- c(order(mss)[1] + shift,min(mss))
    }
    return(out)
}

#2. Utility functions

#counts the number of digits in label
digits <- function(label) {
    label <- label[1]
    count <- 0
    while (label >= 1) {
        count <- count + 1
        label <- label/10
    }
    return(count)
}


# patrick's suggestion...
cutdigits <- function(labels, dig) {
    df <- max(0,digits(labels) - dig)
    return(trunc(labels/(10^df)))
}

#removes trailing zeros from labels
cutzeros <- function(labels) {
    for (i in seq_len(length(labels))) {
        while (trunc(labels[i]/10)*10 == labels[i]) {
            labels[i] <- labels[i]/10
        }
    }
    return(labels)
}

#returns the number of non-zero digits in labels
nonzeros <- function(labels) {
    for (i in seq_len(length(labels))) {
        while (trunc(labels[i]/10)*10 == labels[i]) {
            labels[i] <- labels[i]/10
        }
        labels[i] <- digits(labels[i])
    }
    return(labels)
}


msssplitcluster <- function(clust1, l1, id1, medoid1,
                            med2dist, right, dist1,
                            kmax = 9, khigh = 9,
                            within = "med", between = "med") {
    if (!medoid1)
        warning("Medoid missing - continue to split cluster")
    else{
        if (sum(medoid1 == id1) == 0 & medoid1)
            warning("Medoid not in cluster - continue to split cluster")
    }
    if (is.matrix(clust1)) {
        p1 <- length(clust1[,1])
        n <- length(clust1[1,])
    }
    else p1 <- 1
    if (p1 < 3)
        k1 <- 1
    else{
        l <- length(clust1[,1])
        dissvec <- methods::slot(dist1, "Data")

        kmax <- min(p1 - 1,kmax,na.rm = TRUE)
        khigh <- min(p1 - 1,khigh,na.rm = TRUE)
        k1 <- msscheck(dist1, kmax, khigh, within, between)[1]
        if (k1 > 1) {
            pamobj <- pam(dissvec, k1, diss = TRUE)
            newclussizes <- pamobj$clusinfo[,1]
            newmedoids1 <- id1[pamobj$medoids]
            newlabels1 <- pamobj$clustering
            distnewmedoids <- NULL
            for (j in seq_len(k1))
                distnewmedoids[j] <-
                mean(med2dist[newlabels1 == newlabels1[pamobj$medoids[j]]])
            if (right == 1)
                #ord <- order(distnewmedoids, decreasing = TRUE)
                ord <- rev(order(distnewmedoids))
            else
                ord <- order(distnewmedoids)
            newmedoids1 <- newmedoids1[ord]
            newclussizes <- newclussizes[ord]
            oldlab <- newlabels1
            for (j in seq_len(k1))
                newlabels1[oldlab == ord[j]] <- j
            newlabels1 <- rep(10*l1,l) + newlabels1
        }
    }
    if (k1 == 1) {
        newmedoids1 <- medoid1
        newlabels1 <- rep(10*l1,p1)
        newclussizes <- p1
    }
    for (a in seq_len(length(newmedoids1))) {
        if (sum(newmedoids1[a] == id1) == 0)
            warning("Problem with new medoids after splitting cluster")
    }
    list(k1,newmedoids1,newlabels1,newclussizes)
}



mssnextlevel <- function(data ,prevlevel, dmat, kmax = 9,
                         khigh = 9, within = "med", between = "med") {
    n <- length(data[1,])
    p <- length(data[,1])

    id <- seq_len(p)
    k <- prevlevel[[1]]
    medoids <- prevlevel[[2]]
    labels <- prevlevel[[4]]
    newk <- 0
    newlabels <- newmedoids <- newclussizes <- NULL
    count <- 1
    ordlabels <- sort(unique(labels))
    if (length(ordlabels) != k)
        warning("Number of unique labels not equal
                number of clusters in mssnextlevel()")
    if (sum(is.na(medoids))) {
        warning("Missing values in medoid vector in nextlevel()")
        medoids[is.na(medoids)] <- FALSE
    }
    if (length(unique(medoids)) < k && sum(medoids))
        warning("Medoids not unique in mssnextlevel()")
    checkmeans <- FALSE
    if (length(medoids) == 1 && !medoids) {
        warning("No medoids provided in mssnxtlevel()")
        usemean <- TRUE
    }
    else{
        if (sum(medoids > 1) == k)
            usemean <- FALSE
        else
            checkmeans <- TRUE
    }
    for (j in (seq_len(k))) {
        clust1 <- data[labels == ordlabels[j],]
        id1 <- id[labels == ordlabels[j]]
        if (length(id1) > 1)
            clust1 <- as.matrix(clust1)
        l1 <- ordlabels[j]
        right <- (j < k)
        medoid1 <- ifelse(is.na(medoids[j]),0,medoids[j])
        if (j < k)
            medoid2 <- medoids[j + 1]
        else
            medoid2 <- medoids[j - 1]
        if (length(id1) > 1) {
            med2dist <- rowMeans(as(dmat[labels == ordlabels[j],
                                         labels == labels[medoid2],
                                         drop = FALSE],"matrix"))
        }else{
            med2dist <- mean(dmat[labels == ordlabels[j],
                                  labels == labels[medoid2]])
        }

        splitobj <- msssplitcluster(clust1,l1,id1,medoid1,med2dist,right,
                                    dmat[labels == l1,labels == l1],
                                    kmax,khigh,within,between)
        newlabels[labels == ordlabels[j]] <- splitobj[[3]]
        k1 <- splitobj[[1]]
        start <- count
        end <- count + k1 - 1
        newmedoids[start:end] <- splitobj[[2]]
        newclussizes[start:end] <- splitobj[[4]]
        count <- count + k1
    }
    count <- newk <- count - 1
    newmedoids <- newmedoids[seq_len(newk)]
    newclussizes <- newclussizes[seq_len(newk)]
    final <- 0
    if (count == k)
        final <- 1
    if (max(newclussizes) == 3)
        final <- 1
    list(newk,newmedoids,newclussizes,newlabels,final,
         rbind(prevlevel[[6]],cbind(sort(unique(newlabels)),newmedoids)))
}

orderelements  <-  function(level, data, rel = "own",
                            d = "cosangle", dmat = NULL) {

    idn <- seq_len(length(data[,1]))
    k <- level[[1]]
    labels <- level[[4]]
    medoids <- level[[2]]
    clussizes <- level[[3]]
    ord <- order(labels)
    idnord <- idn[ord]
    subdataord <- data[ord,]

    if (is.null(dmat))
        dmat <- distancematrix(data,d)
    distord <- dmat[ord,]

    labelsord <- labels[ord]
    count <- 1
    for (j in (seq_len(k))) {
        start <- count
        end <- count + clussizes[j] - 1
        if (clussizes[j] > 2) {
            tempid <- idnord[start:end]
            if (rel == "co") {
                distj <- distord[,ord][start:end,start:end]
                idnord[start:end] <- tempid[improveordering(distj)]
            }
            else{
                if (rel == "neighbor") {
                    if (j < k)
                        mednext <- medoids[j + 1]
                    else
                        mednext <- medoids[j - 1]
                }else
                    mednext <- medoids[j]

                dmednext <- distord[start:end,mednext]

                if (rel == "neighbor") {
                    if (j < k)
                        #ordtemp <- order(dmednext, decreasing = TRUE)
                        ordtemp <- rev(order(dmednext))
                    else
                        ordtemp <- order(dmednext)
                }
                else
                    ordtemp <- order(dmednext)
                idnord[start:end] <- tempid[ordtemp]
            }
        }
        else
            idnord[start:end] <- idnord[start:end]
        count <- count + clussizes[j]
    }
    list(data[idnord,],idnord)
}



mssinitlevel <- function(data, kmax=9, khigh=9, d="cosangle", dmat=NULL,
                         within="med", between="med", ord="co",
                         verbose=FALSE)
{
    #print("mssinitlevel")
    if (!is.matrix(data))
        stop("First arg to mssinitlevel() must be a matrix")

    p <- length(data[,1])

    if (is.matrix(dmat))
        dmat <- as.hdist(dmat)

    else
        dmat <- distancematrix(data,d = d)

    if (methods::slot(dmat, "Size") != p)
        stop("Data and distance matrix dimensions
             do not agree in mssinitlevel()")

    m <- msscheck(dmat,kmax,khigh,within,between)
    if (m[1] == 1) {
        if (verbose)
            cat("No strong evidence for clusters in the first level -
          \n continuing to split root node anyway. \n")
        m <- msscheck(dmat,kmax,khigh,within,between,force = TRUE)
    }

    pamobj <- pam(methods::slot(dmat, "Data"), m[1], diss = TRUE)
    rowmedoids <- pamobj$medoids
    final <- ifelse(max(pamobj$clusinfo[,1]) < 3,1,0)
    if (m[1] > 2) {
        medoidsdata <- as.matrix(data[rowmedoids,])
        l <- length(rowmedoids)
        medoidsdist <- dmat[rowmedoids,rowmedoids]
        if (ord == "co")
            medoidsord <- improveordering(medoidsdist)
        if (ord == "clust") {
            mpamobj <- pam(methods::slot(medoidsdist, "Size"), 2, diss = TRUE)
            labelsmed <- mpamobj$clustering
            medmed <- mpamobj$medoids
            clussizes <- mpamobj$clusinfo[,1]

            prevlevel <- mssnextlevel(medoidsdata,list(2,medmed,clussizes,
                                                       labelsmed,0,
                                                       cbind(c(1,2),medmed)),
                                      dmat = medoidsdist,kmax,khigh,within,
                                      between)
            final <- prevlevel[[5]]
            if (final == 0) {
                depth <- (l - 1)
                for (j in seq_len(depth)) {
                    if (final == 0) {
                        clustnext <- mssnextlevel(medoidsdata,prevlevel,
                                                  dmat = medoidsdist,kmax,
                                                  khigh,within,between)
                        final <- clustnext[[5]]
                    }
                    if (final == 1) {
                        j <- depth
                        prevlevel <- clustnext
                    }
                }
            }
            medoidsord <- orderelements(prevlevel,medoidsdata,rel = "neighbor",
                                        d = d,dmat = medoidsdist)[[2]]
        }
        k <- m[1]
        rowmedoids <- rowmedoids[medoidsord]
        labels <- lab2 <- pamobj$clustering
        for (j in (seq_len(k)))
            lab2[labels == medoidsord[j]] <- j
        output <- list(k,rowmedoids,pamobj$clusinfo[,1][medoidsord],lab2,
                       final,cbind(seq_len(k),rowmedoids))
    }
    else
        output <- list(2,pamobj$medoids,pamobj$clusinfo[,1],
                       pamobj$clustering,final,cbind(seq_len(2),pamobj$medoids))
    return(output)
}


paircoll <- function(i,j,data,level,d="cosangle",
                     dmat=NULL, newmed="medsil") {
    p <- length(data[,1])
    k <- level[[1]]
    labels <- level[[4]]
    medoids <- level[[2]]
    clussizes <- level[[3]]
    N <- length(level[[6]][,1])
    block <- level[[6]][(N - k + 1):N,]
    if (N == k)
        prevblock <- NULL
    else
        prevblock <- level[[6]][seq_len(N - k),]
    if (max(i,j) > k)
        stop("The clusters to collapse do not exist in paircoll()")
    labeli <- labels[medoids[i]]
    labelj <- labels[medoids[j]]
    oldlabels <- labels
    labels[labels == labelj] <- labeli
    trunclabels <- trunc(oldlabels/10)
    labelparents <- unique(trunclabels)
    parenti <- order(labelparents)[labelparents == trunc(labeli/10)]
    parentj <- order(labelparents)[labelparents == trunc(labelj/10)]
    if (newmed == "nn")
        fakemed <- (data[medoids[i],]*clussizes[i] +
                        data[medoids[j],]*clussizes[j])/(clussizes[i] +
                                                             clussizes[j])
    if (newmed == "uwnn")
        fakemed <- (data[medoids[i],] + data[medoids[j],])/2
    if (newmed == "nn" || newmed == "uwnn") {
        rowsub <- (seq_len(p))[labels == labeli]
        distsfm <- distancevector(data[rowsub,],as.vector(fakemed),d)
        medoids[i] <- rowsub[order(distsfm)[1]]
    }
    else{
        #colldist <- as.matrix(dmat[labels == labeli,labels == labeli])
        colldist <- as(dmat[labels == labeli,labels == labeli],"matrix")

        rowsub <- (seq_len(p))[labels == labeli]
        if (newmed == "center") {
            sumdist <- rowSums(colldist)
        }
        if (newmed == "medsil") {
            othermed <- medoids[-c(i,j)]
            collp <- length(labels[labels == labeli])
            othern <- length(othermed)
            if (othern == 0)
                stop("Not enough medoids to use newmed='medsil' in paircoll()")
            if (is(dmat,"hdist")) {
                if (othern == 1) {
                    #otherdist <- cbind(dmat[labels == labeli,othermed])
                    otherdist <- dmat[labels == labeli,othermed]
                }else{
                    #otherdist <- rbind(dmat[labels == labeli,othermed])
                    otherdist <- dmat[labels == labeli,othermed]
                }
            }

            if (othern == 1)
                b <- otherdist
            else
                b <- apply(otherdist,1,min)
            b <- matrix(b,nrow=collp,ncol=collp)
            diag(b) <- 0
            b <- abs(b-colldist)/pmax(colldist,b)
            sumdist <- rowSums(b)
            medoids[i] <- rowsub[order(sumdist,decreasing=TRUE) == 1]
            #medoids[i] <- rowsub[rev(order(sumdist)) == 1]
        }
    }
    k <- k-1
    clussizes[i] <- clussizes[i] + clussizes[j]
    block[i,2] <- medoids[i]
    if (j<=k) {
        for (l in (j:k)) {
            medoids[l] <- medoids[l + 1]
            clussizes[l] <- clussizes[l + 1]
            block[l,] <- block[l + 1,]
        }
    }
    medoids <- medoids[seq_len(k)]
    clussizes <- clussizes[seq_len(k)]
    block <- block[seq_len(k),]
    if (labels[labels == labeli][1]/10 == trunclabels[labels == labeli][1]) {
        labels[labels == labeli] <- labels[labels == labeli] + 1
        labels[labels == labelj] <- labels[labels == labelj] + 1
        block[i,1] <- block[i,1] + 1
    }
    return(list(k,medoids,clussizes,labels,level[[5]],rbind(prevblock,block)))
}


collap <- function(data,level,d = "cosangle", dmat = NULL, newmed = "medsil") {
    k <- level[[1]]
    if (k < 3) {
        warning("Not enough medoids to use newmed='medsil' in collap() - \n
            using newmed='nn' instead \n")
        newmed <- "nn"
    }
    medoids <- level[[2]]
    clussizes <- level[[3]]
    if (sum(is.na(clussizes)))
        warning("NA in clussizes")
    medoidsdata <- data[medoids,]
    if (sum(is.na(medoidsdata)) > 0)
        warning("Missing value(s) in medoidsdata in collap()")

    distmed <- dmat[medoids,medoids]
    distv <- methods::slot(distmed, "Data")

    indexmin <- order(distv)[1]
    best <- vectmatrix(indexmin,k)
    clustfinal <- paircoll(best[1],best[2],data,level,d,dmat,newmed)
    return(clustfinal)
}

msscollap <- function(data,level,khigh,d="cosangle", dmat = NULL,
                      newmed = "medsil", within = "med",
                      between = "med", impr = 0) {

    newk <- level[[1]]
    mss1 <- labelstomss(level[[4]],dmat,khigh,within,between)
    maxncoll <- max(0,newk - 2)
    ncoll <- 0
    coll <- 1
    if (newk <= 2)
        coll <- 0
    while ((coll == 1) && (ncoll <= maxncoll)) {
        levelc <- collap(data,level,d,dmat,newmed)
        mss2 <- labelstomss(levelc[[4]],dmat,khigh,within,between)
        if (mss1 == 0)
            r <- 0
        else
            r <- (mss1 - mss2)/mss1
        if (r < impr)
            coll <- 0
        else{
            mss1 <- mss2
            level <- levelc
            ncoll <- ncoll + 1
        }
    }
    return(level)
}

mssmulticollap <- function(data,level,khigh,d="cosangle",
                           dmat = NULL, newmed = "medsil",
                           within = "med", between = "med", impr = 0) {
    if (!is.matrix(data))
        stop("First arg to mssmulticollap() must be a matrix")

    medoids <- level[[2]]
    medoidsdata <- data[medoids,]
    if (sum(is.na(medoidsdata)) > 0)
        warning("Missing value(s) in medoidsdata in mssmulticollap()")

    distmed <- dmat[medoids,medoids]
    k <- level[[1]]
    ord <- order(methods::slot(distmed, "Data"))

    mss1 <- labelstomss(level[[4]],dmat,khigh,within,between)
    maxncoll <- max(0,k*(k - 1)/2)
    ncoll <- 0
    i <- 1
    while (i <= maxncoll) {
        clusts <- vectmatrix(ord[i],k)
        if (k < 3) {
            if (newmed == "medsil")
                warning("Can't use newmed=medsil with less than 3 clusters. \n
                        Substituting newmed=nn")
            newmed <- "nn"
        }
        levelc <- paircoll(clusts[1],clusts[2],data,level,d,dmat,newmed)
        mss2 <- labelstomss(levelc[[4]],dmat,khigh,within,between)
        if (mss1 == 0)
            r <- 0
        else
            r <- (mss1 - mss2)/mss1
        if (r >= impr) {
            mss1 <- mss2
            level <- levelc
            ncoll <- ncoll + 1
            k <- level[[1]]
            maxncoll <- max(0,k*(k - 1)/2)
            i <- 0
            medoids <- level[[2]]
            medoidsdata <- data[medoids,]
            if (sum(is.na(medoidsdata)) > 0)
                warning("Missing value(s) in medoidsdata in mssmulticollap()")

            distmed <- dmat[medoids,medoids]
            ord <- order(methods::slot(distmed, "Data"))
        }
        i <- i + 1
    }
    return(level)
}


mssrundown <- function(data, K=16, kmax=9, khigh=9, d="cosangle",
                       dmat=NULL, initord="co", coll="seq",
                       newmed="medsil", stop=TRUE,
                       finish=FALSE, within="med",
                       between="med",impr=0, verbose=FALSE)
{
    #print("mssrundown")
    if (!is.matrix(data))
        stop("First arg to mssrundown() must be a matrix")

    bestlevel <- level <- mssinitlevel(data, kmax, khigh,
                                       d, dmat, within,
                                       between, initord, verbose)
    bestmss <- mss <- labelstomss(level[[4]],dmat,khigh,within,between)
    bestl <- l <- 1
    ind <- 0
    if (verbose)
        cat("Searching for main clusters... \n")
    if (level[[5]] == 1)
        return(level)
    while((l<=K) && (ind == 0)) {
        if (verbose) cat("Level ",l,"\n")
        if (coll == "seq")
            levelc <- msscollap(data,level,khigh,d,dmat,
                                newmed,within,between,impr)
        if (coll == "all")
            levelc <- mssmulticollap(data,level,khigh,d,
                                     dmat,newmed,within,between,impr)
        mss <- labelstomss(levelc[[4]],dmat,khigh,within,between)
        if (mss>=bestmss & stop == TRUE)
            ind <- 1
        else{
            if (mss<bestmss) {
                bestlevel <- levelc
                bestmss <- mss
                bestl <- l
            }
        }
        l <- l+1
        if (l<=K) {
            level <- mssnextlevel(data,levelc,dmat,kmax,khigh,within,between)
            if (finish == TRUE) {
                if (sum(trunc(level[[4]]/10)*10 == level[[4]]) ==
                    length(level[[4]]) & l<=K) {
                    ind <- 1
                    bestlevel <- levelc
                    bestmss <- mss
                    bestl <- (l-1)
                }
            }
        }
    }
    if (verbose)
        cat("Identified", bestlevel[[1]],
            " main clusters in level",
            bestl, "with MSS =",bestmss,"\n")
    return(bestlevel)
}

msscomplete <- function(level, data, K=16, khigh=9, d="cosangle",
                        dmat=NULL, within="med", between="med", verbose=FALSE)
{
    if (!is.matrix(data))
        stop("First arg to msscomplete() must be a matrix")

    count <- digits(level[[4]][1])
    if (verbose)
        cat("Running down without collapsing from Level",count,"\n")
    while ((max(level[[3]]) > 3) & (count < K)) {
        level <- newnextlevel(data,level,dmat,2,khigh)
        count <- count + 1
        if (verbose) cat("Level",count,"\n")
    }
    return(level)
}



newnextlevel <- function(data, prevlevel, dmat, klow=2, khigh=6) {
    #print("newnextlevel")
    if (!is.matrix(data))
        stop("First arg to newnextlevel() must be a matrix")
    p <- length(data[,1])
    n <- length(data[1,])
    id <- seq_len(p)
    k <- prevlevel[[1]]
    medoids <- prevlevel[[2]]
    labels <- prevlevel[[4]]
    newk <- 0
    newlabels <- newmedoids <- newclussizes <- NULL
    count <- 1
    ordlabels <- sort(unique(labels))
    if (length(ordlabels) != k)
        warning("Number of unique labels not equal
                number of clusters in newnextlevel()")
    if (sum(is.na(medoids))) {
        warning("Missing value(s) in medoid vector in newnextlevel()")
        medoids[is.na(medoids)] <- FALSE
    }
    if (length(unique(medoids)) < k && sum(medoids))
        warning("Medoids in newnextlevel() are not unique")
    checkmeans <- FALSE
    if (length(medoids) == 1 && !medoids) {
        warning("No medoids provided in newnextlevel()")
        usemean <- TRUE
    }
    else{
        if (sum(medoids > 1) == k)
            usemean <- FALSE
        else
            checkmeans <- TRUE
    }
    for (j in (seq_len(k))) {
        clust1 <- data[labels == ordlabels[j],]
        id1 <- id[labels == ordlabels[j]]
        if (length(id1) > 1) {
            kmax <- min(c(khigh,dim(clust1)[1] - 1))
            clust1 <- as.matrix(clust1)
        }
        l1 <- ordlabels[j]
        right <- (j < k)
        medoid1 <- ifelse(is.na(medoids[j]),0,medoids[j])
        if (j < k)
            medoid2 <- medoids[j + 1]
        else
            medoid2 <- medoids[j - 1]
        if (length(id1) > 1)
            med2dist <- rowMeans(as(dmat[labels == ordlabels[j],
                                         labels == labels[medoid2],
                                         drop = FALSE],"matrix"))
        else
            med2dist <- mean(dmat[labels == ordlabels[j],
                                  labels == labels[medoid2]])

        splitobj <- newsplitcluster(clust1,l1,id1,klow,kmax,medoid1,med2dist,
                                    right,dmat[labels == l1,labels == l1])
        newlabels[labels == ordlabels[j]] <- splitobj[[3]]
        k1 <- splitobj[[1]]
        start <- count
        end <- count + k1 - 1
        newmedoids[start:end] <- splitobj[[2]]
        newclussizes[start:end] <- splitobj[[4]]
        count <- count + k1 - 1 + 1
    }
    count <- count - 1
    newk <- count
    newmedoids <- newmedoids[seq_len(newk)]
    newclussizes <- newclussizes[seq_len(newk)]
    final <- 0
    if (count == k)
        final <- 1
    if (max(newclussizes) == 3)
        final <- 1
    return(list(newk,newmedoids,newclussizes,newlabels,
                final,rbind(rbind(prevlevel[[6]],
                                  cbind(sort(unique(newlabels)),newmedoids)))))
}

newsplitcluster <- function(clust1, l1, id1, klow=2, khigh=2, medoid1,
                            med2dist, right, dist1) {
    #print("newsplitcluter")
    if (!medoid1)
        warning("Medoid missing - continue to split cluster")
    else{
        if (sum(medoid1 == id1) == 0 & medoid1)
            warning("Medoid not in cluster - continue to split cluster")
    }
    if (is.matrix(clust1)) {
        p1 <- length(clust1[,1])
        n <- length(clust1[1,])
    }
    else p1 <- 1
    if (p1 < 3) {
        k1 <- 1
        newmedoids1 <- medoid1
        newlabels1 <- rep(10*l1,p1)
        newclussizes <- p1
    }
    else{
        l <- length(clust1[,1])
        dissvec <- methods::slot(dist1, "Data")

        kmax <- min(p1 - 1, khigh)
        a <- rep(0,(kmax - klow + 2))
        best <- 2
        for (j in (klow:kmax)) {
            a[j] <- pam(dissvec,j,diss = TRUE)$silinfo$avg.width
            if (a[j] > a[best]) best <- j
        }
        k1 <- best
        pamobj <- pam(dissvec, k1, diss = TRUE)
        newclussizes <- pamobj$clusinfo[,1]
        newmedoids1 <- id1[pamobj$medoids]
        newlabels1 <- pamobj$clustering
        distnewmedoids <- NULL
        for (j in seq_len(k1))
            distnewmedoids[j] <-
            mean(med2dist[newlabels1 ==
                              newlabels1[pamobj$medoids[j]]])
        if (right == 1)
            ord <- order(distnewmedoids, decreasing = TRUE)
        #ord <- rev(order(distnewmedoids))
        else
            ord <- order(distnewmedoids)
        newmedoids1 <- newmedoids1[ord]
        newclussizes <- newclussizes[ord]
        oldlab <- newlabels1
        for (j in seq_len(k1))
            newlabels1[oldlab == ord[j]] <- j
        newlabels1 <- rep(10*l1,l) + newlabels1
    }
    for (a in seq_len(length(newmedoids1))) {
        if (sum(newmedoids1[a] == id1) == 0)
            warning("Problem with new medoids after splitting cluster")
    }
    return(list(k1,newmedoids1,newlabels1,newclussizes))
}



hopach <- function(data, dmat=NULL, d="cosangle", clusters="best", K=15,
                   kmax=9, khigh=9, coll="seq", newmed="medsil",
                   mss="med", impr=0, initord="co", ord="own", verbose=FALSE) {
    # if (inherits(data,"ExpressionSet"))
    #   data  <-  exprs(data)
    data <- as.matrix(data)

    # Convert to hdist immediately #
    if ( is.null(dmat) ) {
        dmat <- distancematrix(data, d)
    }else if ( is.matrix(dmat) ) {
        dmat <- as(dmat, "hdist")
    }else if ( "dist" %in% is(dmat)) {
        dmat <- hdist(Data = as.numeric(dmat), Size = attr(dmat, "Size"),
                      Labels = seq_len(attr(dmat, "Size")),
                      Call = as.character(attr(dmat, "call"))[3])
    }else if (!is.hdist(dmat)) {
        stop("Distance matrix could not be created into hdist object.")
    }

    if (K > 15) {
        K <- 15
        warning("K set to 15 - can't do more than 15 splits")
    }
    if (K < 1) {
        K <- 1
        warning("K set to 1 - can't do less than 1 level")
    }
    if (clusters != "none") {
        cuttree <- mssrundown(data, K, kmax, khigh, d, dmat,
                              initord, coll, newmed,
                              stop = (clusters == "greedy"), finish = TRUE,
                              within = mss,
                              between = mss, impr, verbose)
        if (cuttree[[1]] > 1)
            cutord <- orderelements(cuttree, data, rel = ord, d, dmat)[[2]]
        else
            cutord <- NULL
        out1 <- list(k = cuttree[[1]], medoids = cuttree[[2]],
                     sizes = cuttree[[3]], labels = cuttree[[4]],
                     order = cutord)
        finaltree <- msscomplete(cuttree, data, K, khigh, d,
                                 dmat, within = mss, between = mss,
                                 verbose)
    }
    else{
        out1 <- NULL
        finaltree <- msscomplete(mssinitlevel(as.matrix(data),
                                              kmax, khigh, d, dmat,
                                              within = mss, between = mss,
                                              initord), data, K, khigh,
                                 d, dmat, within = mss,
                                 between = mss, verbose)
    }
    dimnames(finaltree[[6]]) <- list(NULL,c("label","medoid"))
    out2 <- list(labels = finaltree[[4]],
                 order = orderelements(finaltree, data, rel = ord, d,
                                       dmat)[[2]], medoids = finaltree[[6]])
    return(list(clustering = out1, final = out2, call = match.call(),
                metric = d))
}

distancematrix  <-  function(X, d, na.rm = TRUE)
{
    X  <-  as.matrix(X)
    if (d == "euclid")
        return(hopach::disseuclid(X, na.rm))
    if (d == "cor")
        return(hopach::disscor(X, na.rm))
    if (d == "abscor")
        return(hopach::dissabscor(X, na.rm))
    if (d == "cosangle")
        return(hopach::disscosangle(X, na.rm))
    if (d == "abscosangle")
        return(hopach::dissabscosangle(X, na.rm))
    stop("Distance metric ", d, " not available")
}

