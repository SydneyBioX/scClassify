#' @importFrom graphics plot abline
#' @importFrom stats uniroot qnorm dnorm coef
#' @importFrom mixtools normalmixEM

KNNcor <- function(corMat,
                   subLevelModel,
                   cutoff_method = c("dynamic", "static"),
                   k = 10,
                   prob_threshold = 0.8,
                   cor_threshold_static = 0.5,
                   cor_threshold_high = 0.7,
                   topLevel = F,
                   verbose = T) {


  cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"), several.ok = FALSE)


  # If we use static correlation cutoff
  if (cutoff_method == "static") {
    if (verbose) {
      print("Using static correlation cutoff...")
    }

    # get the KNN labels
    topKNN <- apply(corMat, 2, function(x) subLevelModel$y[order(x, decreasing = T)][1:k])

    # get the KNN correlation
    topKNN_cor <- apply(corMat, 2, function(x) x[order(x, decreasing = T)][1:k])

    # If the KNN labels with correlation less than the threshood,
    # set the labels as -1
    topKNN <- ifelse(topKNN_cor >= cor_threshold_static, topKNN, -1)


    # Get the predicted results
    predRes <- apply(topKNN, 2, function(x) {
      tab <- table(x)/length(x)
      if (max(tab, na.rm = T) < prob_threshold) {
        0
      }else{
        if (names(tab)[which(tab == max(tab, na.rm = TRUE))] == "-1") {
          0
        }else{
          names(tab)[which(tab == max(tab, na.rm = TRUE))]
        }
      }
    })

  }


  # If we use dynamic correlation cutoff based on mixture model
  if (cutoff_method == "dynamic") {

    if (verbose) {
      print("Using static correlation cutoff...")
    }


    # If this is the top level
    if (topLevel) {
      # if(i==2){

      # cat("The first layer correlation cutoff: \n")


      # We select different cutoffs for each cell type
      unique_y <- levels(subLevelModel$y)
      cor_threshold_tmp <- c()

      for (l in 1:length(unique_y)) {
        print(paste("y", l))
        # Fitting the mixture model if it is not unimodal distribution.
        corMat_vec <- corMat_vec[sort(sample(length(corMat_vec), max(round(length(corMat_vec)*0.001), min(200000, length(corMat_vec)))))]
        corMat_vec <- corMat_vec[corMat_vec != min(corMat)]
        suppressMessages(dip_test <- diptest::dip.test(corMat_vec, B = 10000))
        if (dip_test$p.value <= 0.01) {
          quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                               k = length(unique_y), maxit = 2000,
                                                               mu = c(-0.5, rep(0.5, length(unique_y) - 2),  1),
                                                               ECM = TRUE, verb = verbose),
                                         silent = T))

          if (class(mixmdl) != "try-error" ) {
            if (suppressWarnings(min(unique(mixmdl$rho)) == 0)) {
              quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                   k = length(unique_y), maxit = 2000,
                                                                   ECM = TRUE, verb = verbose),
                                             silent = T))
            }
          }else{
            quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                 k = length(unique_y), maxit = 2000,
                                                                 ECM = TRUE, verb = verbose),
                                           silent = T))
          }
        }

        # Caculate the threshold for this branch
        if (dip_test$p.value > 0.01) {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else if (class(mixmdl) == "try-error") {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else{
          # plot(mixmdl,which = 2)
          t_G2 <- getThreshold(mixmdl, verbose = verbose)
          cor_threshold_tmp = c(cor_threshold_tmp, t_G2)
        }



        # if (verbose) {
        #   cat("Correlation threshold: \n")
        #   print(cor_threshold_tmp)
        # }

      }
      names(cor_threshold_tmp) <- unique_y


      topKNN <- apply(corMat, 2, function(x) subLevelModel$y[order(x, decreasing = T)][1:k])
      topKNN_cor <- apply(corMat, 2, function(x) x[order(x, decreasing = T)][1:k])

      # get the threshold
      topKNN_threshold <- apply(corMat, 2,
                                function(x) cor_threshold_tmp[as.numeric((subLevelModel$y[order(x, decreasing = T)][1:k]))])

      topKNN <- ifelse(topKNN_cor >= topKNN_threshold, topKNN, -1)

      predRes <- apply(topKNN, 2, function(x){
        # x <- stats::na.omit(x)

        tab <- table(x)/length(x)
        if (max(tab, na.rm = T) < prob_threshold) {
          0
        }else{
          if (names(tab)[which(tab == max(tab, na.rm = TRUE))] == "-1") {
            0
          }else{
            names(tab)[which(tab == max(tab, na.rm = TRUE))]
          }

        }
      })


    }else{

      unique_y <- levels(subLevelModel$y)
      cor_threshold_tmp <- c()

      for (l in 1:length(unique_y)) {


        corMat_vec <- c(corMat[subLevelModel$y == unique_y[l],])

        # Fitting the mixture model if it is not unimodal distribution.
        corMat_vec <- corMat_vec[sort(sample(length(corMat_vec), max(round(length(corMat_vec)*0.001), min(200000, length(corMat_vec)))))]
        corMat_vec <- corMat_vec[corMat_vec != min(corMat)]
        suppressMessages(dip_test <- diptest::dip.test(corMat_vec, B = 10000))
        if (dip_test$p.value <= 0.01) {
          quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                               k = length(unique_y), maxit = 2000,
                                                               mu = c(-0.5, rep(0.5, length(unique_y) - 2),  1),
                                                               ECM = TRUE, verb = verbose),
                                         silent = T))

          if (class(mixmdl) != "try-error" ) {
            if (suppressWarnings(min(unique(mixmdl$rho)) == 0)) {
              quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                   k = length(unique_y), maxit = 2000,
                                                                   ECM = TRUE, verb = verbose),
                                             silent = T))
            }
          }else{
            quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                 k = length(unique_y), maxit = 2000,
                                                                 ECM = TRUE, verb = verbose),
                                           silent = T))
          }
        }

        # Caculate the threshold for this branch
        if (dip_test$p.value > 0.01) {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else if (class(mixmdl) == "try-error") {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else{
          # plot(mixmdl,which = 2)
          t_G2 <- getThreshold(mixmdl, verbose = verbose)
          cor_threshold_tmp = c(cor_threshold_tmp, t_G2)
        }


      }
      names(cor_threshold_tmp) <- unique_y

      topKNN <- apply(corMat, 2, function(x) subLevelModel$y[order(x, decreasing = T)][1:k])
      topKNN_cor <- apply(corMat, 2, function(x) x[order(x, decreasing = T)][1:k])

      # get the threshold
      # the threshold_tmp is order based on the levels of factor. no need to change the character
      topKNN_threshold <- apply(corMat, 2,
                                function(x) cor_threshold_tmp[as.numeric((subLevelModel$y[order(x, decreasing = T)][1:k]))])

      topKNN <- ifelse(topKNN_cor >= topKNN_threshold, topKNN, -1)


      predRes <- apply(topKNN, 2, function(x){
        # x <- stats::na.omit(x)
        tab <- table(x)/length(x)
        if (max(tab, na.rm = T) < prob_threshold) {
          0
        }else{
          if (names(tab)[which(tab == max(tab, na.rm = TRUE))] == "-1") {
            0
          }else{
            names(tab)[which(tab == max(tab, na.rm = TRUE))]
          }

        }
      })

    }


  }

  return(list(predRes = predRes))
}

WKNNcor <- function(corMat,
                    subLevelModel,
                    cutoff_method = c("dynamic", "static"),
                    k = 10,
                    prob_threshold = 0.8,
                    cor_threshold_static = 0.5,
                    cor_threshold_high = 0.7,
                    topLevel = F,
                    verbose = T){

  cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"), several.ok = FALSE)

  if (cutoff_method == "static") {

    if (verbose) {
      print("Using static correlation cutoff...")
    }


    topKNN <- apply(corMat, 2, function(x) subLevelModel$y[order(x, decreasing = T)][1:k])
    topKNN_cor <- apply(corMat, 2, function(x) x[order(x, decreasing = T)][1:k])
    topKNN <- ifelse(topKNN_cor >= cor_threshold_static, topKNN, -1)

    topKNN_weight <- apply(topKNN_cor, 2, function(x){
      x <- stats::na.omit(x) # add in 20190613
      if (x[1] == x[length(x)]) {
        rep(1, length(x))
      } else{
        (x - x[length(x)])/(x[1] - x[length(x)])
      }
    })



    predRes <- sapply(1:ncol(corMat),
                      function(i){
                        vote <- stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum)
                        maxIdx <- which.max(vote$x)
                        if (max(vote$x/sum(vote$x)) < prob_threshold) {
                          0
                        }else{
                          if (vote$Group.1[maxIdx] == "-1") {
                            0
                          }else{
                            vote$Group.1[maxIdx]
                          }
                        }


                      })

    names(predRes) <- colnames(corMat)


  }

  if (cutoff_method == "dynamic") {

    if (verbose) {
      print("Using dynamic correlation cutoff...")
    }


    # if(i==2){
    if (topLevel) {
      # cat("The first layer correlation cutoff: \n")



      unique_y <- levels(subLevelModel$y)
      cor_threshold_tmp <- c()
      for (l in 1:length(unique_y)) {
        # print(paste("y", l))

        # Fitting the mixture model if it is not unimodal distribution.
        corMat_vec <- corMat_vec[sort(sample(length(corMat_vec), max(round(length(corMat_vec)*0.001), min(200000, length(corMat_vec)))))]
        corMat_vec <- corMat_vec[corMat_vec != min(corMat)]
        suppressMessages(dip_test <- diptest::dip.test(corMat_vec, B = 10000))
        if (dip_test$p.value <= 0.01) {
          quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                               k = length(unique_y), maxit = 2000,
                                                               mu = c(-0.5, rep(0.5, length(unique_y) - 2),  1),
                                                               ECM = TRUE, verb = verbose),
                                         silent = T))

          if (class(mixmdl) != "try-error" ) {
            if (suppressWarnings(min(unique(mixmdl$rho)) == 0)) {
              quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                   k = length(unique_y), maxit = 2000,
                                                                   ECM = TRUE, verb = verbose),
                                             silent = T))
            }
          }else{
            quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                 k = length(unique_y), maxit = 2000,
                                                                 ECM = TRUE, verb = verbose),
                                           silent = T))
          }
        }

        # Caculate the threshold for this branch
        if (dip_test$p.value > 0.01) {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else if (class(mixmdl) == "try-error") {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else{
          # plot(mixmdl,which = 2)
          t_G2 <- getThreshold(mixmdl, verbose = verbose)
          cor_threshold_tmp = c(cor_threshold_tmp, t_G2)
        }


        if (verbose) {
          graphics::hist(corMat_vec, breaks = 100)
          graphics::abline(v = cor_threshold_tmp, col = "red", lwd = 3)

        }
      }


      names(cor_threshold_tmp) <- unique_y

      # if (verbose) {
      #   cat("Correlation threshold: \n")
      #   print(cor_threshold_tmp)
      # }

      ### calcualte the Weighted KNN
      topKNN <- apply(corMat, 2, function(x)subLevelModel$y[order(x, decreasing = T)][1:k])
      topKNN_cor <- apply(corMat, 2, function(x)x[order(x, decreasing = T)][1:k])
      topKNN_threshold <- apply(corMat, 2, function(x)cor_threshold_tmp[as.numeric(subLevelModel$y[order(x, decreasing = T)][1:k])])
      topKNN <- ifelse(topKNN_cor >= topKNN_threshold, topKNN, min(corMat))

      topKNN_weight <- apply(topKNN_cor, 2, function(x){
        x <- stats::na.omit(x) # add 20190613
        if (x[1] == x[length(x)]) {
          rep(1, length(x))
        }else{
          (x - x[length(x)])/(x[1] - x[length(x)])
        }
      })

      # print(topKNN_weight[,1:6])

      # if (verbose) {
      #   votes <- sapply(1:ncol(corMat),
      #                   function(i){
      #                     vote <- (stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum))
      #                     vote$x/sum(vote$x)
      #
      #
      #                   })
      #   print(head(votes))
      #   cat("Summary of the votes \n")
      #   print(summary(unlist(votes)))
      #   graphics::hist(unlist(votes), main = "Votes", breaks = 100)
      # }



      predRes <- sapply(1:ncol(corMat),
                        function(i){
                          vote <- stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum)
                          maxIdx <- which.max(vote$x)
                          if (max(vote$x/sum(vote$x)) < prob_threshold) {
                            0
                          }else{
                            if (vote$Group.1[maxIdx] == "-1") {
                              0
                            }else{
                              vote$Group.1[maxIdx]
                            }
                          }
                        })
      names(predRes) <- colnames(corMat)



    }else{

      unique_y <- levels(subLevelModel$y)
      cor_threshold_tmp <- c()
      for (l in 1:length(unique_y)) {

        corMat_vec <- c(as.matrix(corMat[subLevelModel$y == unique_y[l],]))

        # Fitting the mixture model if it is not unimodal distribution.
        corMat_vec <- corMat_vec[sort(sample(length(corMat_vec), max(round(length(corMat_vec)*0.001), min(200000, length(corMat_vec)))))]
        corMat_vec <- corMat_vec[corMat_vec != min(corMat)]
        suppressMessages(dip_test <- diptest::dip.test(corMat_vec, B = 10000))
        if (dip_test$p.value <= 0.01) {
          quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                               k = length(unique_y), maxit = 2000,
                                                               mu = c(-0.5, rep(0.5, length(unique_y) - 2),  1),
                                                               ECM = TRUE, verb = verbose),
                                         silent = T))

          if (class(mixmdl) != "try-error" ) {
            if (suppressWarnings(min(unique(mixmdl$rho)) == 0)) {
              quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                   k = length(unique_y), maxit = 2000,
                                                                   ECM = TRUE, verb = verbose),
                                             silent = T))
            }
          }else{
            quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                 k = length(unique_y), maxit = 2000,
                                                                 ECM = TRUE, verb = verbose),
                                           silent = T))
          }
        }

        # Caculate the threshold for this branch
        if (dip_test$p.value > 0.01) {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else if (class(mixmdl) == "try-error") {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else{
          # plot(mixmdl,which = 2)
          t_G2 <- getThreshold(mixmdl, verbose = verbose)
          cor_threshold_tmp = c(cor_threshold_tmp, t_G2)
        }

      }
      names(cor_threshold_tmp) <- unique_y

      # if (verbose) {
      #   cat("Correlation threshold: \n")
      #   print(cor_threshold_tmp)
      # }

      ### calcualte the Weighted KNN
      topKNN <- apply(corMat, 2, function(x)subLevelModel$y[order(x, decreasing = T)][1:k])
      topKNN_cor <- apply(corMat, 2, function(x)x[order(x, decreasing = T)][1:k])
      topKNN_threshold <- apply(corMat, 2, function(x)cor_threshold_tmp[as.numeric(subLevelModel$y[order(x, decreasing = T)][1:k])])
      topKNN <- ifelse(topKNN_cor >= topKNN_threshold, topKNN, -1)

      topKNN_weight <- apply(topKNN_cor, 2, function(x){
        x <- stats::na.omit(x) # add 20190613
        if (x[1] == x[length(x)]) {
          rep(1, length(x))
        }else{
          (x - x[length(x)])/(x[1] - x[length(x)])
        }
      })

      # if (verbose) {
      #   votes <- sapply(1:ncol(corMat),
      #                   function(i){
      #                     vote <- (stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum))
      #                     vote$x/sum(vote$x)
      #
      #
      #                   })
      #
      #   cat("Summary of the votes \n")
      #   print(summary(unlist(votes)))
      #   graphics::hist(unlist(votes), main = "Votes", breaks = 100)
      # }




      predRes <- sapply(1:ncol(corMat),
                        function(i){
                          vote <- stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum)
                          maxIdx <- which.max(vote$x)
                          if (max(vote$x/sum(vote$x)) < prob_threshold) {
                            0
                          }else{
                            if (vote$Group.1[maxIdx] == "-1") {
                              0
                            }else{
                              vote$Group.1[maxIdx]
                            }
                          }


                        })
      names(predRes) <- colnames(corMat)

    }


  }

  return(list(predRes = predRes))
}






DWKNNcor <- function(corMat,
                     subLevelModel,
                     cutoff_method = c("dynamic", "static"),
                     k = 10,
                     prob_threshold = 0.8,
                     cor_threshold_static = 0.5,
                     cor_threshold_high = 0.7,
                     topLevel = F,
                     verbose = T){

  cutoff_method <- match.arg(cutoff_method, c("dynamic", "static"), several.ok = FALSE)

  if (cutoff_method == "static") {

    if (verbose) {
      print("Using static correlation cutoff...")
    }


    topKNN <- apply(corMat, 2, function(x) subLevelModel$y[order(x, decreasing = T)][1:k])
    topKNN_cor <- apply(corMat, 2, function(x) x[order(x, decreasing = T)][1:k])
    topKNN <- ifelse(topKNN_cor >= cor_threshold_static, topKNN, -1)



    topKNN_weight <- apply(topKNN_cor, 2, function(x){
      x <- stats::na.omit(x) # add in 20190613
      if (x[1] == x[length(x)]) {
        rep(1, length(x))
      }else{
        (x - x[length(x)])*(2 - x[length(x)] - x[1])/((x[1] - x[length(x)])*(2 - x[length(x)] - x))
      }
    })


    predRes <- sapply(1:ncol(corMat),
                      function(i){
                        vote <- stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum)
                        maxIdx <- which.max(vote$x)
                        if (max(vote$x/sum(vote$x)) < prob_threshold) {
                          0
                        }else{
                          if (vote$Group.1[maxIdx] == "-1") {
                            0
                          }else{
                            vote$Group.1[maxIdx]
                          }
                        }


                      })

    names(predRes) <- colnames(corMat)


  }

  if (cutoff_method == "dynamic") {

    if (verbose) {
      print("Using dynamic correlation cutoff...")
    }


    # if(i==2){
    if (topLevel) {
      # cat("The first layer correlation cutoff: \n")



      unique_y <- levels(subLevelModel$y)
      cor_threshold_tmp <- c()
      for (l in 1:length(unique_y)) {
        # print(paste("y", l))

        # Fitting the mixture model if it is not unimodal distribution.
        corMat_vec <- corMat_vec[sort(sample(length(corMat_vec), max(round(length(corMat_vec)*0.001), min(200000, length(corMat_vec)))))]
        corMat_vec <- corMat_vec[corMat_vec != min(corMat)]
        suppressMessages(dip_test <- diptest::dip.test(corMat_vec, B = 10000))
        if (dip_test$p.value <= 0.01) {
          quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                               k = length(unique_y), maxit = 2000,
                                                               mu = c(-0.5, rep(0.5, length(unique_y) - 2),  1),
                                                               ECM = TRUE, verb = verbose),
                                         silent = T))

          if (class(mixmdl) != "try-error" ) {
            if (suppressWarnings(min(unique(mixmdl$rho)) == 0)) {
              quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                   k = length(unique_y), maxit = 2000,
                                                                   ECM = TRUE, verb = verbose),
                                             silent = T))
            }
          }else{
            quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                 k = length(unique_y), maxit = 2000,
                                                                 ECM = TRUE, verb = verbose),
                                           silent = T))
          }
        }

        # Caculate the threshold for this branch
        if (dip_test$p.value > 0.01) {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else if (class(mixmdl) == "try-error") {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else{
          # plot(mixmdl,which = 2)
          t_G2 <- getThreshold(mixmdl, verbose = verbose)
          cor_threshold_tmp = c(cor_threshold_tmp, t_G2)
        }
      }

      names(cor_threshold_tmp) <- unique_y

      # if (verbose) {
      #   cat("Correlation threshold: \n")
      #   print(cor_threshold_tmp)
      # }

      ### calcualte the Weighted KNN
      topKNN <- apply(corMat, 2, function(x)subLevelModel$y[order(x, decreasing = T)][1:k])
      topKNN_cor <- apply(corMat, 2, function(x)x[order(x, decreasing = T)][1:k])
      topKNN_threshold <- apply(corMat, 2, function(x)cor_threshold_tmp[as.numeric(subLevelModel$y[order(x, decreasing = T)][1:k])])
      topKNN <- ifelse(topKNN_cor >= topKNN_threshold, topKNN, min(corMat))

      topKNN_weight <- apply(topKNN_cor, 2, function(x){
        x <- stats::na.omit(x) # add in 20190613
        if (x[1] == x[length(x)]) {
          rep(1, length(x))
        }else{
          (x - x[length(x)])*(2 - x[length(x)] - x[1])/((x[1] - x[length(x)])*(2 - x[length(x)] - x))
        }
      })


      # if (verbose) {
      #   votes <- sapply(1:ncol(corMat),
      #                   function(i){
      #                     vote <- (stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum))
      #                     vote$x/sum(vote$x)
      #
      #
      #                   })
      #   # print(head(votes))
      #   cat("Summary of the votes \n")
      #   print(summary(unlist(votes)))
      #   graphics::hist(unlist(votes), main = "Votes", breaks = 100)
      # }



      predRes <- sapply(1:ncol(corMat),
                        function(i){
                          vote <- stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum)
                          maxIdx <- which.max(vote$x)
                          if (max(vote$x/sum(vote$x)) < prob_threshold) {
                            0
                          }else{
                            if (vote$Group.1[maxIdx] == "-1") {
                              0
                            }else{
                              vote$Group.1[maxIdx]
                            }
                          }
                        })
      names(predRes) <- colnames(corMat)



    }else{

      unique_y <- levels(subLevelModel$y)
      cor_threshold_tmp <- c()
      for (l in 1:length(unique_y)) {

        corMat_vec <- c(corMat[subLevelModel$y == unique_y[l],])



        # Fitting the mixture model if it is not unimodal distribution.
        corMat_vec <- corMat_vec[sort(sample(length(corMat_vec), max(round(length(corMat_vec)*0.001), min(200000, length(corMat_vec)))))]
        corMat_vec <- corMat_vec[corMat_vec != min(corMat)]
        suppressMessages(dip_test <- diptest::dip.test(corMat_vec, B = 10000))
        if (dip_test$p.value <= 0.01) {
          quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                               k = length(unique_y), maxit = 2000,
                                                               mu = c(-0.5, rep(0.5, length(unique_y) - 2),  1),
                                                               ECM = TRUE, verb = verbose),
                                         silent = T))

          if (class(mixmdl) != "try-error" ) {
            if (suppressWarnings(min(unique(mixmdl$rho)) == 0)) {
              quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                   k = length(unique_y), maxit = 2000,
                                                                   ECM = TRUE, verb = verbose),
                                             silent = T))
            }
          }else{
            quiet(mixmdl <- try(mixtools::normalmixEM(corMat_vec, fast = T, maxrestarts = 100,
                                                                 k = length(unique_y), maxit = 2000,
                                                                 ECM = TRUE, verb = verbose),
                                           silent = T))
          }
        }

        # Caculate the threshold for this branch
        if (dip_test$p.value > 0.01) {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else if (class(mixmdl) == "try-error") {
          cor_threshold_tmp = c(cor_threshold_tmp, 0)
        }else{
          # plot(mixmdl,which = 2)
          t_G2 <- getThreshold(mixmdl, verbose = verbose)
          cor_threshold_tmp = c(cor_threshold_tmp, t_G2)
        }


      }
      names(cor_threshold_tmp) <- unique_y

      # if (verbose) {
      #   cat("Correlation threshold: \n")
      #   print(cor_threshold_tmp)
      # }

      ### calcualte the Weighted KNN
      topKNN <- apply(corMat, 2, function(x)subLevelModel$y[order(x, decreasing = T)][1:k])
      topKNN_cor <- apply(corMat, 2, function(x)x[order(x, decreasing = T)][1:k])
      topKNN_threshold <- apply(corMat, 2, function(x)cor_threshold_tmp[as.numeric(subLevelModel$y[order(x, decreasing = T)][1:k])])
      topKNN <- ifelse(topKNN_cor >= topKNN_threshold, topKNN, -1)

      topKNN_weight <- apply(topKNN_cor, 2, function(x){
        x <- stats::na.omit(x) # add in 20190613
        if (x[1] == x[length(x)]) {
          rep(1, length(x))
        }else{
          (x - x[length(x)])*(2 - x[length(x)] - x[1])/((x[1] - x[length(x)])*(2 - x[length(x)] - x))
        }
      })

      # if (verbose) {
      #   votes <- sapply(1:ncol(corMat),
      #                   function(i){
      #                     vote <- (stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum))
      #                     vote$x/sum(vote$x)
      #
      #
      #                   })
      #   cat("Summary of the votes \n")
      #   print(summary(unlist(votes)))
      #   graphics::hist(unlist(votes), main = "Votes", breaks = 100)
      # }




      predRes <- sapply(1:ncol(corMat),
                        function(i){
                          vote <- stats::aggregate(topKNN_weight[,i], by = list(topKNN[,i]), sum)
                          maxIdx <- which.max(vote$x)
                          if (max(vote$x/sum(vote$x)) < prob_threshold) {
                            0
                          }else{
                            if (vote$Group.1[maxIdx] == "-1") {
                              0
                            }else{
                              vote$Group.1[maxIdx]
                            }
                          }


                        })
      names(predRes) <- colnames(corMat)

    }


  }

  return(list(predRes = predRes))
}

# Function to generate the mixture model based on the mixtools normalmixEM results

funMixModel <- function(x, mu1, mu2, sd1, sd2, rho1, rho2) {

  dnorm(x, mean = mu1, sd = sd1) * rho1 - dnorm(x, mean = mu2, sd = sd2) * rho2


}



# Function to get the threshold for correlation based on mixture model

getThreshold <- function(mixmdl, verbose = FALSE){

  # if (verbose) {
  #   plot(mixmdl, which = 2)
  # }

  membership <- apply(mixmdl$posterior, 1, which.max)
  m_list <- sort(unique(membership))

  mu_list <- mixmdl$mu
  names(mu_list) <- c(1:length(mu_list))
  mu_list <- mu_list[m_list]

  if (length(mu_list) > 1) {
    idx1 <- as.numeric(names(mu_list)[order(mu_list)][1])
    idx2 <- as.numeric(names(mu_list)[order(mu_list)][2])

    root <- try(uniroot(funMixModel, interval = c(mixmdl$mu[idx1] - mixmdl$sigma[idx1], mixmdl$mu[idx2] + mixmdl$sigma[idx2]),
                        mu1 = mixmdl$mu[idx1], mu2 = mixmdl$mu[idx2],
                        sd1 = mixmdl$sigma[idx1], sd2 = mixmdl$sigma[idx2],
                        rho1 = mixmdl$lambda[idx1], rho2 = mixmdl$lambda[idx2]),
                silent = T)


    if (class(root) != "try-error") {
      # if (verbose) {
      #   abline(v = root$root, col = "red")
      #   abline(v = mixmdl$mu[idx1] + qnorm(0.99) * mixmdl$sigma[idx1], col = "blue")
      # }
      t <- root$root
    }else{
      t <- 0
    }

  }else{
    t <- 0
  }

  return(t)
}


