
#' Fit learning curve for accuracy matrix
#' @param accMat Matrix of accuracy rate where column indicate
#' different sample size
#' @param n Vector indicates the sample size
#' @param auto_initial whether automatical intialise
#' @param a input the parameter a starting point
#' @param b input the parameter a starting point
#' @param c input the parameter a starting point
#' @param d_list range of d
#' @param fitmodel "nls", "nls_mix", "gam"
#' @param plot indicates whether plot or not
#' @param verbose indicates whether verbose or not
#' @return list of results
#' @examples
#' set.seed(2019)
#' n <- seq(20, 10000, 100)
#' accMat <- do.call(cbind, lapply(1:length(n), function(i){
#' tmp_n <- rep(n[i], 50)
#' y <- -2/(tmp_n^0.8) + 0.95 + rnorm(length(tmp_n), 0, 0.02)
#' }))
#' res <- learningCurve(accMat = accMat, n)
#' N <- getN(res, acc = 0.9)
#'
#' @author Yingxin Lin
#' @import ggplot2
#' @importFrom minpack.lm nlsLM
#' @importFrom stats quantile coef predict
#' @importFrom mgcv gam
#' @export

learningCurve <- function(accMat, n, auto_initial = TRUE,
                          a = NULL, b = NULL, c = NULL, d_list = NULL,
                          fitmodel = c("nls", "nls_mix", "gam"),
                          plot = TRUE, verbose = TRUE){

    ## TODO: to make several ok TRUE
    fitmodel <- match.arg(fitmodel, c("nls", "nls_mix", "gam"),
                          several.ok = FALSE)



    if (any(c("matrix", "data.frame") %in% is(accMat))) {
        if (ncol(accMat) != length(n)) {
            stop("Number of column doesn't match with the length of n")
        }
        dat <- data.frame(y = colMeans(accMat),
                          y_25 = apply(accMat, 2,
                                       function(x) quantile(x, 0.25)),
                          y_75 = apply(accMat, 2,
                                       function(x) quantile(x, 0.75))
        )
    }

    if ( "list" %in% is(accMat)) {
        if (length(accMat) != length(n)) {
            stop("Number of column doesn't match with the length of n")
        }
        dat <- data.frame(y = unlist(lapply(accMat, mean)),
                          y_25 = unlist(lapply(accMat,
                                               function(x) quantile(x, 0.25))),
                          y_75 = unlist(lapply(accMat,
                                               function(x) quantile(x, 0.75))))
    }



    # Fitting model
    model <- list()
    for (i in seq(ncol(dat))) {
        if (fitmodel == "nls_mix") {
            model[[i]] <- fitLC_mixture(dat[,i], n,
                                        b = b, d_list = d_list,
                                        verbose = verbose)

        } else if (fitmodel == "gam") {
            model[[i]] <- mgcv::gam(acc ~ s(n, bs = "cr"),
                                    data = data.frame(acc = dat[,i], n = n),
                                    method = "REML")
        }

        else{
            model[[i]] <- fitLC(dat[,i], n, auto_initial = auto_initial,
                                a = a, b = b, c = c,
                                verbose = verbose, b_max = 1)
            if (i == 1) {
                b_max <- coef(model[[i]])["b"]
            }
        }

    }


    # Get fitted value
    dat$n <- n
    new_n <- data.frame(n = seq(min(n), max(n), by = 0.1))
    fit <- lapply(model, function(x) predict(x, newdata = new_n))

    names(model) <- names(fit) <- c("mean", "quantile_25", "quantile_75")
    fit[["n"]] <- new_n$n
    dat_fit <- data.frame(do.call(cbind, fit))

    if (plot) {


        cols <- c("Mean" = "#c8133b","Quantile25/75" = "#ea8783")
        g <-  ggplot2::ggplot(dat, ggplot2::aes(x = n, y = y))  +
            xlab("N") + ylab("Accuracy Rate") +
            ggplot2::geom_point(alpha = 0.7) +
            ggplot2::geom_line(data = dat_fit, aes(x = n,
                                                   y = mean,
                                                   color = "Mean"),
                               linetype = "solid", size = 1) +
            ggplot2::geom_line(data = dat_fit, aes(x = n, y = quantile_25,
                                                   color = "Quantile25/75"),
                               linetype = "dashed") +
            ggplot2::geom_line(data = dat_fit, aes(x = n, y = quantile_75,
                                                   color = "Quantile25/75"),
                               linetype = "dashed") +
            ggplot2::scale_color_manual(values = cols) +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = "bottom")

        return(list(model = model, fit = fit, plot = g))

    }else{
        return(list(model = model, fit = fit))
    }
}


init_fit <- function(acc, n){
    fit <- lm(log(n) ~ log(max(acc) + 0.001 - acc))

    c <- -coef(fit)[2]
    a <- -exp(coef(fit)[1])
    return(list(a = a, c = c))
}

# Function to fit the learning curve by inverse power function

fitLC <- function(acc, n, auto_initial = TRUE,
                  a = NULL, b = NULL, c = NULL, verbose = TRUE,
                  b_max = NULL){
    dat_train <- data.frame(n = n, acc = acc)

    if (is.null(b)) {
        b = 0.9
    }

    if (is.null(b_max)) {
        b_max = 1
    }

    if (auto_initial) {
        init <- init_fit(acc, n)

        lower_bound <- c(-max(max(dat_train$acc) + 0.01, 1)*(10^init$c),
                         -Inf, 0)

        learning_curve <-  try({minpack.lm::nlsLM(acc ~ I(1 / n^(c) * a) + b,
                                                  data = dat_train,
                                                  start = list(a = init$a,
                                                               c = init$c,
                                                               b = b),
                                                  upper = c(10, Inf, b_max),
                                                  lower = lower_bound
        )
        }, silent = TRUE)


    }else if (!is.null(a) & !is.null(b) & !is.null(c)) {
        # For the case that user supplies starting point
        learning_curve <- stats::nls(acc ~ I(1 / n^(c) * a) + b,
                                     data = dat_train,
                                     start = list(a = a, c = c, b = b))
    }else{
        # this starting point may not work all the time
        learning_curve <- stats::nls(acc ~ I(1 / n^(c) * a) + b,
                                     data = dat_train,
                                     start = list(a = -10, c = 1, b = 1))
    }

    para <- coef(learning_curve)
    names(para) <- c("a1", "c1", "b1")

    if (verbose) {
        print(summary(learning_curve))
    }

    learning_curve$para <- para
    learning_curve$fit_model <- "nls"

    return(learning_curve)
}




fitLC_mixture <- function(acc, n,  b = NULL,
                          d_list = NULL,
                          verbose = TRUE){
    dat_train <- data.frame(n = n, acc = acc)

    if (is.null(b)) {
        b = 0.9
    }

    if (is.null(d_list)) {
        d_list <- seq(min(n) + 10, round(max(n)/2), 10)
    }

    init <- init_fit(acc, n)

    rss <- c()
    for (i in seq_len(length(d_list))) {
        d <- d_list[i]
        learning_curve <-  try({minpack.lm::nlsLM(acc ~ H1(d - n) *
                                                      I(1 / n^(c) * a) +
                                                      H(n - d)  *
                                                      I(1 / n^(c1) *
                                                            a*d^(c1 - c))  + b,
                                                  data = dat_train,
                                                  start = list(a = init$a,
                                                               c = init$c,
                                                               b = b,
                                                               c1 = init$c),
                                                  upper = c(Inf, Inf, 1, Inf),
                                                  lower = c(-Inf, -Inf, 0,
                                                            -Inf),
                                                  control = list(maxiter =
                                                                     1000))},
                               silent = TRUE)

        if (!"try-error" %in% is(learning_curve)) {
            rss <- c(rss, summary(learning_curve)$sigma)
        }else{
            rss <- c(rss, Inf)
        }

    }


    d_opt <- d_list[which.min(rss)]

    learning_curve <-  minpack.lm::nlsLM(acc ~ H1(d_opt - n) *
                                             I(1 / n^(c) * a) +
                                             H(n - d_opt)  *
                                             I(1 / n^(c1) * a*d_opt^(c1-c))  +
                                             b,
                                         data = dat_train,
                                         start = list(a = init$a, c = init$c,
                                                      b = 0.95,
                                                      c1 = init$c),
                                         upper = c(Inf, Inf, 1, Inf),
                                         lower = c(-Inf, -Inf, 0, -Inf),
                                         control = list(maxiter = 1000)

    )

    para <- coef(learning_curve)
    para <- c(para, para[1]*d_opt^(para[4] - para[2]), d_opt)
    names(para) <- c("a1", "c1", "b1", "c2", "a2", "d")

    if (verbose) {
        print(summary(learning_curve))
    }
    learning_curve$para <- para
    learning_curve$fit_model <- "nls_mix"
    return(learning_curve)
}


#' Function to get the required N given by the accuracy and
#' the learning curve model
#' @param res model results returned by \code{learning_curve} function
#' @param acc accuracy that are quired
#' @return sample size that are required
#'
#' @examples
#' set.seed(2019)
#' n <- seq(20, 10000, 100)
#' accMat <- do.call(cbind, lapply(1:length(n), function(i){
#' tmp_n <- rep(n[i], 50)
#' y <- -2/(tmp_n^0.8) + 0.95 + rnorm(length(tmp_n), 0, 0.02)
#' }))
#' res <- learningCurve(accMat = accMat, n)
#' N <- getN(res, acc = 0.9)
#'
#' @export

getN <- function(res, acc = 0.9) {
    para <- coef(res$model$mean)
    names(para) <- c("a", "c", "b")
    if (acc > para["b"]) {
        stop("Required accuracy is too high to achieve :(")
    }
    N <- round(exp(1/para["c"]*log(para["a"]/(acc - para["b"]))))
    names(N) <- NULL

    return(N)
}

H <- function(x) as.numeric(x > 0)
H1 <- function(x) as.numeric(x >= 0)

