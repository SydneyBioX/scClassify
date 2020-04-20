
setClassUnion("numeric_OR_NULL", members = c("numeric","NULL"))


#' An S4 class to stored training model for scClassify
#'
#' @slot name Name of the training dataset
#' @slot cellTypeTrain A vector of cell type in training dataset
#' @slot cellTypeTree A list indicate a cell type tree
#' @slot features A vector of character indicates the
#' features that are trained for this data
#' @slot model A list stored the training model,
#' including the features that are selected
#' and the cell expression matrix that are used for training
#' @slot modelweights A vector of numeric indicates the weights of each model
#' @slot metaData A DataFrame stored meta data of training model
#'
#' @importClassesFrom S4Vectors character_OR_NULL
#' @importClassesFrom S4Vectors DataFrame
#'
#'
#' @export
setClass("scClassifyTrainModel",
         slots = c(
             name = "character_OR_NULL",
             cellTypeTrain = "character",
             cellTypeTree = "list",
             features = "character",
             model = "list",
             modelweights = "numeric_OR_NULL",
             metaData = "DataFrame"))



#' The scClassifyTrainModel class
#'
#' The scClassifyTrainModel class is designed to stored
#' training model for scClassify
#' @param name Name of the training dataset
#' @param cellTypeTrain A vector of cell type in training dataset
#' @param cellTypeTree A list indicate a cell type tree
#' @param features A vector of character indicates the
#' features that are trained for this data
#' @param model A list stored the training model,
#' including the features that are selected
#' and the cell expression matrix that are used for training
#' @param modelweights A vector of numeric indicates the weights of each model
#' @param metaData A DataFrame stored meta data of training model
#'
#' @return A scClassifyTrainModel object
#'
#' @author Yingxin Lin
#'
#'
#' @importFrom methods new

.scClassifyTrainModel <- function(name,
                                  cellTypeTree,
                                  cellTypeTrain,
                                  features,
                                  model,
                                  modelweights,
                                  metaData) {

    scClassify_train_obj <- new("scClassifyTrainModel",
                                name = name,
                                cellTypeTree = cellTypeTree,
                                cellTypeTrain = cellTypeTrain,
                                features = features,
                                model = model,
                                modelweights = modelweights,
                                metaData = metaData)
    names(scClassify_train_obj@modelweights) <- names(modelweights)
    return(scClassify_train_obj)

}




setMethod("show", "scClassifyTrainModel",
          function(object) {
              cat(paste("Class:", class(object), "\n"))
              cat(paste("Model name:", object@name), "\n")
              cat("Feature selection methods: ")
              cat(object@features, "\n")
              cat("Number of cells in the training data: ")
              cat(length(object@cellTypeTrain), "\n")
              cat("Number of cell types in the training data: ")
              cat(length(object@cellTypeTree[[length(object@cellTypeTree)]]),
                  "\n")
          })



#' An S4 class to stored a list of training models from scClassify
#'
#'
#' @importClassesFrom S4Vectors SimpleList
#' @export

setClass("scClassifyTrainModelList",
         contains = "SimpleList")



#' The scClassifyTrainModelList class
#'
#' @param ... scClassifyTrainModel objects
#'
#' @return A scClassifyTrainModelList object
#'
#' @examples
#'
#' data("trainClassExample_xin")
#' data("trainClassExample_wang")
#' trainClassExampleList <- scClassifyTrainModelList(trainClassExample_xin,
#' trainClassExample_wang
#' )
#'
#' @importFrom S4Vectors SimpleList
#'
#' @export
scClassifyTrainModelList <- function(...) {
    l <- new("scClassifyTrainModelList", S4Vectors::SimpleList(...))
    name_list <- lapply(l@listData, function(x) x@name)
    if (any(duplicated(name_list))) {
        dup_names <- unique(name_list[duplicated(name_list)])
        for (i in seq_along(dup_names)) {
            idx_to_change <- which(name_list %in% dup_names[i])
            name_list[idx_to_change] <- paste(name_list[idx_to_change],
                                              seq_along(idx_to_change),
                                              sep = "_")
        }
    }
    names(l@listData) <- name_list
    l
}




setMethod("show", "scClassifyTrainModelList", function(object) {
    cat(paste("Class:", class(object), "\n"))
    cat("Number of Trainin model: ")
    cat(length(object), "\n")

    cat("Training data name:",
        names(object@listData), "\n")
})


#' Accessors of cellTypeTree for scClassifyTrainModel
#'
#' Methods to access various components of the `scClassifyTrainModel` object.
#'
#' @usage cellTypeTree(x)
#'
#' @param x A `scClassifyTrainModel` object.
#'
#' @examples
#'
#' data(trainClassExample_xin)
#' cellTypeTree(trainClassExample_xin)
#'
#' @return cellTypeTree of the scClassifyTrainModel slot
#'
#' @aliases
#' cellTypeTree,scClassifyTrainModel-method
#' cellTypeTree
#'
#' @export
setGeneric("cellTypeTree", function(x)
    standardGeneric("cellTypeTree"))
setMethod("cellTypeTree", "scClassifyTrainModel",
          function(x) {
              x@cellTypeTree
          })




setGeneric("cellTypeTree<-", function(x, value)
    standardGeneric("cellTypeTree<-"))
setReplaceMethod("cellTypeTree", "scClassifyTrainModel",
                 function(x, value) {
                     x@cellTypeTree <- value
                     x
                 })



#' Accessors of features for scClassifyTrainModel
#'
#' Methods to access various components of the `scClassifyTrainModel` object.
#'
#' @usage features(x)
#'
#' @param x A `scClassifyTrainModel` object.
#'
#' @examples
#'
#' data(trainClassExample_xin)
#' features(trainClassExample_xin)
#'
#' @return features of the scClassifyTrainModel slot
#'
#' @aliases
#' features,scClassifyTrainModel-method
#' features
#'
#' @export
setGeneric("features", function(x) standardGeneric("features"))
setMethod("features", "scClassifyTrainModel", function(x) {
    x@features
})



setGeneric("features<-", function(x, value)
    standardGeneric("features<-"))
setReplaceMethod("features", "scClassifyTrainModel",
                 function(x, value) {
                     x@features <- value
                     x
                 })








#' Accessors of name for scClassifyTrainModel
#'
#' Methods to access various components of the `scClassifyTrainModel` object.
#'
#' @usage name(x)
#'
#' @param x A `scClassifyTrainModel` object.
#'
#' @examples
#'
#' data(trainClassExample_xin)
#' name(trainClassExample_xin)
#'
#' @return name of the scClassifyTrainModel slot
#'
#' @aliases
#' name,scClassifyTrainModel-method
#' name
#'
#'
#' @export
setGeneric("name", function(x) standardGeneric("name"))
setMethod("name", "scClassifyTrainModel", function(x) {
    x@name
})



setGeneric("name<-", function(x, value)
    standardGeneric("name<-"))
setReplaceMethod("name", "scClassifyTrainModel",
                 function(x, value) {
                     x@name <- value
                     x
                 })








#' Accessors of cellTypeTrain for scClassifyTrainModel
#'
#' Methods to access various components of the `scClassifyTrainModel` object.
#'
#' @usage cellTypeTrain(x)
#'
#' @param x A `scClassifyTrainModel` object.
#'
#' @examples
#'
#' data(trainClassExample_xin)
#' cellTypeTrain(trainClassExample_xin)
#'
#' @return cellTypeTrain of the scClassifyTrainModel slot
#'
#' @aliases
#' cellTypeTrain,scClassifyTrainModel-method
#' cellTypeTrain
#'
#'
#' @export
setGeneric("cellTypeTrain", function(x) standardGeneric("cellTypeTrain"))
setMethod("cellTypeTrain", "scClassifyTrainModel", function(x) {
    x@cellTypeTrain
})



setGeneric("cellTypeTrain<-", function(x, value)
    standardGeneric("cellTypeTrain<-"))
setReplaceMethod("cellTypeTrain", "scClassifyTrainModel",
                 function(x, value) {
                     x@cellTypeTrain <- value
                     x
                 })








#' Accessors of modelweights for scClassifyTrainModel
#'
#' Methods to access various components of the `scClassifyTrainModel` object.
#'
#' @usage modelweights(x)
#'
#' @param x A `scClassifyTrainModel` object.
#'
#' @examples
#'
#' data(trainClassExample_xin)
#' modelweights(trainClassExample_xin)
#'
#' @return modelweights of the scClassifyTrainModel slot
#'
#'
#' @aliases
#' modelweights,scClassifyTrainModel-method
#' modelweights
#'
#'
#' @export
setGeneric("modelweights", function(x) standardGeneric("modelweights"))
setMethod("modelweights", "scClassifyTrainModel", function(x) {
    x@modelweights
})



setGeneric("modelweights<-", function(x, value)
    standardGeneric("modelweights<-"))
setReplaceMethod("modelweights", "scClassifyTrainModel",
                 function(x, value) {
                     x@modelweights <- value
                     x
                 })






#' Accessors of model for scClassifyTrainModel
#'
#' Methods to access various components of the `scClassifyTrainModel` object.
#'
#' @usage model(x)
#'
#' @param x A `scClassifyTrainModel` object.
#'
#' @return model of the scClassifyTrainModel slot
#'
#' @examples
#'
#' data(trainClassExample_xin)
#' model(trainClassExample_xin)
#'
#' @aliases
#' model,scClassifyTrainModel-method
#' model
#'
#'
#' @export
setGeneric("model", function(x) standardGeneric("model"))
setMethod("model", "scClassifyTrainModel", function(x) {
    x@model
})



setGeneric("model<-", function(x, value)
    standardGeneric("model<-"))
setReplaceMethod("model", "scClassifyTrainModel",
                 function(x, value) {
                     x@model <- value
                     x
                 })

