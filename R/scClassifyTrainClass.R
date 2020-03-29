
setClassUnion("numeric_OR_NULL", members = c("numeric","NULL"))


#' An S4 class to stored training model for scClassify
#'
#' @slot name Name of the training dataset
#' @slot cellTypeTrain A vector of cell type in training dataset
#' @slot cellTypeTree A list indicate a cell type tree
#' @slot features A vector of character indicates the features that are trained for this data
#' @slot model A list stored the training model, including the features that are selected
#' and the cell expression matrix that are used for training
#' @slot modelweights A vector of numeric indicates the weights of each model
#' @slot metaData A DataFrame stored meta data of training model
#'
#' @importClassesFrom S4Vectors character_OR_NULL
#' @importClassesFrom S4Vectors DataFrame
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
#' The scClassifyTrainModel class is designed to stored training model for scClassify
#' @param name Name of the training dataset
#' @param cellTypeTrain A vector of cell type in training dataset
#' @param cellTypeTree A list indicate a cell type tree
#' @param features A vector of character indicates the features that are trained for this data
#' @param model A list stored the training model, including the features that are selected
#' and the cell expression matrix that are used for training
#' @param modelweights A vector of numeric indicates the weights of each model
#' @param metaData A DataFrame stored meta data of training model
#'
#' @return A scClassifyTrainMode object
#'
#' @author Yingxin Lin
#'
#' @docType class
#' @importFrom methods new
#' @export

scClassifyTrainModel <- function(name,
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
            cat(length(object@cellTypeTree[[length(object@cellTypeTree)]]), "\n")
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
#' @export

scClassifyTrainModelList <- function(...) {
  l <- new("scClassifyTrainModelList", S4Vectors::SimpleList(...))
  names(l@listData) <- lapply(l@listData, function(x) x@name)
  l
}




setMethod("show", "scClassifyTrainModelList", function(object) {
  cat(paste("Class:", class(object), "\n"))
  cat("Number of Trainin model: ")
  cat(length(object), "\n")

  cat("Training data name:",
      names(object@listData), "\n")
})


