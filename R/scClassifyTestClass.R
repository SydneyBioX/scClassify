#' An S4 class to stored training model for scClassify
#'
#' @slot testingDataName Name of the testing dataset
#' @slot trainingModelName Name of training dataset
#' @slot ensembleLabel A vector of character indicates the ensemble predicted label
#' @slot ensembleScore A vector of numeric indicates scores of the ensemble predicted label
#' @slot ensembleModelWeights A vector of numeric indicates the weights of each model
#' @slot predictedLabel A data.frame indicates the predicted labels for each base model
#' (each column indicates one base model)
#' @slot predictedMatrix A list indicates the labels of all level predictions for
#'  each base model (each list indicates one base model)
#' @slot resCategories A data.frame indicates the predicted results for each base model
#' (if true test label is provided)
#' @slot featureCombination A vector of character indicates the features and
#' similarity combinations that are trained for this data
#'
#'
#' @importClassesFrom S4Vectors character_OR_NULL




setClass("scClassifyTestRes",
         slots = c(
           testingDataName = "character",
           trainingModelName = "character",
           ensembleLabel = "character_OR_NULL",
           ensembleScore = "numeric",
           ensembleModelWeights = "numeric",
           predictedLabel = "data.frame",
           predictedMatrix = "list",
           resCategories = "data.frame",
           featureCombination = "character"))


