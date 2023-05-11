#' @title Dichotomize via Recursive Partitioning
#' 
#' @description
#' Dichotomize a single \link[base]{numeric} predictor 
#' using recursive partitioning \link[rpart]{rpart}.
#' 
#' @param y a \link[survival]{Surv} object, a \link[base]{logical} \link[base]{vector}, 
#' or a \link[base]{numeric} \link[base]{vector}, the endpoint
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, the predictor
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' 
#' The \link[rpart]{labels.rpart} of the first node of 
#' the recursive partitioning \link[rpart]{rpart} tree
#' is considered as the dichotomizing branch of the \link[base]{numeric} predictor \code{x}.
#' 
#' The term *branch* indicates the combination of an inequality sign
#' (\link[base]{>}, \link[base]{>=}, \link[base]{<} and \link[base]{<=}) 
#' and a \link[base]{numeric} threshold \eqn{a}.
#' 
#' This branch is further processed, such that
#' \itemize{
#' \item {\eqn{<a} is returned as \eqn{\geq a}}
#' \item {\eqn{\leq a} is returned as \eqn{>a}}
#' \item {\eqn{> a} and \eqn{\geq a} are returned as is.}
#' }
#' This step is necessary for a narrative of *higher than the threshold*.
#' 
#' @note
#' 
#' In the recursive partitioning \link[rpart]{rpart} tree, we use the following tuning parameters 
#' in \link[rpart]{rpart.control}, 
#' \itemize{
#' \item {\code{maxdepth = 2L}, because only the first node is needed}
#' \item {\code{cp = .Machine$double.eps}, so that a split is guaranteed 
#' no matter how small improvement in overall \eqn{R^2} is}
#' }
#' 
#' 
#' @returns 
#' \link{rpartD} returns a \link[base]{language} object.
#' 
#' @examples
#' # see ?rpart::rpart examples
#' data(kyphosis, package = 'rpart')
#' rpartD(y = kyphosis$Kyphosis, x = kyphosis$Start)
#' rpartD(y = kyphosis$Kyphosis, x = -kyphosis$Start)
#' 
#' @importFrom rpart rpart rpart.control
#' @export
rpartD <- function(y, x, ...) {
  tree <- rpart(
    formula = y ~ x, 
    data = data.frame(y = y, x = x), 
    # method = ???,
    control = rpart.control(
      cp = .Machine$double.eps, # to force a split even if the overall lack of fit is not decreased
      maxdepth = 2L # only the first node is needed
    ))
  if (!length(tree$splits)) stop('we must force a split')
  
  labs <- labels(tree) # ?rpart:::labels.rpart
  ret <- str2lang(labs[2L]) # first node!!!
  if (ret[[1L]] == '<=') {
    ret[[1L]] <- quote(`>`)
  } else if (ret[[1L]] == '<') {
    ret[[1L]] <- quote(`>=`)
  } # else do nothing
  
  ret[[2L]] <- alist(x = )$x # 2nd element is 'empty' :)
  # keep 2nd element (i.e., do NOT use `call('>', 3)`)
  # it is easer to construct ?base::eval later :)
  
  ret[[3L]] <- tree$splits[1L, 4L] # threshold, in case `labels` are truncated due to `digits`
  return(ret)

}



