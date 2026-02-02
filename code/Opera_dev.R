blockToSeries <- function(X){
  if (is.null(dim(X)) || length(dim(X)) > 3){
    stop("X must be of dimension 2 or 3")
  } 
  # if X is a matrix, we convert it to a vector 
  if (length(dim(X)) == 2) {
    return(c(t(X)))
  } 
  if (length(dim(X)) == 3) {
    return(array(c(aperm(X,c(2,1,3))),dim = c(dim(X)[1]*dim(X)[2],dim(X)[3])))
  }
}




#' Compute oracle predictions
#'
#' The function \code{oracle} performs a strategie that cannot be defined online
#' (in contrast to \link{mixture}). It requires in advance the knowledge of the whole
#' data set \code{Y} and the expert advice to be well defined.
#' Examples of oracles are the best fixed expert, the best fixed convex
#' combination rule, the best linear combination rule, or the best expert
#' that can shift a few times.
#'
#' @param Y  A vector containing the observations
#' to be predicted.
#' @param experts A matrix containing the experts
#' forecasts. Each column corresponds to the predictions proposed by an expert
#' to predict \code{Y}.  It has as many columns as there are experts.
#' @param model A character string specifying the oracle to use or a list with a component \code{name} specifying the oracle and any additional parameter needed.
#' Currently available oracles are:
#' \describe{
#'    \item{'expert'}{The best fixed (constant over time) expert oracle.}
#'    \item{'convex'}{The best fixed convex combination (vector of non-negative weights that sum to 1)}
#'    \item{'linear'}{The best fixed linear combination of expert}
#'    \item{'shifting'}{It computes for all number $m$ of stwitches the
#' sequence of experts with at most $m$ shifts that would have performed the
#' best to predict the sequence of observations in \code{Y}.}
#' }
#' @param loss.type A string or a list with a component 'name' specifying
#' the loss function considered to evaluate the performance. It can be
#' 'square', 'absolute', 'percentage', or 'pinball'. In the case of the pinball loss, the quantile 
#' can be provided by assigning to loss.type a list of two elements: 
#' \describe{
#'      \item{name}{A string defining the name of the loss function (i.e., 'pinball')}
#'      \item{tau}{ A number in \code{[0,1]} defining the quantile to be predicted. The default value is 0.5 to predict the median.}
#' } 
#' 
#' @param lambda A positive number used by the 'linear' oracle only. 
#' A possible $L_2$ regularization parameter for computing the linear oracle 
#' (if the design matrix is not identifiable)
#' @param niter A positive integer for 'convex' and 'linear' oracles 
#' if direct computation of the oracle is not implemented. 
#' It defines the number of optimization steps to perform in 
#' order to approximate the oracle (default value is 3).
#' 
#' @param awake A matrix specifying the
#' activation coefficients of the experts. Its entries lie in \code{[0,1]}.
#' Possible if some experts are specialists and do not always form and suggest
#' prediction. If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}. Remark that to compute the best expert oracle, 
#' the performance of unactive (or partially active) experts is computed by using 
#' the prediction of the uniform average of active experts.
#' 
#' @param ... Additional parameters
#' that are passed to \code{\link{optim}} function is order to perform convex optimization 
#' (see parameter \code{niter}).
#'
#' @return An object of class 'oracle' that contains:
#' \item{loss}{ The average loss suffered by the oracle. For the 'shifting' oracle,
#' it is a vector of length \code{T} where
#' \code{T} is the number of instance to be predicted (i.e., the length of the
#' sequence \code{Y}). The value of $loss(m)$ is the loss
#' (determined by the parameter \code{loss.type}) suffered by the
#' best sequence of expert with at
#' most $m-1$ shifts.
#' }
#' \item{coefficients}{ Not for the 'shifting' oracle. A vector containing the best weight vector corresponding to the oracle. }
#' \item{prediction}{ Not for the 'shifting' oracle. A vector containing the
#' predictions of the oracle.  }
#' \item{rmse}{If loss.type is the square loss (default) only.
#' The root mean square error (i.e., it is the square root of \code{loss}.}
#' 
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @export oracle

oracle <- function(Y, experts, model = "convex", loss.type = "square", awake = NULL, 
                   lambda = NULL, niter = NULL, ...) UseMethod("oracle")


#' @export 
oracle.default <- function(Y, experts, model = "convex", loss.type = "square", awake = NULL, 
                           lambda = NULL, niter = NULL, ...) {
  
  # Test that Y and experts have correct dimensions
  if (is.null(Y) || is.null(experts)) {
    stop("Y and experts should not be null")
  }
  if (length(Y) == 1) {
    experts <- as.matrix(experts)
    if (nrow(experts) == 1 || ncol(experts) == 1) {
      experts <- matrix(experts, nrow = 1)
    } else {
      stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
    }
  }
  d = 1
  # We convert the data to 1-dimensional data if needed
  if (length(dim(Y)) == 2 && length(dim(experts)) == 3) {
    d = dim(Y)[2]
    T = dim(Y)[1]
    
    Y = blockToSeries(Y)
    experts = blockToSeries(experts)
    if (!is.null(awake)) {
      awake = blockToSeries(awake)
    }
  }
  
  if (!(length(Y) == nrow(experts))) {
    stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
  }
  
  if (is.null(loss.type)) {
    loss.type <- list(name = "square")
  }
  if (!is.list(loss.type)) {
    loss.type <- list(name = loss.type)
  }
  if (!(loss.type$name %in% c("pinball", "square", "percentage", "absolute"))) {
    stop("loss.type should be one of these: 'absolute', 'percentage', 'square', 'pinball'")
  }
  if (!is.null(loss.type$tau) && loss.type$name != "pinball") {
    warning("Unused parameter tau (loss.type != 'pinball')")
  }
  if (!is.null(lambda) && model != "linear") {
    warning("Unused lambda parameter (model != 'linear')")
  }
  if (is.null(lambda) && model == "linear") 
    lambda <- 0 
  
  if (!is.null(niter) && model != "convex" && model != "linear") {
    warning("Unused niter parameter (model should be 'convex' or 'linear')")
  }
  if (is.null(niter)) 
    niter <- 3
  
  if ((!is.null(awake) || sum(is.na(experts) > 0)) && model != "convex" && model != 
        "shifting") {
    if (model != "expert") {
      stop(paste("Sleeping or missing values not allowed for best", model, "oracle."))
    }
    else {
      warning("When experts are unactive (or sleeping), their prediction are replaced with the uniform average of active experts")
    }
  }
  
  
  if (!(model %in% c("convex", "linear", "shifting", "expert"))) {
    stop("Wrong model specification")
  }
  if (min(Y) <= 0 && loss.type$name == "percentage") {
    stop("Y should be non-negative for percentage loss function")
  }
  names.experts <- colnames(experts)
  experts <- matrix(as.numeric(as.matrix(experts)), nrow = length(Y))
  colnames(experts) <- names.experts
  
  # if we are looking for the best convex combination of experts
  if (model == "convex") {
    res <- bestConvex(Y, experts, awake = awake, loss.type = loss.type, niter = niter, 
                      ...)
  }

  if (model == "linear") {
    res <- tryCatch(
      bestLinear(Y, experts, lambda = lambda, loss.type = loss.type),
      error = function(err) {
        bestLinear(Y, experts, lambda = 1e-14, loss.type = loss.type) 
      })
  }
  
  if (model == "shifting") {
    res <- bestShifts(Y, experts, awake = awake, loss.type = loss.type)
  }
  
  if (!is.null(awake)) {
    pond <- apply(awake,1,mean)
    pred.unif <- apply(experts * awake, 1,mean) /pond
    experts.pred <- experts * awake + pred.unif * (1-awake)
  } else {
    experts.pred <- experts
  }
  
  loss.experts <- apply(apply(experts.pred, 2, function(x) {
    loss(x, Y, loss.type = loss.type)
  }), 2, mean,na.rm=T)
  
  if (model == "expert") {
    best.loss <- min(loss.experts)
    coefficients <- (loss.experts == best.loss)/sum(loss.experts == best.loss)
    best.expert <- which(coefficients > 0)[1]
    res <- list(loss = best.loss, coefficients = coefficients, prediction = experts.pred[, 
                                                                                         best.expert], loss.experts = loss.experts)
  }
  
  res$d <- d
  res$loss.experts <- loss.experts
  res$model <- model
  res$loss.type <- loss.type
  res$call <- match.call()
  if (model != "shifting") {
    res$residuals <- Y - res$prediction
    res$loss <- mean(loss(res$prediction, Y, loss.type))
    res$rmse <- sqrt(mean(loss(res$prediction, Y, "square")))
  }
  
  # we convert the data back to d-dimensional series if needed
  if (d>1){
    Y <- seriesToBlock(Y,d)
    experts <- seriesToBlock(experts,d)
    if (!is.null(res$residuals)) {
      res$residuals <- seriesToBlock(res$residuals,d)
    }
    res$prediction <- seriesToBlock(res$prediction,d)
    if (!is.null(awake)) {
      awake <- seriesToBlock(awake,d)
    }
  }
  res$Y <- Y
  res$experts <- experts
  res$awake <- awake
  
  
  class(res) <- "oracle"
  return(res)
} 




# best convex oracle
#' @importFrom stats runif optim 
bestConvex <- function(y, experts, awake = NULL, loss.type = list(name = "square"), 
                       niter = 1, ...) {
  experts <- matrix(as.numeric(as.matrix(experts)), nrow = length(y))
  N <- ncol(experts)

  # if there are no NA and if awake is null we can perform an exact resolution for
  # the square loss
  idx.na <- which(is.na(experts))
  
  if (length(idx.na) == 0 && is.null(awake) && loss.type$name == "square") {
    y.na <- is.na(y)
    y <- y[!y.na]
    x <- experts[!y.na, ]
    eq <- paste("y ~ x-1")
    
    
    Q <- crossprod(x)
    c <- crossprod(x, y)
    A <- cbind(1, diag(nrow(Q)))
    b <- c(1, rep(0, nrow(Q)))
    m <- 1
    res <- tryCatch({
      if (!requireNamespace("quadprog", quietly = TRUE)) {
        warning("The quadprog package must be installed to use this functionality")
        #Either exit or do something without quadprog
        return(NULL)
      } else {
        quadprog::solve.QP(Dmat = Q, dvec = c, Amat = A, bvec = b, meq = m)
      }
    }, error = function(e) {
      NULL
    })
    if (!is.null(res)) {
      coefficients <- matrix(res$solution, ncol = N)
      prediction <- experts %*% t(coefficients)
      bestLoss <- mean(loss(x = prediction, y))
    }
  } else {
    res <- NULL
  }
  if (is.null(res)) {
    warning("The best convex oracle is only approximated (using optim).")
    if (is.null(awake)) {
      awake <- as.matrix(array(1, dim(experts)))
    }
    awake[idx.na] <- 0
    experts[idx.na] <- 0
    if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
      y[which(is.na(y))]<- 0
      awake[idx.na] <- 1 #to avoid 0 division
      experts[idx.na] <- 0
    }
    
    lossp <- function(p) {
      return(lossConv(p, y, experts, awake, loss.type))
    }
    
    best_p <- rep(0, N)
    bestLoss <- exp(700)
    
    for (i in 1:niter) {
      # Random initialization
      p <- runif(N, 0, 1)
      p <- p/sum(p)

      # Convex optimization
      w <- optim(p, lossp, gr = NULL, lower = 1e-20, method = "L-BFGS-B", ...)

      # Projection on the simplex
      w <- pmax(w$par, 0)
      l <- lossp(w)
      if (bestLoss > l) {
        bestLoss <- l
        best_p <- w
      }
    }
    coefficients <- matrix(best_p, ncol = N)
    coefficients <- coefficients/apply(coefficients, 1, sum)
    pond <- awake %*% t(coefficients)
    prediction <- ((as.numeric(experts) * awake) %*% t(coefficients))/pond
  }
  res <- list(coefficients = coefficients, prediction = prediction)
  return(res)
} 


# best linear oracle
#' @importFrom stats rnorm optim 
bestLinear <- function(y, experts, lambda = 0, loss.type = list(name = "square"), 
                       niter = 1, ...) {
  experts <- as.matrix(experts)
  N <- ncol(experts)
  
  coefficients <- NULL
  if (loss.type$name == "square") {
    coefficients <- 
      tryCatch(
        solve(lambda * diag(1, ncol(experts)) + t(experts) %*% experts, t(experts) %*% y),
        error = function(err){
          lambda = 1e-14
          solve(lambda * diag(1, ncol(experts)) + t(experts) %*% experts, t(experts) %*% y)
          warning("Ill conditioned problem. Regularized with lambda = 1e-14.")
        })
    
  } else if (loss.type$name == "pinball") {
    if (is.null(loss.type$tau)) {
      loss.type$tau <- 0.5
    }
    if (!requireNamespace("quantreg", quietly = TRUE)) {
      warning("The quantreg package must be installed to use this functionality")
      #Either exit or do something without quantreg
      return(NULL)
    } else {
      coefficients <- tryCatch({
        quantreg::rq(y ~ experts - 1, tau = loss.type$tau)$coefficients
      }, error = function(e) {
        NULL
      })
    }
  }
  if (is.null(coefficients)) {
    warning("The best linear oracle is only approximated (using optim).")
    lossu <- function(u) {
      return(mean(loss(x = experts %*% matrix(u, nrow = ncol(experts)), y = y, 
                       loss.type = loss.type)))
    }
    
    best_u <- rep(0, N)
    bestLoss <- exp(700)
    
    for (i in 1:niter) {
      # Random initialization
      u <- rnorm(N, 0, 1)
      
      # Convex initialization
      w <- optim(u, lossu, gr = NULL, ...)$par
      l <- lossu(w)
      if (bestLoss > l) {
        bestLoss <- l
        best_u <- w
      }
    }
    coefficients <- matrix(best_u, nrow = N)
  }
  
  prev <- experts %*% coefficients
  return(list(coefficients = c(coefficients), prediction = prev))
} 


#best sequence of experts oracle
bestShifts <- function(y, experts, awake = NULL, loss.type = list(name = "square")) {
  N <- ncol(experts)
  T <- nrow(experts)
  INF <- exp(700)
  # m-1 shifts, expert
  
  # Full activation if unspecified
  if (is.null(awake)) {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    y[which(is.na(y))]<- 0
    awake[idx.na] <- 1 #to avoid 0 division
    experts[idx.na] <- 0
  }
  
  L <- array(INF, dim = c(T, N))
  L[1, ] <- 0
  
  for (t in 1:T) {
    Et1 <- which(awake[t - 1, ] > 0)
    Et <- which(awake[t, ] > 0)
    # for (l in 1:3) {
      instanceLoss <- loss(x = experts[t, ], y = y[t], loss.type = loss.type) * awake[t, ]
    # }
    L[1:t, -Et] <- INF
    if (t > 1) {
      L1 <- L[1:t, Et1]
      idx_min <- apply(L1[, ], 1, order)[1:2, ]
      for (m in t:2) {
        for (i in Et) {
          if (idx_min[1, m - 1] == i) 
            aux <- idx_min[2, m - 1] else aux <- idx_min[1, m - 1]
            
            if (L[m, i] < L1[m - 1, aux]) 
              L[m, i] <- L[m, i] + instanceLoss[i] else L[m, i] <- L1[m - 1, aux] + instanceLoss[i]
        }
      }
    }
    L[1, ] <- L[1, ] + instanceLoss
  }
  loss.experts <- L[, ]/T
  loss <- apply(loss.experts, 1, min)
  for (i in 2:T) {
    if (loss[i] > loss[i-1]) {
      loss[i] <- loss[i-1]
    }
  }
  res <- list(loss = loss)
  if (is.list(loss.type) && loss.type$name == "square") {
    res <- list(loss = loss, rmse = sqrt(loss))
  }
  return(res)
} 




loss <- function(x, y, loss.type = "square") {
  
  if (!is.list(loss.type)) {
    loss.type <- list(name = loss.type)
  }
  if (is.null(loss.type$tau) && loss.type$name == "pinball") {
    loss.type$tau <- 0.5
  }
  
  if (loss.type$name == "square") 
    l <- (x - y)^2 else if (loss.type$name == "absolute") 
    l <- abs(x - y) else if (loss.type$name == "percentage") 
    l <- abs(x - y)/y else if (loss.type$name == "pinball") 
    l <- (loss.type$tau - (y < x)) * (y - x) else stop("loss.type should be one of these: 'absolute', 'percentage', 'square', 'pinball'")
  return(l)
} 



lossConv <- function(p, y, experts, awake = NULL, loss.type = "square") {
  
  experts <- matrix(as.numeric(as.matrix(experts)), nrow = length(y))
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  p <- matrix(as.numeric(as.matrix(p)), nrow = N)
  
  # Experts are always active if awake is unspecified
  if (is.null(awake)) {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  awake <- matrix(as.numeric(as.matrix(awake)), nrow = length(y))
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    y[which(is.na(y))]<- 0
    awake[idx.na] <- 1 #to avoid 0 division
    experts[idx.na] <- 0
  }
  
  pond <- awake %*% p
  pred <- ((experts * awake) %*% p)/pond
  l <- mean(loss(pred, y, loss.type = loss.type))
  return(l)
} 

lossPred <- function(x, y, pred = NULL, loss.type = "square", loss.gradient = FALSE) {
  #third argument "pred" is needed to estimate the sub-gradient in the case of loss.gradient=true
  if (!is.list(loss.type)) {
    loss.type <- list(name = loss.type)
  }
  if (is.null(loss.type$tau) && loss.type$name == "pinball") {
    loss.type$tau <- 0.5
  }
  npred <- length(pred)
  nx <- length(x)
  if (npred > 1 && nx > 1) {
    if (!loss.gradient) {
      if (loss.type$name == "square") 
        l <- matrix(rep((x - y)^2, npred), ncol = npred) else if (loss.type$name == "absolute") 
          l <- matrix(rep(abs(x - y), npred), ncol = npred) else if (loss.type$name == "percentage") 
            l <- matrix(rep(abs(x - y)/y, npred), ncol = npred) else if (loss.type$name == "pinball") 
              l <- matrix(rep(((y < x) - loss.type$tau) * (x - y), npred), ncol = npred)
    } else {#gradient trick : we linearize the loss x->l(x) with x->l'(x)*x
      if (loss.type$name == "square") 
        l <- 2 * t(matrix(rep(pred - y, nx), ncol = nx)) * matrix(rep(x, 
                                                                      npred), ncol = npred) else if (loss.type$name == "absolute") 
                                                                        l <- t(matrix(rep(sign(pred - y), nx), ncol = nx)) * matrix(rep(x, 
                                                                                                                                        npred), ncol = npred) else if (loss.type$name == "percentage") 
                                                                                                                                          l <- matrix(rep(x, npred), ncol = npred)/y * t(matrix(rep(sign(pred - 
                                                                                                                                                                                                           y), nx), ncol = nx)) else if (loss.type$name == "pinball") 
                                                                                                                                                                                                             l <- t(matrix(rep((y < pred) - loss.type$tau, nx), ncol = nx)) * 
                                                                                                                                                                                                               matrix(rep(x, npred), ncol = npred)
    }
  } else {
    if (!loss.gradient) {
      if (loss.type$name == "square") 
        l <- (x - y)^2 else if (loss.type$name == "absolute") 
          l <- abs(x - y) else if (loss.type$name == "percentage") 
            l <- abs(x - y)/y else if (loss.type$name == "pinball") 
              l <- ((y < x) - loss.type$tau) * (x - y)
    } else {#gradient trick : we linearize the loss x->l(x) with x->l'(x)*x
      if (loss.type$name == "square") 
        l <- 2 * (pred - y) * x else if (loss.type$name == "absolute") 
          l <- sign(pred - y) * x else if (loss.type$name == "percentage")
            l <- x/y * sign(pred - y) else if (loss.type$name == "pinball") 
              l <- ((y < pred) - loss.type$tau) * x
    }
  }
  return(l)
} 



mixture <- function(Y = NULL, experts = NULL, model = "MLpol", loss.type = "square", 
  loss.gradient = TRUE, coefficients = "Uniform", awake = NULL, parameters = list()) UseMethod("mixture")


mixture.default <- function(Y = NULL, experts = NULL, model = "MLpol", loss.type = "square", 
  loss.gradient = TRUE, coefficients = "Uniform", awake = NULL, parameters = list()) {

  if (!is.list(loss.type)) {
    loss.type <- list(name = loss.type)
  }
  if (!(loss.type$name %in% c("pinball", "square", "percentage", "absolute"))) {
    stop("loss.type should be one of these: 'absolute', 'percentage', 'square', 'pinball'")
  }
  
  
  object <- list(model = model, loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = coefficients, parameters = parameters, Y = NULL, experts = NULL, 
    awake = NULL, training = NULL, names.experts = colnames(experts), T = 0, d = "unknown")
  
  class(object) <- "mixture"
  

  # Test that Y and experts have correct dimensions
  if ((is.null(Y) && !is.null(experts)) || (!is.null(Y) && is.null(experts))) {
    stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
  }
  
  if (!is.null(Y)) {
    
    # Test the dimension of Y: if Y is a matrix, the number of columns is the space of prediction
    if (is.null(dim(Y))) {
      d = 1
      T = length(Y)
    } else {
      d = ncol(Y)
      T = nrow(Y)
      if (d > 1 && T > 1 && length(dim(experts)) < 3) {
        stop("Bad dimensions: nrow(experts) should be equal to dim(experts)[3]")
      } 
      if (length(dim(experts)) == 3) {
        if ((dim(experts)[1] != T) || (dim(experts)[2] != d)){
          stop("Bad dimensions between Y and experts")
        }
      }
      if (T == 1 && d>1) {
        if (length(dim(experts)) == 2) {
          if (dim(experts)[1] != d) {
            stop("Bad dimensions between Y and experts")
          }
        }
      }
    }
    
    if (T == 1 && d == 1) {
      experts <- as.matrix(experts)
      if (nrow(experts) == 1 || ncol(experts) == 1) {
        experts <- matrix(as.numeric(experts), nrow = 1)
      } else {
        stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
      }
    }
    
    if (dim(experts)[1] != T) {
      stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
    }
    object$d <- d
    object <- predict(object, newY = Y, newexperts = experts, awake = awake, 
      type = "model")
    
  }
  return(object)
} 

predict.mixture <- function(object, newexperts = NULL, newY = NULL, awake = NULL, 
                            online = TRUE, type = c("model", "response", "weights", "all"), ...) {
  result <- object
  d <- object$d
  if ((d == 1) || (d == "unknown" && is.null(dim(newY)))) {
    object$d <- 1    
    return(predictReal(object, newexperts, newY, awake, 
                online, type, ...))
  } else {
    if (d == "unknown") {
      d = dim(newY)[2]
      T = dim(newY)[1]
      # Bad dimension for experts
      if (T > 1 && length(dim(newexperts)) < 3) {
        stop("Bad dimensions: nrow(experts) should be equal to dim(experts)[3]")
      } 
      if (length(dim(newexperts)) == 3) {
        if ((dim(newexperts)[1] != T) || (dim(newexperts)[2] != d)){
          stop("Bad dimensions between Y and experts")
        }
      }      
      if (T == 1) {
        if (length(dim(newexperts)) == 2) {
          if (dim(newexperts)[1] != d) {
            stop("Bad dimensions between Y and experts")
          } else {
            newexperts = array(newexperts, dim = c(1,dim(newexperts)))
          }
        }
      }
    }
    result$d <- d
    awakei <- NULL
    for (i in 1:nrow(newY)) {
      if (!online){
        stop("Batch prediction are currently not supported for dimension > 1")
      }
      if (!is.null(awake)){
         awakei <- as.matrix(awake[i,,])
         #awakei <- as.matrix(awake[i])
      }
      result <- predictReal(result, newexperts = as.matrix(newexperts[i,,]), newY = c(newY[i,]), awake = awakei, 
                            online = FALSE, type, ...)
    }
  }
  result$weights <- matrix(result$weights, nrow = result$T)
  
  return(result)
}

predictReal <- function(object, newexperts = NULL, newY = NULL, awake = NULL, 
                            online = TRUE, type = c("model", "response", "weights", "all"), ...) {

  type <- match.arg(type)

  # Number of instant and number of experts
  if (!is.null(newY)) {
    T <- length(newY)
    if (is.null(object$names.experts)) {
      if (T == 1) {
        object$names.experts <- names(newexperts)
      } else {
        object$names.experts <- colnames(newexperts)
      }
    }
    newexperts <- matrix(newexperts, nrow = T)
    N <- ncol(newexperts)
  } else if (object$coefficients[1] != "Uniform") {
    N <- length(object$coefficients)
    newexperts <- matrix(newexperts, ncol = N)
    T <- nrow(newexperts)
  } else {
    warning("You should provide observations to train non trivial model")
    N <- ncol(newexperts)
    T <- nrow(newexperts)
    
    if (is.null(newexperts)) {
      result <- switch(type, model = object, response = NULL, weights = NULL, 
                       all = list(model = object, response = NULL, weights = NULL))
      return(result)
    }
  }

  if (!is.null(awake)) {
    awake <- matrix(awake, nrow = T)
    if (!identical(dim(awake), dim(newexperts))) {
      stop("Bad dimensions: awake and newexperts should have same dimensions")
    }
  } else {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  if ((object$model=="OGD") && (object$model=="ridge") && (object$model=="ridgeCalib")) {#pour pouvoir gerer les obs manquantes dans les agregations
    idx.na <- which(is.na(newexperts))
    awake[idx.na] <- 0
    newexperts[idx.na] <- 0
  }

  # test possible warnings and errors
  if (is.null(object$training) && (object$coefficients != "Uniform") && (object$model == 
                                                                         "MLpol")) {
    stop(paste(object$model, "cannot handle non-uniform prior weight vector"))
  }
  
  if (object$coefficients[1] == "Uniform") {
    object$coefficients <- rep(1/N, N)
  }
  
  if (length(object$coefficients) != N) {
    stop("Bad number of experts: (length(object$coefficients) != nrow(newexperts))")
  }
  
  if (!is.null(object$loss.type$tau) && object$loss.type$name != "pinball") {
    warning("Unused parameter tau (loss.type$name != 'pinball)")
  }
  
  if (object$loss.type$name != "square" && object$model == "Ridge") {
    stop(paste("Square loss is require for", object$model, "model."))
  }
  
  if (!is.null(awake) && !identical(awake, matrix(1, nrow = T, ncol = N)) && (object$model == 
                                                                              "Ridge" || object$model == "OGD")) {
    stop(paste("Sleeping or missing values not allowed for", object$model, "model."))
  }

  # if no expert advice is provided, it returns the fitted object
  if (is.null(newexperts)) {
    if (!is.null(newY)) {
      stop("Expert advice should be provided if newY is non null")
    }
    result <- switch(type, model = object, response = object$prediction, weights = object$weights, 
                     all = list(model = object, response = object$prediction, weights = object$weights))
    return(result)
  }
  
  if (!is.null(newexperts) && is.null(newY) && online) {
    stop("newY cannot be null to perform online prediction. Provide newY or set online = FALSE")
  }
  
  newexperts <- matrix(as.numeric(as.matrix(newexperts)), nrow = T)
  colnames(newexperts)<-object$names.experts

  # Batch prediction and weights
  if (!online) {
    w <- matrix(object$coefficients, ncol = 1)
    pond <- c(awake %*% w)
    newpred <- c(((newexperts * awake) %*% w)/pond)
    newweights <- (t(t(awake) * c(w)) / pond)[seq(1,T,by=object$d),]
  }  

  # Online prediction and weights if newY is provided
  if (!is.null(newY)) {
    
    if (min(newY) <= 0 && object$loss.type$name == "percentage") {
      stop("Y should be non-negative for percentage loss function")
    }
    
    ## If averaged is true, the models do not use coefficient as the next weight
    if (!is.null(object$parameters$averaged) && object$parameters$averaged && !is.null(object$training)) {
      object$coefficients <- object$training$next.weights
    }
    

    if (object$model == "MLpol") {
      newobject <- MLpol(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
                         loss.gradient = object$loss.gradient, training = object$training)
      newobject$parameters <- list(eta = rbind(object$parameters$eta, newobject$parameters$eta))
    }
    
    if (object$model == "BOA") {
      algo <- eval(parse(text = 'object$model'))
      newobject <- eval(parse(text=paste(algo,"(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
                        loss.gradient = object$loss.gradient, w0 = object$coefficients, training = object$training)",sep="")))
      newobject$parameters <- list(eta = rbind(object$parameters$eta, newobject$parameters$eta))
    }
    
    if (object$model == "EWA") {
      if (is.null(object$parameters$eta) || !is.null(object$parameters$grid.eta)) {
        if (is.null(object$parameters$grid.eta)) {
          object$parameters$grid.eta <- 1
        }
        
        newobject <- ewaCalib(y = newY, experts = newexperts, awake = awake, 
                              loss.type = object$loss.type, loss.gradient = object$loss.gradient, 
                              w0 = object$coefficients, gamma = object$parameters$gamma, grid.eta = sort(object$parameters$grid.eta), 
                              training = object$training)
        newobject$parameters$eta <- c(object$parameters$eta, newobject$parameters$eta)
      } else {
        newobject <- ewa(y = newY, experts = newexperts, eta = object$parameters$eta, 
                         awake = awake, loss.type = object$loss.type, loss.gradient = object$loss.gradient, 
                         w0 = object$coefficients, training = object$training)
      }
    }

    if (object$model == "NBM") {
      newobject <- nbm(y = newY, experts = newexperts, eta = object$parameters$eta, 
                        awake = awake, loss.type = object$loss.type, 
                        w0 = object$coefficients, training = object$training)
    }

    newobject$Y <- rbind(object$Y, matrix(newY, ncol = object$d))
    newobject$experts <- rbind(object$experts, newexperts)
    newobject$names.experts <- object$names.experts
    if (is.null(newobject$names.experts)) {
      if (!is.null(colnames(newexperts))) {
        newobject$names.experts <- colnames(newexperts)
      } else {
        if (!is.null(names(newexperts))) {
          newobject$names.experts <- names(newexperts)
        } else {
          newobject$names.experts <- paste("X",1:N,sep="")
        }
      }
    }
    newobject$awake <- rbind(object$awake, awake)

    colnames(newobject$experts) <- newobject$names.experts
    colnames(newobject$weights) <- newobject$names.experts
    colnames(newobject$awake) <- newobject$names.experts
    
    # Averaging of the weights if asked by the averaged parameter
    if (is.null(object$parameters$averaged)) {
      newobject$parameters$averaged = FALSE
    } else {
      newobject$parameters$averaged = object$parameters$averaged
    }
    
    newobject$weights = newobject$weights[seq(1,T,by=object$d),]
    
    if (newobject$parameters$averaged) {
      if (object$T == 0) {
        newweights.avg <- apply(newobject$weights, 2, cumsum) / (1:(T/object$d))
      } else {
        newweights.avg <- (object$training$sumweights + apply(newobject$weights, 2, cumsum)) / (object$T + 1:(T/object$d))
      }
      
      newobject$training$sumweights <- (object$T + T/object$d) * newweights.avg[T/object$d,] + newobject$coefficients
      newobject$training$next.weights <- newobject$coefficients
      newobject$coefficients <- newobject$training$sumweights / (object$T + T/object$d + 1)
    }
    
    # If online is true, we use online prediction
    if (online) {
      if (newobject$parameters$averaged) {
        newweights <- newweights.avg
        newpred <- apply(newweights.avg * newexperts,1,sum)
      } else {
        newweights <- newobject$weights
        newpred <- newobject$prediction
      }
    }
    newobject$prediction <- rbind(object$prediction, matrix(newpred, ncol = object$d))
    newobject$weights <- rbind(object$weights, newweights)
    rownames(newobject$weights) <- NULL
    newobject$loss <- mean(loss(c(newobject$prediction), c(newobject$Y), loss.type = newobject$loss.type))
    newobject$T <- object$T + T/object$d
    newobject$d <- object$d
  } else {
    newobject <- object
  }
  class(newobject) <- "mixture"
  
  
  
  result <- switch(type, model = newobject, response = matrix(newpred, ncol = object$d), weights = newweights, 
                   all = list(model = newobject, response = newpred, weights = newweights))
  
  return(result)
} 

ewaCalib <- function(y, experts, grid.eta = 1, awake = NULL, loss.type = "square", 
  loss.gradient = TRUE, w0 = NULL, trace = F, gamma = 2, training = NULL) {
  experts <- as.matrix(experts)

  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  if (is.null(gamma)) {
    gamma <- 2
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    y[which(is.na(y))]<- 0
    awake[idx.na] <- 1 #to avoid 0 division
    experts[idx.na] <- 0
  }
  
  T0 <- 0  # number of observations in previous runs
  neta <- length(grid.eta)  # Initial number of learning parameter in the grid to be optimized
  besteta <- floor(neta)/2 + 1  # We start with the parameter eta in the middle of the grid
  eta <- rep(grid.eta[besteta], T)  # Vector of calibrated learning rates (will be filled online by the algorithm)
  cumulativeLoss <- rep(0, neta)  # Cumulative loss suffered by each learning rate
  iLoss <- rep(0, neta)

  R <- array(0, c(N, neta))  # Matrix of regret suffered by each learning rate (columns) against each expert (rows)
  weta <- array(w0, dim = c(N, neta))  # Weight matrix assigned by each algorithm EWA(eta) 
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Predictions formed by the mixture
  
  if (!is.null(training)) {
    weta <- training$weta[[length(training$weta)]]
    w0 <- training$w0
    R <- training$R
    cumulativeLoss <- training$cumulativeLoss[nrow(training$cumulativeLoss)]
    besteta <- training$besteta
    eta[1] <- grid.eta[besteta]
    T0 <- training$T
  }
  
  R.w0 <- R
  for (k in 1:neta) {
    R.w0[, k] <- R[, k] + log(w0)/grid.eta[k]
  }  # We initialize the regrets so that w0 is the initial weight vector
  
  
  for (t in 1:T) {
    # Display the state of progress of the algorithm
    if (!(t%%floor(T/10)) && trace) 
      cat(floor(10 * t/T) * 10, "% -- ")
    
    # Weights, prediction formed by EWA(eta[t]) where eta[t] is the learning rate
    # calibrated online
    weights[t, ] <- weta[, besteta] * awake[t, ]/sum(weta[, besteta] * awake[t, 
      ])
    prediction[t] <- experts[t, ] %*% weights[t, ]
    eta[t] <- grid.eta[besteta]
    
    # Weights, predictions formed by each EWA(eta) for eta in the grid 'grid.eta'
    pred <- experts[t, ] %*% t(t(weta * awake[t, ])/apply(weta * awake[t, ], 
      2, sum))
    iLoss <- loss(pred, y[t], loss.type)
    cumulativeLoss <- cumulativeLoss + iLoss  # cumulative loss without gradient trick
    if (neta == 1){
      lpred <- lossPred(pred, y[t], pred, loss.type, loss.gradient)
    } else {
      lpred <- diag(lossPred(pred, y[t], pred, loss.type, loss.gradient))  # gradient loss suffered by each eta on the grid
    }
    lexp <- lossPred(experts[t, ], y[t], pred, loss.type, loss.gradient)  # gradient loss suffered by each expert
    
    # Regret update
    R.w0 <- R.w0 + awake[t, ] * t(c(lpred) - t(lexp))
    
    # Update of the best parameter
    besteta <- order(cumulativeLoss)[1]
 
    weta <- truncate1(exp(t(t(matrix(R.w0, ncol = neta)) * grid.eta)))
  }
  # Next weights
  w <- weta[, besteta]/sum(weta[, besteta])
  
  R <- R.w0
  for (k in 1:neta) {
    R[, k] <- R.w0[, k] - log(w0)/grid.eta[k]
  }  
  
  object <- list(model = "EWA", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)
  
  object$parameters <- list(eta = eta[1:T], grid.eta = grid.eta)
  object$weights <- weights
  object$prediction <- prediction

  object$training <- list(T = T0 + T, weta = c(training$weta,list(weta)), w0 = w0, R = R, besteta = besteta,
      grid.loss = cumulativeLoss/(T0 + T), oldexperts = rbind(training$oldexperts, experts),
      oldY = c(training$oldY, y), oldawake = rbind(training$oldawake,awake), iLoss = rbind(training$iLoss, iLoss),
      cumulativeLoss = rbind(training$cumulativeLoss, cumulativeLoss) )

  class(object) <- "mixture"
  if (trace) 
    cat("\n")
  return(object)
}




fixedshareBoaCalib <- function(y, experts, grid.alpha = 10^(-4:-1), awake = NULL, 
  loss.type = "square", loss.gradient = TRUE, w0 = NULL, trace = F, gamma = 2, 
  training = NULL) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    y[which(is.na(y))]<- 0
    awake[idx.na] <- 1 #to avoid 0 division
    experts[idx.na] <- 0
  }
  
  T0 <- 0  # number of observations in previous runs
  nalpha <- length(grid.alpha)  # Number of mixing rates in the grid
  bestpar <- c(floor(nalpha)/2) + 1  # We start with the parameters in the middle of the grids
  parEta <- array(0, dim = c(0,N) )
  parAlpha <- c()
  
  cumulativeLoss <- array(0, dim = c(nalpha))
  iLoss <- array(0, dim = c(nalpha))
  
  R.reg <- array(0, dim = c(N, nalpha))
  R <- array(0, dim = c(N, nalpha))
  B <- array(0, dim = c(N, nalpha))
  V <- array(0, dim = c(N, nalpha))
  r <- array(0, dim = c(N, nalpha))
  r.reg <- array(0, dim = c(N, nalpha))
  eta <- array(exp(350), dim=c(T+1, N, nalpha))
  
  wpar <- array(w0, dim = c(N, nalpha))  # Weight array (in 3 dimensions) assigned by each algorithm FixedShare(eta,alpha)
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Predictions formed by the mixture
  
  if (!is.null(training)) {
    bestpar <- training$bestpar
    w0 <- training$w0
    R <- training$R
    R.reg <- training$R.reg #TODO verifier bonnes dimensions (N*nalpha)
    wpar <- training$wpar
    eta[1, ,] <- training$eta
    B <- training$B
    V <- training$V
    cumulativeLoss <- training$cumulativeLoss
    T0 <- training$T
  }

  for (t in 1:T) {
    # Display the state of progress of the algorithm
    if (!(t%%floor(T/10)) && trace) 
      cat(floor(10 * t/T) * 10, "% -- ")
    # Weights, prediction forme by FixedShare(eta[t],alpha[t]) where par[t,] =
    # c(eta[t],alpha[t]) are the parameters calibrated online
    weights[t, ] <- awake[t, ]*wpar[, bestpar]/sum(awake[t, ]*wpar[, bestpar])
    prediction[t] <- experts[t, ] %*% weights[t, ]
    parEta <- rbind(parEta, eta[t, , bestpar])
    parAlpha <- c(parAlpha, grid.alpha[bestpar])

    # Loop over the mixing rates alpha in the grid 'grid.alpha'
    for (k in 1:nalpha) {
      # Weights, prediction, and losses formed by FixedShare for grid.alpha[k]
      # in the grid
      # waux <- t(t(t(truncate1(exp(log(eta[t, , k]*w0) + eta[t, , k] * R.reg[, k])) * t(awake[t, ])))
      #           /apply(as.matrix(t(truncate1(exp(log(eta[t, , k]*w0) + eta[t, , k] * R.reg[, k])) * t(awake[t, ]))), 2, sum))
      waux <- awake[t, ]*wpar[, k]/sum(awake[t, ]*wpar[, k])
      pred <- experts[t, ] %*% waux

      iLoss[k] <- loss(pred, y[t], loss.type)
      cumulativeLoss[k] <- cumulativeLoss[k] + iLoss[k]  # Non gradient cumulative losses
      lpred <- lossPred(pred, y[t], pred, loss.type, loss.gradient)#TODO verifier que tout ok dimension (cf FSCALIB normal avec if neta=1)
      lexp <- lossPred(experts[t, ], y[t], pred, loss.type, loss.gradient)

      # Vector, instantaneous regret
      r[,k] <-  awake[t, ] * (lpred - lexp) # -l_{j,t} in wintenberger
      # Update the learning rates
      B[, k] <- pmax(B[, k], abs(r[,k])) #vector, estimation of the loss bound for the M experts
      V[, k] <- V[, k] + r[,k]^2 #vector, cumulated quadratic loss for each experts

      eta[t + 1, , k] <- pmin(pmin(1/B[, k], sqrt(log(1/w0)/V[, k])),exp(350)) #vector, learning rate for the M experts

      r.reg[,k] <- 1/2 * (r[,k] - eta[t + 1, , k] * r[,k]^2 + B[, k] * (eta[t + 1, , k] * r[,k] > 1/2))

      # Update the regret and the regularized regret used by BOA
      #R <- t(t(log(wpar[, k]))/eta[t + 1, , k]) + awake[t, ] * t(c(lpred) - t(lexp))
      R[, k] <- R[, k] + r[,k]
      R.reg[, k] <- R.reg[, k] + r.reg[,k] # TODO est-ce qu'il ne faudrait pas multiplier r.reg par awake comme dans FSCALIB normal ?
      
      # Weight update
      v <- truncate1(exp(log(eta[t+1, , k]*w0) + eta[t+1, , k] * R.reg[, k]))
      v <- v/sum(v)  # Renormalization

      wpar[, k] <- grid.alpha[k]/N + (1 - grid.alpha[k]) * v
    }
    # Grid update *************** find the best parameter in the grid
    bestpar <- which(cumulativeLoss == min(cumulativeLoss), arr.ind = TRUE)[1, ]
  }

  # Next weights
  w <- wpar[, bestpar]/sum(wpar[, bestpar])

  object <- list(model = "FSBOA", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)

  object$parameters <- list(eta = parEta, alpha = parAlpha, grid.alpha = grid.alpha)
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(T = T0 + T, wpar = wpar, w0 = w0, R=R, R.reg = R.reg, V=V, B=B,
    cumulativeLoss = cumulativeLoss, iLoss = rbind(training$iLoss,iLoss), bestpar = bestpar, grid.loss = cumulativeLoss/(T0 + T),
      eta = eta[T + 1, ,],r=c(training$r,list(r)),r.reg=c(training$r.reg,list(r.reg)))
  
  class(object) <- "mixture"
  
  if (trace) 
    cat("\n")
  return(object)
}




fixedshareSOCOBoaCalib <- function(y, experts, grid.alpha = 10^(-4:-1), awake = NULL, 
  loss.type = "square", loss.gradient = TRUE, w0 = NULL, trace = F, gamma = 2, 
  training = NULL) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    y[which(is.na(y))]<- 0
    awake[idx.na] <- 1 #to avoid 0 division
    experts[idx.na] <- 0
  }
  
  T0 <- 0  # number of observations in previous runs
  nalpha <- length(grid.alpha)  # Number of mixing rates in the grid
  bestpar <- c(floor(nalpha)/2) + 1  # We start with the parameters in the middle of the grids
  parEta <- array(0, dim = c(0,N) )
  parAlpha <- c()
  
  cumulativeLoss <- array(0, dim = c(nalpha))
  iLoss <- array(0, dim = c(nalpha))
  
  R.reg <- array(0, dim = c(N, nalpha))
  R <- array(0, dim = c(N, nalpha))
  B <- array(0, dim = c(N, nalpha))
  V <- array(0, dim = c(N, nalpha))
  r <- array(0, dim = c(N, nalpha))
  r.reg <- array(0, dim = c(N, nalpha))
  eta <- array(exp(350), dim=c(T+1, N, nalpha))
  eta_inv_2 <- array(exp(-350*2), dim=c(N, nalpha))
  
  wpar <- array(w0, dim = c(N, nalpha))  # Weight array (in 3 dimensions) assigned by each algorithm FixedShare(eta,alpha)
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Predictions formed by the mixture
  
  if (!is.null(training)) {
    bestpar <- training$bestpar
    w0 <- training$w0
    R <- training$R
    R.reg <- training$R.reg #TODO verifier bonnes dimensions (N*nalpha)
    wpar <- training$wpar
    eta[1, ,] <- training$eta
    B <- training$B
    V <- training$V
    cumulativeLoss <- training$cumulativeLoss
    T0 <- training$T
  }

  for (t in 1:T) {
    # Display the state of progress of the algorithm
    if (!(t%%floor(T/10)) && trace) 
      cat(floor(10 * t/T) * 10, "% -- ")
    # Weights, prediction forme by FixedShare(eta[t],alpha[t]) where par[t,] =
    # c(eta[t],alpha[t]) are the parameters calibrated online
    weights[t, ] <- awake[t, ]*wpar[, bestpar]/sum(awake[t, ]*wpar[, bestpar])
    prediction[t] <- experts[t, ] %*% weights[t, ]
    parEta <- rbind(parEta, eta[t, , bestpar])
    parAlpha <- c(parAlpha, grid.alpha[bestpar])

    # Loop over the mixing rates alpha in the grid 'grid.alpha'
    for (k in 1:nalpha) {
      waux <- awake[t, ]*wpar[, k]/sum(awake[t, ]*wpar[, k])
      pred <- experts[t, ] %*% waux

      iLoss[k] <- loss(pred, y[t], loss.type)
      cumulativeLoss[k] <- cumulativeLoss[k] + iLoss[k]  # Non gradient cumulative losses
      lpred <- lossPred(pred, y[t], pred, loss.type, loss.gradient)#TODO verifier que tout ok dimension (cf FSCALIB normal avec if neta=1)
      lexp <- lossPred(experts[t, ], y[t], pred, loss.type, loss.gradient)

      # Vector, instantaneous regret
      r[,k] <-  awake[t, ] * (lpred - lexp) # -l_{j,t} in wintenberger

      # Update the learning rates
      eta_inv_2[,k] = 1/(eta[t, ,k]*eta[t, ,k]) + r[,k]^2 / log(1/w0)
      eta[t + 1, ,k] <- pmin(sqrt(1/eta_inv_2[t+1, ,k]),exp(350)) #vector, learning rate for the M experts

      # Update the regret and the regularized regret used by BOA
      r.reg[,k] <- r[,k] - eta[t + 1, , k] * r[,k]^2
      R[, k] <- R[, k] + r[,k]
      R.reg[, k] <- R.reg[, k] + r.reg[,k] # TODO est-ce qu'il ne faudrait pas multiplier r.reg par awake comme dans FSCALIB normal ?
      
      # Weight update
      v <- truncate1(exp(log(eta[t+1, , k]*w0) + eta[t+1, , k] * R.reg[, k]))
      v <- v/sum(v)  # Renormalization

      wpar[, k] <- grid.alpha[k]/N + (1 - grid.alpha[k]) * v
    }
    # Grid update *************** find the best parameter in the grid
    bestpar <- which(cumulativeLoss == min(cumulativeLoss), arr.ind = TRUE)[1, ]
  }

  # Next weights
  w <- wpar[, bestpar]/sum(wpar[, bestpar])

  object <- list(model = "FSSOCOBOA", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)

  object$parameters <- list(eta = parEta, alpha = parAlpha, grid.alpha = grid.alpha)
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(T = T0 + T, wpar = wpar, w0 = w0, R=R, R.reg = R.reg, V=V, B=B,
    cumulativeLoss = cumulativeLoss, iLoss = rbind(training$iLoss,iLoss), bestpar = bestpar, grid.loss = cumulativeLoss/(T0 + T),
      eta = eta[T + 1, ,],r=c(training$r,list(r)),r.reg=c(training$r.reg,list(r.reg)))
  
  class(object) <- "mixture"
  
  if (trace) 
    cat("\n")
  return(object)
}




originalfixedshareBoaCalib <- function(y, experts, grid.alpha = 10^(-4:-1), awake = NULL, 
  loss.type = "square", loss.gradient = TRUE, w0 = NULL, trace = F, gamma = 2, 
  training = NULL) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  if (is.null(gamma)) {
    gamma <- 2
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    y[which(is.na(y))]<- 0
    awake[idx.na] <- 1 #to avoid 0 division
    experts[idx.na] <- 0
  }
  
  T0 <- 0  # number of observations in previous runs
  nalpha <- length(grid.alpha)  # Number of mixing rates in the grid
  bestpar <- c(floor(nalpha)/2) + 1  # We start with the parameters in the middle of the grids
  parEta <- NULL
  parAlpha <- NULL
  
  cumulativeLoss <- array(0, dim = c(nalpha))
  iLoss <- array(0, dim = c(nalpha))
  
  R.reg <- array(0, dim = c(N, nalpha))
  B <- array(0, dim = c(N, nalpha))
  V <- array(0, dim = c(N, nalpha))
  r <- array(0, dim = c(N, nalpha))
  r.reg <- array(0, dim = c(N, nalpha))
  lpred <- array(0, dim = c(nalpha))
  eta <- array(exp(350), dim=c(T+1, N, nalpha))
  
  wpar <- array(w0, dim = c(N, nalpha))  # Weight array (in 3 dimensions) assigned by each algorithm FixedShare(eta,alpha)
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Predictions formed by the mixture
  
  if (!is.null(training)) {
    bestpar <- training$bestpar
    wpar <- training$wpar[[length(training$wpar)]]
    w0 <- training$w0
    R.reg <- training$R.reg
    B <- training$B
    V <- training$V
    cumulativeLoss <- training$cumulativeLoss
    T0 <- training$T
    #eta[1, ,] <- training$eta[, ]
  }

  for (t in 1:T) {
    # Display the state of progress of the algorithm
    if (!(t%%floor(T/10)) && trace) 
      cat(floor(10 * t/T) * 10, "% -- ")
    # Weights, prediction forme by FixedShare(eta[t],alpha[t]) where par[t,] =
    # c(eta[t],alpha[t]) are the parameters calibrated online
    weights[t, ] <- wpar[, bestpar] * awake[t, ]/sum(wpar[,bestpar] * awake[t, ])
    prediction[t] <- experts[t, ] %*% weights[t, ]
    parEta <- rbind(parEta, eta[t, , bestpar])
    parAlpha <- c(parAlpha, grid.alpha[bestpar])
    
    # Loop over the mixing rates alpha in the grid 'grid.alpha'
    for (k in 1:nalpha) {
      
      # Weights, prediction, and losses formed by FixedShare(eta,alpha) for (eta,alpha)
      # in the grid
      waux <- t(t(wpar[, k] * awake[t, ])/apply(as.matrix(wpar[, k] * awake[t, 
        ]), 2, sum))
      pred <- experts[t, ] %*% waux

      iLoss[k] <- loss(pred, y[t], loss.type)
      cumulativeLoss[k] <- cumulativeLoss[k] + iLoss[k]  # Non gradient cumulative losses
      lpred[k] <- lossPred(pred, y[t], pred, loss.type, loss.gradient)
      lexp <- lossPred(experts[t, ], y[t], pred, loss.type, loss.gradient)

      # Vector, instantaneous regret
      r[,k] <-  awake[t, ] * (lpred[k] - lexp) # -l_{j,t} in wintenberger
      # Update the learning rates
      B[, k] <- pmax(B[, k], abs(r[,k])) #vector, estimation of the loss bound for the M experts
      V[, k] <- V[, k] + r[,k]^2 #vector, cumulated quadratic loss for each experts  
      eta[t + 1, , k] <- pmin(pmin(1/B[, k], sqrt(log(1/w0)/V[, k])),exp(350)) #vector, learning rate for the M experts
      if (max(eta[t + 1, , k]) > exp(300)) {
        # if some losses still have not been observed, exp(300)=>loss~0
        r.reg[,k] <- r[,k] + B[,k]/2
      } else {
        r.reg[,k] <- (1/2)*(r[,k] - eta[t + 1, , k] * r[,k]^2)
      }

      # Update the regret and the regularized regret used by BOA
      #R <- t(t(log(wpar[, k]))/eta[t + 1, , k]) + awake[t, ] * t(c(lpred) - t(lexp))
      R.reg[, k] <- t(t(log(wpar[, k]))/eta[t + 1, , k]) + awake[t, ]*r.reg[,k]
      
      # Weight update
      v <- truncate1(eta[t + 1, , k]*exp(log(w0) + eta[t + 1, , k] * R.reg[, k]))
      #v <- truncate1(exp(t(t(matrix(R.reg, ncol = neta)) * grid.eta)))
      v <- v/sum(v)  # Renormalization of each column
      wpar[, k] <- grid.alpha[k]/N + (1 - grid.alpha[k]) * v
    }

    # Grid update *************** find the best parameter in the grid
    bestpar <- which(cumulativeLoss == min(cumulativeLoss), arr.ind = TRUE)[1, ]
  }
  
  # Next weights
  w <- wpar[, bestpar]/sum(wpar[, bestpar])

  object <- list(model = "FSBOA", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)

  object$parameters <- list(eta = parEta, alpha = parAlpha, grid.alpha = grid.alpha)
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(T = T0 + T, wpar = c(training$wpar,list(wpar)), w0 = w0, R.reg = R.reg, V=V, B=B,
    cumulativeLoss = cumulativeLoss, iLoss = rbind(training$iLoss,iLoss), bestpar = bestpar, grid.loss = cumulativeLoss/(T0 + T),
      oldexperts = rbind(training$oldexperts, experts), oldY = c(training$oldY, y), oldawake = rbind(training$oldawake, 
        awake), eta = eta[T + 1, , ],r=c(training$r,list(r)),r.reg=c(training$r.reg,list(r.reg)))
  
  class(object) <- "mixture"
  
  if (trace) 
    cat("\n")
  return(object)
}



# Ridge aggregation rule with automatic calibration of smoothing parameters
ridgeCalib <- function(y, experts, grid.lambda = 1, w0 = NULL, trace = FALSE, gamma = 2, 
                       training = NULL) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform intial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- matrix(1/N, ncol = N)
  }
  if (is.null(grid.lambda)) {
    grid.lambda <- 1
  }
  if (is.null(gamma)) {
    gamma <- 2
  }
  
  # Smoothing parameter grid
  nlambda <- length(grid.lambda)
  grid.lambda <- matrix(grid.lambda, nrow = nlambda)
  cumulativeLoss <- rep(0, nlambda)
  
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Vector of predictions formed by the mixing algorithm
  pred.lambda <- matrix(0, ncol = nlambda, nrow = T)  # Prediction of mixture algorithm with different learning rates eta
  
  if (!is.null(training)) {
    At <- training$At
    bt <- training$bt
    bestlambda <- training$bestlambda
    wlambda <- training$wlambda
    w0 <- training$w0
    cumulativeLoss <- training$cumulativeLoss
    T0 <- training$T
  } else {
    At <- diag(0, N)
    bt <- rep(0, N)
    bestlambda <- floor(nlambda)/2 + 1  # We start with the parameter in the middle of the grid
    wlambda <- array(w0, dim = c(N, nlambda))  # Weight matrix proposed by each Ridge(lambda) where lambda is a parameter of the grid
    T0 <- 0
  }
  
  lambda <- rep(grid.lambda[bestlambda], T)
  for (t in 1:T) {
    # Display the state of progress of the algorithm
    if (!(t%%floor(T/10)) && trace) 
      cat(floor(10 * t/T) * 10, "% -- ")
    
    # Weights, prediction forme by Ridge(lambda[t]) where lambda[t] is the learning
    # rate calibrated online
    weights[t, ] <- wlambda[, bestlambda]
    prediction[t] <- experts[t, ] %*% weights[t, ]
    lambda[t] <- grid.lambda[bestlambda]
    
    # Weights, predictions formed by each Ridge(lambda) for lambda in the grid
    # 'grid.lambda'
    pred.lambda[t, ] <- experts[t, ] %*% wlambda
    cumulativeLoss <- cumulativeLoss + (pred.lambda[t, ] - y[t])^2
    
    # Parameter update
    At <- At + experts[t, ] %*% t(experts[t, ])
    bt <- bt + y[t] * experts[t, ]
    
    # Grid update **************
    bestlambda <- order(cumulativeLoss)[1]  # find the best smoothing rate lambda in the grid
    
    # Expand the grid if the best parameter lies on an extremity
    if (bestlambda == nlambda) {
      if (trace) 
        cat(" + ")
      newlambda <- grid.lambda[bestlambda] * gamma^(1:3)
      grid.lambda <- c(grid.lambda, newlambda)
      nlambda <- nlambda + length(newlambda)
      for (k in 1:length(newlambda)) {
        perfnewlambda <- tryCatch(ridge(y = c(training$oldY, y[1:t]), experts = rbind(training$oldexperts, 
                                                                                      matrix(experts[1:t, ], ncol = N)), lambda = newlambda[k], w0 = w0), 
                                  error = function(e) {
                                    list(prediction = rep(0, t))
                                  })
        newcumulativeLoss <- sum((perfnewlambda$prediction - y[1:t])^2)
        cumulativeLoss <- c(cumulativeLoss, newcumulativeLoss)
        pred.lambda <- cbind(pred.lambda, c(perfnewlambda$prediction, rep(0, 
                                                                          (T - t))))
      }
    }
    if (bestlambda == 1) {
      if (trace) 
        cat(" - ")
      newlambda <- grid.lambda[bestlambda]/gamma^(1:3)
      nlambda <- nlambda + length(newlambda)
      bestlambda <- bestlambda + length(newlambda)
      for (k in 1:length(newlambda)) {
        grid.lambda <- c(newlambda[k], grid.lambda)
        perfnewlambda <- tryCatch(y = ridge(c(training$oldY, y[1:t]), experts = rbind(training$oldexperts, 
                                                                                      matrix(experts[1:t, ], ncol = N)), lambda = newlambda[k], w0 = w0), 
                                  error = function(e) {
                                    list(prediction = rep(NA, t))
                                  })
        newcumulativeLoss <- sum((perfnewlambda$prediction - y[1:t])^2)
        cumulativeLoss <- c(newcumulativeLoss, cumulativeLoss)
        pred.lambda <- cbind(c(perfnewlambda$prediction, rep(NA, (T - t))), 
                             pred.lambda)
      }
    }
    wlambda <- matrix(0, nrow = N, ncol = nlambda)
    for (k in 1:nlambda) {
      wlambda[, k] <- tryCatch(solve(grid.lambda[k] * diag(1, N) + At, matrix(grid.lambda[k] * 
                                                                                w0, nrow = N) + bt), error = function(e) {
                                                                                  NA
                                                                                })
    }
    # the smoothing parameter has to big large enough in order to have invertible
    # design matrix
    lambda.min <- which(!(is.na(wlambda[1, ])))[1]
    bestlambda <- max(lambda.min, bestlambda)
  }
  
  object <- list(model = "Ridge", loss.type = list(name = "square"), coefficients = wlambda[, 
                                                                                            bestlambda])
  
  object$parameters <- list(lambda = c(lambda[1:T]), grid.lambda = c(grid.lambda))
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(T = T0 + T, wlambda = wlambda, w0 = w0, At = At, bt = bt, 
                          bestlambda = bestlambda, cumulativeLoss = cumulativeLoss, grid.loss = cumulativeLoss/(T0 + 
                                                                                                                  T), oldexperts = rbind(training$oldexperts, experts), oldY = c(training$oldY, 
                                                                                                                                                                                 y))
  
  if (trace) 
    cat("\n")
  return(object)
} 


print.mixture <- function(x, ...) {
  cat("Aggregation rule: ")
  cat(x$model, "\n")
  cat("Loss function: ", x$loss.type$name, "loss", "\n")
  cat("Gradient trick: ", x$loss.gradient, "\n")
  cat("Coefficients: ")
  if (x$coefficients[1] != "Uniform") {
    cat("\n")
    x$coefficients <- data.frame(signif(matrix(as.matrix(x$coefficients), nrow = 1)))
    names(x$coefficients) <- colnames(x$experts)
    rownames(x$coefficients) <- ""
    print(signif(x$coefficients, digits = 3))
  } else {
    print("Uniform mixture")
  }
}

seriesToBlock <- function(X, d) {
  f <- function(Y){
    matrix(Y, byrow = TRUE, ncol = d)
  }
  if (is.null(dim(X)) || dim(X)[2] == 1) {
    if ((length(X) %% d) != 0) {
      stop("length(X) must be a multiple of d")
    }
    return(f(X))
  } else {
    n <- dim(X)[1]
    K <- dim(X)[2]
    if ((n %% d) != 0) {
      stop("dim(X)[1] should be a multiple of d")
    }
    M <- array(apply(X,2,f),dim = c(n/d,d,K))
    return(M)
  }
}


blockToSeries <- function(X){
  if (is.null(dim(X)) || length(dim(X)) > 3){
    stop("X must be of dimension 2 or 3")
  } 
  # if X is a matrix, we convert it to a vector 
  if (length(dim(X)) == 2) {
    return(c(t(X)))
  } 
  if (length(dim(X)) == 3) {
    return(array(c(aperm(X,c(2,1,3))),dim = c(dim(X)[1]*dim(X)[2],dim(X)[3])))
  }
}
simplexProj <- function(u, U = 1) {
  N <- length(u)
  v = sort(u,decreasing = TRUE)
  rho = max(which((v + (U-cumsum(v)) / (1:N)) >0))
  if (rho >= 1) {
    lambda = 1/rho * (U - sum(v[1:rho]))
    return(pmax(u+ lambda,0))
  } else {
    return(rep(1/N,N))
  }
}

norme1 = function(x){sum(abs(x))}

projectionL1 = function(x,U=1) {
  if (norme1(x)<=U){
    return(x)
  } else {
    return(sign(x)*simplexProj(abs(x),U))
  }
}

#' Plot an aggregation procedure
#' @describeIn oracle \code{plot}. It has one optional arguments. 
#' @param x An object of class \code{oracle}. 
#' @param sort if set to TRUE (default), it sorts the experts by performance before the plots.
#' @param col colors
#' @importFrom graphics axis box mtext par plot
#' @export 
plot.oracle <- function(x, sort = TRUE, col = NULL, ...) {
  def.par <- par(no.readonly = TRUE)
  
  if (x$d > 1) {
    x$experts <- blockToSeries(x$experts)
    x$Y <- blockToSeries(x$Y)
  }
  x$experts <- data.frame(x$experts)
  
  T <- nrow(x$experts)
  K <- ncol(x$experts)
  
  if (is.null(col)) {
    if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
      print("The RColorBrewer package must be installed to get better colors\n")
      col.palette <- 2:min((K+1),7)
    } else{
      col.palette <- RColorBrewer::brewer.pal(n = 3,name = "Set1")    
    }
    col <- col.palette
  }
  
  if (!is.null(x$names.experts)) {
    names(x$experts) <- x$names.experts
  } else {
    if (is.null(names(x$experts))) {
      names(x$experts) <- x$names.experts <- paste("X", 1:K,sep="")
    }
  }
  
  if (x$model == "expert") {
    err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
    if (sort) {
      idx.sorted <- order(c(x$loss.experts, err.unif))
      i.min <- 1
    } else {
      idx.sorted = c(K+1,1:K)
      i.min <- order(x$loss.experts)[1]
    }
    my.col <- rep(1, K + 1)
    my.col[which(idx.sorted == K + 1)] <- col[1]
    my.col[which(idx.sorted != K + 1)[i.min]] <- col[2]
    
    par(mar = c(4.5, 4, 2, 2))
    plot(c(x$loss.experts, err.unif)[idx.sorted], xlab = "", ylab = paste(x$loss.type$name, 
                                                                          "loss"), main = "Average loss suffered by the experts", axes = F, pch = 3, 
         col = my.col, lwd = 2,...)
    axis(1, at = 1:(K + 1), labels = FALSE)
    mtext(at = 1:(K + 1), text = c(names(x$experts), "Uniform")[idx.sorted], 
          side = 1, las = 2, col = my.col, line = 0.8)
    axis(2)
    box()
    
  }
  
  if (x$model == "convex" || x$model == "linear") {
    err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
    idx.sorted <- order(c(x$loss.experts, err.unif, x$loss))
    my.col <- rep(1, K + 2)
    my.col[which(idx.sorted == K + 1)] <- col[1]
    my.col[which(idx.sorted == K + 2)] <- col[2]
    my.col[which(!(idx.sorted %in% c(K + 1, K + 2)))[1]] <- col[3]
    y.max <- c(x$loss.experts, err.unif, x$loss)[idx.sorted]
    
    par(mar = c(4.5, 4, 2, 2))
    plot(c(x$loss.experts, err.unif, x$loss)[idx.sorted], xlab = "", ylab = paste(x$loss.type$name, 
                                                                                  "loss"), main = "Average loss suffered by the experts", axes = F, pch = 3, 
         col = my.col, lwd = 2)
    axis(1, at = 1:(K + 2), labels = FALSE)
    mtext(at = 1:(K + 2), text = c(names(x$experts), "Uniform", 
                                   x$model)[idx.sorted], 
          side = 1, las = 2, col = my.col, line = 0.8)
    axis(2)
    box()
  }
  
  if (x$model == "shifting") {
    L <- x$loss
    if (x$loss.type == "square") {
      L <- sqrt(x$loss)
      y.lab <- "rmse"
    } else if (x$loss.type == "absolute") {
      y.lab <- "mean absolute error"
    } else if (x$loss.type == "percentage") {
      y.lab <- "mape"
    }
    
    plot(0:(length(L) - 1), L, xlab = "Number of shifts", ylab = y.lab, type = "o", 
         pch = 20, cex = 0.6, lwd = 2, main = "Error suffered by the shifting oracle")
  }
  par(def.par)
}


#' Plot an object of class mixture
#' 
#' provides different diagnostic plots for an aggregation procedure.
#' @param x an object of class mixture. If awake is provided (i.e., some experts are unactive), 
#' their residuals and cumulative losses are computed by using the predictions of the mixture.
#' @param pause if set to TRUE (default) displays the plots separately, otherwise on a single page
#' @param col the color to use to represent each experts, if set to NULL (default) use R\code{RColorBrewer::brewer.pal(...,"Spectral"}
#' @param ... additional plotting parameters
#' 
#' 
#' @return plots representing: plot of weights of each expert in function of time, boxplots of these weights,
#' cumulative loss \eqn{L_T=\sum_{t=1}^T l_{i,t}} of each expert in function of time, cumulative residuals \eqn{\sum_{t=1}^T (y_t-f_{i,t})} of each 
#' expert's forecast in function of time, average loss suffered by the experts and the contribution of each expert to the aggregation 
#' \eqn{p_{i,t}f_{i,t}} in function of time.
#' 
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @author Yannig  Goude <yannig.goude@edf.fr>
#' @seealso See \code{\link{opera-package}} and opera-vignette for a brief example about how to use the package.
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics axis box boxplot layout legend lines matplot mtext par plot polygon text
#' @importFrom stats lowess var
#' @export 
#' 
#'
#$$plot.mixture <- function(x, pause = FALSE, col = NULL, ...) {
plot.mixture <- function(x, pause = FALSE, col = NULL,...) {

  #$$def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if (pause) par(ask=TRUE)
  x$experts <- data.frame(x$experts)
  K <- length(x$experts)
  w.order <- order(apply(x$weights,2,mean),decreasing = TRUE)
  
  if (is.null(col)) {
    if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
      print("The RColorBrewer package must be installed to get better colors\n")
      col <- 2:min((K+1),7)
    } else{
      col <- rev(RColorBrewer::brewer.pal(n = max(min(K,11),4),name = "Spectral"))[1:min(K,11)]
    }
  }

  #$$
  if (K<=11) col <- c("springgreen2","springgreen4","gold2","goldenrod4","hotpink","deeppink","blue","deepskyblue","orange","brown","red")
  isHere=rep(TRUE,length(liste_exp_with))
  for (iexp in 1:length(liste_exp_with)) {
    if(length(which(names(x$experts)==liste_exp_with[iexp]))==0){
      isHere[iexp]=FALSE
    }
  }
  col=col[isHere]
  nonOrderedCol=col
  #$$

  my.colors <- col
  
  col <- numeric(K)
  if (K <= length(my.colors)) {
    col[w.order] <- my.colors[1:K]
  } else {
    col[w.order] <- c(my.colors, rep(my.colors[length(my.colors)],K-length(my.colors)))
  }
  #$$if (!pause) {
    #$$layout(matrix(c(1,2,3,4,5,6),nrow = 3,ncol =  2, byrow = TRUE))  
  #$$}
  #$$par(mar = c(3, 3, 1.6, 0.1), mgp = c(2, 0.5, 0))
  T <- x$T
  d <- x$d
  x$experts <- data.frame(x$experts)
  names(x$experts)[which(names(x$experts)=="mos.aro")]="ppm.aro"
  names(x$experts)[which(names(x$experts)=="mos.arp")]="ppm.arp"
  names(x$experts)[which(names(x$experts)=="mos.cep")]="ppm.ifs"
  x$names.experts[which(x$names.experts=="mos.aro")]="ppm.aro"
  x$names.experts[which(x$names.experts=="mos.arp")]="ppm.arp"
  x$names.experts[which(x$names.experts=="mos.cep")]="ppm.ifs"
  x$Y <- c(t(x$Y))
  x$prediction <- c(t(x$prediction))
  x$weights <- data.frame(x$weights)

  if (!is.null(x$awake_transition)) x$awake <- x$awake_transition#$$
  normalized.weights <-x$weights* x$awake #$$
  normalized.weights <- normalized.weights /rowSums(normalized.weights ) #$$
  x$weights <- normalized.weights #$$
  
  if (!is.null(x$names.experts)) {
    names(x$weights) <- names(x$experts) <- x$names.experts
  } else {
    if (is.null(names(x$experts))) {
      names(x$weights) <- names(x$experts) <- x$names.experts <- paste("X", 1:K,sep="")
    }
  }

  l.names <- max(nchar(names(x$experts))) / 3 + 1.7

  #plot weights
    weights=x$weights
    names(weights) <- names(x$experts)
    iteration=1:T
    weights=cbind(iteration,weights)
    data=melt(weights,id='iteration')#on passe au format long
    colnames(data)=c("iteration","expert","weight")
    #ypos=(1-cumsum(data[which(data$iteration==T),]$weight))+data[which(data$iteration==T),"weight"]/2
    # xpos2=c()
    # ypos2=c()
    # weightsBis=weights[,-which(colnames(weights)=="iteration")]
    # for (iexp in 1:length(names.experts)){
    #   newx=which(weightsBis[,iexp]==max(weightsBis[,iexp]))
    #   newx=newx[1]
    #   xpos2=c(xpos2,newx)
    #   #newy=(1-cumsum(weightsBis[newx,])[iexp])+weightsBis[newx,iexp]/2
    #   newy=(1-cumsum(data[which(data$iteration==newx),]$weight))[iexp]+data[which(data$iteration==newx & data$expert==names.experts[iexp]),]$weight/2
    #   ypos2=c(ypos2,newy)
    # }
    # ypos2=unlist(ypos2)
    # label2id=which(ypos2>=0.1)
    # xpos2=xpos2[label2id]
    # ypos2=ypos2[label2id]
    ggplot(data, aes(x=iteration, y=weight, fill=expert)) +
      theme_classic() + #to remove grey background and grid
      geom_area()+
      scale_fill_manual(values=nonOrderedCol)+ #colours
      #geom_text(data=data[which(data$iteration==T),],aes(y=ypos,label=expert),x=T,hjust=0,color=nonOrderedCol,size=5)+
      theme(axis.text.x = element_text(size=22,color = "black"), axis.text.y = element_text(size=22,color = "black"),
            axis.title = element_text(size = 22,color = "black",face = "bold"),
            legend.text=element_text(size=22),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-25,-25),
            legend.title=element_text(size=22,,face = "bold"))+
      geom_vline(xintercept=624, linetype='dashed', color='grey20',size=1.5)+#start of the event in chamonix
      geom_vline(xintercept=634, linetype='dashed', color='grey20',size=1.5)#end of the event in chamonix
      
      #geom_text(data=data[which(data$iteration==T),],x=xpos2,y=ypos2,label=names.experts,color="black")
    ggsave(file=paste(dir_plot,dir_mixture,"weights",".pdf",sep=""), width=10, height=8, dpi=600)

  # pdf(paste(dir_plot,dir_mixture,"weights",".pdf",sep="")) #$$ pdf(paste(dir_plot,dir_mixture,"weights",".pdf",sep=""), width = 11, height = 8) #$$
  # if (x$model=="Boa") agregName="BOA"
  # else agregName= x$model
  # if (x$model == "Ridge") {
  #   # Linear aggregation rule
  #   #$$par(mar = c(3, 3, 2, l.names/2), mgp = c(1, 0.5, 0))
  #   matplot(x$weights, type = "l", xlab = "", ylab = "", lty = 1:5, main = agregName, col = col,...)
  #   mtext(side = 2, text = "Weights", line = 1.8, cex = 1.1)
  #   mtext(side = 1, text = "Ieration", line = 1.8, cex = 1.1)
  #   # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
  #   mtext(side = 4, text = names(x$experts), at = x$weights[T,], las = 2, col = col, cex= 0.7, line = 0.1)
  # } else {
  #   # Convex aggregation rule
  #   #$$par(mar = c(3, 3, 2, l.names/2), mgp = c(1, 0.5, 0))
  #   plot(c(1), type = "l", col = 1:8, lwd = 2, axes = F, xlim = c(1, T),
  #       ylim = c(0,1), ylab = "", xlab = "")#expression("BOA"^"s"))
  #   mtext(side = 2, text = "Weights", line = 1.8, cex = 1.1)
  #   mtext(side = 1, text = "Iteration", line = 1.8, cex = 1.1)
  #   # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
  #   x.idx <- c(1, 1:T, T:1)
  #   w.summed <- rep(1,T)
  #   w.summed[which(is.na(x$Y))]=0 #$$ iterations where the observations are missing
  #   i.remaining = rep(TRUE,K)
  #   i.order <- rep(0,K)
  #   for (i in 1:K) {
  #     if (i <K){
  #       j <- which(i.remaining)[which.min(apply(x$weights[,i.remaining],2,function(p){sqrt(var(w.summed-p,na.rm=T))}))]
  #     } else {
  #       j <- which(i.remaining)
  #     }

  #     i.order[i] <- j
  #     y.idx <- c(0, w.summed, rep(0, T))
  #     #polygon(x = x.idx, y = y.idx, col = col[j], border=NA)
  #     polygon(x = x.idx, y = y.idx, col = my.colors[j], border=NA) #to have the same colors for the models on the different graphs
  #     w.summed.old <- w.summed
  #     w.summed <- w.summed - x$weights[,j]
  #     w.summed[which(is.na(x$Y))]=0 #$$ iterations where the observations are missing

  #     i.remaining[j] <- FALSE
  #     writeLegend(f = w.summed.old,w.summed,name = names(x$experts)[j],cex=1)
  #   }
  #   names.toWrite <- names(x$experts)
  #   names.toWrite[w.order[-(1:min(K,15))]] <- ""
  #   mtext(side = 4, text = names.toWrite[i.order], 
  #         #at = (1-cumsum(c(x$weights[T,i.order])))  + x$weights[T,i.order]/2, las = 2, col = col[i.order], cex= 0.5, line = 0.3)
  #         at = (1-cumsum(c(x$weights[T,i.order])))  + x$weights[T,i.order]/2, las = 2, col = my.colors[i.order], cex= 0.95, line = -0.8) #to have the same colors for the models on the different graphs
  #   axis(1)
  #   axis(2)
  #   #box() #to draw a box around the weights
  # }
  # dev.off() #$$

  # Box plot
  pdf(paste(dir_plot,dir_mixture,"boxplots",".pdf",sep="")) #$$
  if (!is.null(x$awake)) {
    pond <- apply(x$awake[d*(1:T),],1,sum)
    normalized.weights <- x$weights * pond / (K*x$awake[d*(1:T),])
    normalized.weights[x$awake[d*(1:T),] == pond] <- NaN
    normalized.weights[normalized.weights>exp(700)] <-NaN #to deal with the case where the expert is asleep (is.infinite() doesn't work here !?)
  } else {
    normalized.weights <- x$weights 
  }
  
  i.order <- w.order[1:min(K,20)]
  par(mar = c(l.names, 3, 1.6, 0.1))
  #boxplot(normalized.weights[,i.order], main = "Weights associated with the experts", col = col[i.order], axes = FALSE, pch='.')
  boxplot(normalized.weights[,i.order], main = "Weights associated with the experts", col = my.colors[i.order], axes = FALSE, pch='.') #to have the same colors for the models on the different graphs
  mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
  axis(1, at = 1:(min(K,20)), labels = FALSE)
  #mtext(at = 1:min(K,20), text = names(x$weights)[i.order], side = 1, las = 2, col = col[i.order], line = 0.8)
  mtext(at = 1:min(K,20), text = names(x$weights)[i.order], side = 1, las = 2, col = my.colors[i.order], line = 0.8) #to have the same colors for the models on the different graphs
  axis(2)
  box()
  dev.off() #$$
  
  #note: always pass alpha on the 0-255 scale
  makeTransparent<-function(someColor, alpha=220)
  {
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
  }
  
  # Cumulative loss
  pdf(paste(dir_plot,dir_mixture,"squareloss",".pdf",sep="")) #$$
  par(mar = c(1.5, 3, 2.5, l.names/2), mgp = c(1, 0.5, 0))
  pred.experts <- (x$experts * x$awake + x$prediction * (1-x$awake))

  x$Y[which(is.na(x$Y))]=0 #to deal with missing observations
  
  cumul.losses <- apply(loss(pred.experts, x$Y, x$loss.type), 2, cumsum)[seq(d,T*d,by=d),]
  cumul.exploss <- cumsum(loss(x$prediction, x$Y, x$loss.type))[seq(d,T*d,by=d)]

  matplot(cumul.losses, type = "l", lty = 1, xlab = "", ylab = "", 
          #main = paste("Cumulative", x$loss.type$name, "loss"), col = makeTransparent(col), ylim = range(c(cumul.losses,cumul.exploss)))
          main = paste("Cumulative", x$loss.type$name, "loss"), col = makeTransparent(my.colors), ylim = range(c(cumul.losses,cumul.exploss))) #to have the same colors for the models on the different graphs
  lines(cumul.exploss, col = 1, lwd = 2)
  mtext(side = 2, text = "Cumulative loss", line = 1.8, cex = 1)
  # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
  #mtext(side = 4, text = x$names.experts, at = cumul.losses[T,], las = 2, col = makeTransparent(col), cex= 0.5, line = 0.3)
  mtext(side = 4, text = x$names.experts, at = cumul.losses[T,], las = 2, col = makeTransparent(my.colors), cex= 0.5, line = 0.3) #to have the same colors for the models on the different graphs
  legend("topleft", c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
  dev.off() #$$
  
  # Cumulative residuals
  pdf(paste(dir_plot,dir_mixture,"residuals",".pdf",sep="")) #$$
  par(mar = c(1.5, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
  # cumul.residuals <- apply(x$Y - pred.experts, 2, cumsum)[seq(d,T*d,by=d),]
  # cumul.expres <- cumsum(x$Y - x$prediction)[seq(d,T*d,by=d)]
  #on prend pas obs-pred comme dans opera mais pred-obs pour avoir un biais positif lorsque le modele est trop chaud
  cumul.residuals <- apply(pred.experts - x$Y, 2, cumsum)[seq(d,T*d,by=d),]
  cumul.expres <- cumsum(x$prediction - x$Y)[seq(d,T*d,by=d)]
  matplot(cumul.residuals, type = "l", lty = 1, xlab = "", ylab = "", 
          #main = paste("Cumulative residuals"), col = makeTransparent(col), ylim = range(c(cumul.residuals,cumul.expres)))
          main = paste("Cumulative residuals"), col = makeTransparent(my.colors), ylim = range(c(cumul.residuals,cumul.expres))) #to have the same colors for the models on the different graphs
  abline(h=0, col="grey", lty="longdash")
  lines(cumul.expres, col = 1, lwd = 2)
  mtext(side = 2, text = "Cumulative residuals", line = 1.8, cex = 1)
  # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
  if (max(cumul.residuals) > abs(min(cumul.residuals))) {
    place = "topleft"
  } else {
    place = "bottomleft"
  }
  #mtext(side = 4, text = x$names.experts, at = cumul.residuals[T,], las = 2, col = col, cex= 0.5, line = 0.3)
  mtext(side = 4, text = x$names.experts, at = cumul.residuals[T,], las = 2, col = my.colors, cex= 0.5, line = 0.3) #to have the same colors for the models on the different graphs
  legend(place, c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
  dev.off() #$$

  #losses
  pdf(paste(dir_plot,dir_mixture,"averageloss",".pdf",sep="")) #$$
  l.names <- max(max(nchar(names(x$experts))) / 3 + 1.7,4)
  x$loss.experts <- apply(loss(x = pred.experts,y = x$Y,loss.type = x$loss.type),2,mean)
  err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
  err.mixt <- x$loss
  idx.sorted <- order(c(x$loss.experts, err.unif, err.mixt))
  #my.col <- c(col,1,1)[idx.sorted]
  my.col <- c(my.colors,1,1)[idx.sorted] #to have the same colors for the models on the different graphs
  my.pch <- c(rep(20, K),8,8)[idx.sorted]
  
  par(mar = c(l.names, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
  plot(c(x$loss.experts, err.unif, err.mixt)[idx.sorted], xlab = "", ylab = "", main = "Average loss", axes = F, pch = my.pch, 
       col = my.col, lwd = 2,type='b')
  mtext(side = 2, text = paste(x$loss.type$name,"loss"), line = 1.8, cex = 2)
  axis(1, at = 1:(K + 2), labels = FALSE)
  mtext(at = 1:(K + 2), text = c(names(x$experts), "Uniform", x$model)[idx.sorted], 
        side = 1, las = 2, col = my.col, line = 0.8,cex = 1.1)
  axis(2)
  box()
  dev.off() #$$

  # cumulative plot of the series
  pdf(paste(dir_plot,dir_mixture,"contribution",".pdf",sep="")) #$$
  par(mar = c(2, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
  if (x$d ==1) {
    #cumulativePlot(W = x$weights,X = x$experts, Y = x$Y,smooth = TRUE,alpha = 0.01,plot.Y = TRUE, col.pal = rev(my.colors))
    cumulativePlot(W = x$weights,X = x$experts, Y = x$Y,smooth = TRUE,alpha = 0.01,plot.Y = TRUE, col.pal = my.colors[w.order]) #to have the same colors for the models on the different graphs
  } else {
    X <- apply(seriesToBlock(X = x$experts,d = x$d),c(1,3),mean)
    Y <- apply(seriesToBlock(x$Y,d = x$d),1,mean)
    colnames(X) <- x$names.experts
    cumulativePlot(W = x$weights,X = X, Y = Y,smooth = TRUE,alpha = 0.01,plot.Y = TRUE, col.pal = col)    
  }
  dev.off() #$$
  #$$par(def.par)
} 



writeLegend <- function(f,g,name,Y.lim=c(0,1),cex=0.8, ...) {
  tau = Y.lim[2]/20
  Tab = matrix(0,ncol = 2, nrow = 100)
  y.seq <- seq(Y.lim[1],Y.lim[2],length.out = 100)
  for (i in 1:100) {
    x = y.seq[i]
    sel = which(g < x & f > x + tau)
    temp <- cumsum(c(1, diff(sel) - 1))
    temp2 <- rle(temp)
    Tab[i,1] <- max(temp2$lengths)
    Tab[i,2] <- sel[which(temp == with(temp2, values[which.max(lengths)]))][1]
  }
  id = which.max(Tab[,1])
  x <- y.seq[id]
  l <- Tab[id,1]
  v <- Tab[id,2]
  if (l > length(f)/20){
    j = floor(60 *l/length(f))
    text(v+l/2,x+tau/2,substr(name,1,j),cex = cex,...)
  }
}

cumulativePlot<-function(W,X,Y,col.pal=NULL, smooth = FALSE, plot.Y = FALSE, alpha = 0.1)
{
  time<-c(1:nrow(X))
  #$$active.experts<-which(colMeans(W)>0)
  active.experts<-which(colMeans(W,na.rm=TRUE)>0)#$$
  W<-W[,active.experts]  
  X<-X[,active.experts]
  no.observations<-which(is.nan(rowSums(W))) #$$ to deal with missing observations
  W[no.observations,]<-1/ncol(W) #$$ to deal with missing observations
  X[no.observations,]<-0 #$$ to deal with missing observations
  
  K <- ncol(X)
  
  
  if(is.null(col.pal)) col.pal <- RColorBrewer::brewer.pal(n = min(K,9),name = "Spectral")
  if (length(col.pal) < K) col.pal <- c(rep(col.pal[1],K-length(col.pal)),col.pal)
  
  
  o<-order(colSums(W),decreasing = F)
  mat<-W[,o]*X[,o]
  Agg<-apply(mat,1,sum)
  colnames(mat)<-colnames(X)[o]
  
  if (!smooth)Y.lim = range(Agg,Y,mat)
  if (smooth) 
  {
    y.lo<-lowess(x = time,y = Y,f = alpha)$y
    Agg.lo<-lowess(x = time,y = Agg,f = alpha)$y
    
    mat.lo<-apply(mat,2,function(z){lowess(x = time,y = z,f = alpha)$y})
    Y.lim = range(Agg.lo,mat.lo)
  }
  
  plot(x = NULL,y = NULL,col=col.pal[1], type='l', xaxt='n',ylim=Y.lim,lty='dotted',
       yaxt='n',xlab="",ylab="",lwd=3,xlim = range(time),
       main = paste("Contribution of each expert to prediction"))
  y.summed <- Agg
  for(i in rev(c(1:ncol(mat))))
  {
    if (!smooth) addPoly(time,y.summed,col=col.pal[i])
    if (smooth) addPoly(time,lowess(y.summed,f = alpha)$y,col=col.pal[i])
    y.summed.old <- y.summed
    y.summed <- y.summed - mat[,i]
    if (!smooth) writeLegend(f=y.summed.old,g= y.summed, name = colnames(mat)[i],Y.lim,col='black')
    if (smooth) writeLegend(f=lowess(y.summed.old,f=alpha/10)$y,g=lowess(y.summed,f=alpha/10)$y, name = colnames(mat)[i],Y.lim,col='black')
  }
  if (plot.Y && !smooth) lines(time,Y,col=1,lwd=2,lty='dotted')
  if (plot.Y && smooth) lines(lowess(x = time,y = Y,f = alpha)$y,col=1,lwd=2,lty='dotted')
  axis(1)
  axis(2)
}


addPoly<-function(x,y,col)
{
  xx <- c(x, rev(x))
  yy <- c(rep(0, length(x)), rev(y))
  polygon(xx, yy, col=col, border=NA)
}


summary.mixture <- function(object, ...) {
  if (is.null(object$Y)) {
    K <- "Unknown"
    T <- 0
    d <- "Unknown"
    TAB <- c("No losses yet")
  } else {
    T <- object$T
    K <- length(object$coefficients)
    d <- object$d
    
    rmse.algo <- sqrt(mean(loss(c(object$prediction), c(object$Y), loss.type = "square")))
    mape.algo <- mean(loss(c(object$prediction), c(object$Y), loss.type = "percentage"))
    rmse.unif <- sqrt(lossConv(rep(1/K, K), c(t(object$Y)), object$experts, awake = object$awake))
    mape.unif <- lossConv(rep(1/K, K), c(t(object$Y)), object$experts, awake = object$awake, 
                          loss.type = "percentage")
    
    TAB <- data.frame(rmse = c(rmse.algo, rmse.unif), mape = c(mape.algo, mape.unif))
    rownames(TAB) <- c(object$model, "Uniform")
  }
  
  res <- list(object = object, coefficients = object$coefficients, losses = TAB, 
              n.experts = K, n.observations = T, n.dimension = d)
  class(res) <- "summary.mixture"
  res
}


print.summary.mixture <- function(x, ...) {
  print(x$object)
  cat("\nNumber of experts: ", x$n.experts)
  cat("\nNumber of observations: ", x$n.observations)
  cat("\nDimension of the data: ", x$n.dimension, "\n\n")
  
  if (!is.null(dim(x$losses))) {
    print(signif(x$losses, digits = 3))
  }
}

# truncate1 a real number The function \code{truncate1} is used to avoid
# instability in \code{R} due to very smal or very big numbers.  It truncate1s
# numbers so that they lie between \code{exp(-700)} and \code{exp(700)}.
# @param x A number to be truncate1d @return The truncated value of \code{x}
# @author Pierre Gaillard <pierre@@gaillard.me> @keywords ~kwd1 ~kwd2
truncate1 <- function(x) {
  #cf commit git hub de pierre gaillard : Update function truncate1 to make it faster.
  #pmin and pmax are run only when necessary (it is significantly longer than min and max to execute).
  #the following code should be faster than just pmin(pmax(x, exp(-700)), exp(700))
  is_sup <- max(x) > exp(700)
  is_inf <- min(x) <= exp(-700)

  res <- x

  if (is_sup) {
    res <- pmin(res, exp(700))
  }
  if (is_inf) {
    res <- pmax(res, exp(-700))
  }

  return(res)
}
