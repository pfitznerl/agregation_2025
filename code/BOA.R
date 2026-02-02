#soco BOA
BOA <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                w0 = NULL, training = NULL) {
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)

  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }

  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    awake[which(is.na(y)),] <- 1 #to avoid 0 division
    experts[which(is.na(y)),] <- 0
    y[which(is.na(y))]<- 0 #has to be done after awake[which(is.na(y)),] <- 1 and experts[which(is.na(y)),] <- 0
  }
  
  R <- rep(0, N)
  R.reg <- R
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  w <- w0
  eta <- matrix(exp(350), ncol = N, nrow = T + 1)
  V <- 0

  if (!is.null(training$R)) {
    w0 <- training$w0
    R <- training$R
    R.reg <- training$R.reg
    w <- truncate1(training$eta*exp(log(w0) + training$eta * R.reg))
    eta[1, ] <- training$eta
    V <- training$V
  }

  for (t in 1:T) {
    #to make the uniform aggregation
    #w=rep(1/N,N)

    p <- awake[t, ] * w/sum(awake[t, ] * w)
    pred <- experts[t, ] %*% p
    
    
    weights[t, ] <- p
    prediction[t] <- pred

    #loss of the prediction
    lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    #vector, losses for each of the M experts
    lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)

    #to make the best compound expert oracle
    # idx=which(abs(experts-y)==min(abs(experts-y),na.rm=T))
    # prediction[t] <- experts[idx]

    # Vector, instantaneous regret
    r <-  awake[t, ] * (lpred - lexp) # -l_{j,t} in wintenberger

    # Update the learning rates
    V <- V + 2.2*r^2 #vector, cumulated quadratic loss for each experts
    if(length(which(V!=0))==N) {
      eta[t + 1, ] <- pmin(sqrt(1/V),exp(350)) 
    } else {
      eta[t + 1, which(V!=0)]= pmin(sqrt(1/V),exp(350)) [which(V>0)]
      eta[t + 1, which(V==0)]= 0
    }

    r.reg <- r - eta[t+1, ] * r^2

    # Update the regret and the regularized regret used by BOA
    R <- R + r
    R.reg <- R.reg + r.reg
    
    w <- truncate1(eta[t + 1, ]*exp(log(w0) + eta[t + 1, ] * R.reg))
    
  }
  
  object <- list(model = "BOA", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w/sum(w))
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction

  object$training <- list(eta = eta[T + 1, ], R = R, w0 = w0, R.reg = R.reg, V= V, r.reg=rbind(training$r.reg,r.reg),r=rbind(training$r,r))
  class(object) <- "mixture"

  return(object)
}