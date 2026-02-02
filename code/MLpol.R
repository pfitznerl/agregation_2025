#MLpol
MLpol <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
  training = NULL) {

  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    y[which(is.na(y))]<- 0
    awake[idx.na] <- 1 #to avoid 0 division
    experts[idx.na] <- 0
  }
  
  # weight assigned to each expert
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  
  # Initialization or the learning parameter
  B <- 0
  eta <- matrix(exp(700), ncol = N, nrow = T + 1)
  # regret suffered by each expert
  R <- rep(0, N)
  
  #if (!is.null(training)) {
  if (!is.null(training$R)) {
    eta[1, ] <- training$eta
    R <- training$R
    B <- training$B
  } else {
    training$eta <- list(eta = eta[1, ])
  }
  sta<-training$sta
  ech<-training$ech
  window<-training$window
  iteration<-training$iteration
  
  w <- rep(0, N)

  for (t in 1:T) {
    if(iteration==1)print("MLpol")

    # We check if there is at least one expert with positive weight
    if (max(awake[t, ] * R) > 0) {
      w <- eta[t, ] * pmax(R, 0)/sum(eta[t, ] * pmax(R, 0))
    } else {
      w <- rep(1, N)
    }
    
    # form the mixture and the prediction
    p <- awake[t, ] * w/sum(awake[t, ] * w)
    pred <- experts[t, ] %*% p
    
    # save the mixture and the prediction
    weights[t, ] <- p
    prediction[t] <- pred
    
    # Observe losses
    lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    
    # Update the regret and the weight
    r <- awake[t, ] * (lpred - lexp)
    R <- R + r
    
    # Update the learning rate
    newB <- max(B, max(r^2))
    eta[t + 1, ] <- 1/(1/eta[t, ] + r^2 + newB - B)
    B <- newB
  }

  # We check if there is at least one expert with positive weight
  if (max(R) > 0 & !is.na(y)) {
    w <- eta[T + 1, ] * pmax(R, 0)/sum(eta[T + 1, ] * pmax(R, 0))
  } else {
    w <- rep(1/N, N)
  }
  object <- list(model = "MLpol", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, B = B,r=rbind(training$r,r))
  class(object) <- "mixture"
  return(object)
}
