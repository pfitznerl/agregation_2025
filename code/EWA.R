#EWA
ewa <- function(y, experts, eta, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                w0 = NULL, training = NULL) {
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
  
  lossMax=exp(-700)
  V=exp(-700)
  R <- rep(0, N)  # Regret vector
  pred <- rep(0, T)  # Prediction vector
  cumulativeLoss <- rep(0,T+1)  # Cumulative loss of the mixture (needed for EWACALIB), T+1 because we need t=0
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  iLoss <- rep(0, T)

  if (!is.null(training$R)) {
    w0 <- training$w0
    R <- training$R
    cumulativeLoss[1] <- training$cumulativeLoss[length(training$cumulativeLoss)]
    lossMax=training$lossMax
    eta<-training$eta
    V<-training$V
  }
  iteration<-training$iteration

  for (t in 1:T) {
    # Weight update
    weights[t, ] <- awake[t, ]*truncate1(exp(log(w0) + eta * R))
    weights[t, ] <- weights[t, ]/sum(weights[t, ])
    
    # Prediction and losses
    pred[t] <- experts[t, ] %*% weights[t, ]
    iLoss[t] <- loss(pred[t], y[t], loss.type)
    cumulativeLoss[t+1] <- cumulativeLoss[t] + iLoss[t]
    lpred <- lossPred(pred[t], y[t], pred[t], loss.type, loss.gradient)
    lexp <- lossPred(experts[t, ], y[t], pred[t], loss.type, loss.gradient)

    # Regret update
    r <- c(c(lpred) - lexp)
    R <- R + awake[t, ] * r
    v=sum(weights[t,]*(awake[t, ]*lexp^2))- (sum(weights[t,]*(awake[t, ]*lexp)))^2
    v=max(v,exp(-700))
    V=V+v

    lossMax=max(lossMax,abs(max(lexp,na.rm=T)-min(lexp,na.rm=T)),na.rm=T)
    eta=min(1/lossMax,1.0739*sqrt(log(N)/V)) #Bianchi 2007 theorem 6
  }
  w <- truncate1(exp(log(w0) + eta * R))/sum(truncate1(exp(log(w0)+eta * R)))
  
  object <- list(model = "EWA", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w)
  
  R <- R
  object$parameters <- list(eta = eta)
  object$weights <- weights
  object$prediction <- pred
  
  object$training <- list(R = R, w0 = w0, r=rbind(training$r,r), iLoss=c(training$iLoss,iLoss), lpred=c(training$lpred,lpred),
                          cumulativeLoss = c(training$cumulativeLoss,cumulativeLoss[-1]),eta=eta,lossMax=lossMax,V=V, v=c(training$v,v))
  
  return(object)
}