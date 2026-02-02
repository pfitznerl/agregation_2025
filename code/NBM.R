nbm <- function(y, experts, eta=0.025, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                w0 = NULL, training = NULL) {
  experts <- as.matrix(experts)
  bExperts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  bExperts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    y[which(is.na(y))]<- 0
    awake[idx.na] <- 1 #to avoid 0 division
    experts[idx.na] <- 0
    bExperts[idx.na] <- 0
  }
  
  R <- rep(0, N)  # Regret vector
  B <- rep(0, N)  # Bias vector
  MSE <- rep(exp(-350), N)  # MSE vector
  pred <- rep(0, T)  # Prediction vector
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture

  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    weights[t,] <- rep(1, N)
  } else weights[1,] <- as.matrix(w0,ncol=N,nrow=1)

  #if (!is.null(training)) {
  if (!is.null(training$R)) {
    w0 <- training$w0
    R <- training$R
    eta<-training$eta
    MSE <- training$MSE
    B <- training$B
  }
  iteration<-training$iteration
  date_valid<-training$date_valid

  for (t in 1:T) {
    # Weight update
    weights[t,]=1/MSE #weights
    weights[t,]=c(weights[t,]/sum(weights,na.rm=TRUE)) #normalized weights
    # Prediction and losses
    bExperts[t, ]=experts[t, ]-B #debiased predictions
    pred[t]=weights[t, ]%*%bExperts[t, ] #NBM prediction
    lpred <- lossPred(pred[t], y[t], pred[t], loss.type, loss.gradient)
    lexp <- lossPred(experts[t, ], y[t], pred[t], loss.type, loss.gradient)

    # Regret update
    r <- c(c(lpred) - lexp)
    R <- R + awake[t, ] * r

    if (as.numeric(strsplit(date_valid, "-")[[1]][2]) >=6 & as.numeric(strsplit(date_valid, "-")[[1]][2]) <= 8) {
        eta=0.05
    } else eta=0.025
    b=experts[t,] - y[t]
    mse=(experts[t,] - y[t])^2

    if(iteration>1){
      for (i in 1:N){
        NROW=nrow(training$b)
        B[i]=sum((1-eta)^(NROW:1)*eta*training$b[,i])+eta*b[i]
        MSE[i]=sum((1-eta)^(NROW:1)*eta*training$mse[,i])+eta*mse[i]
      }
    } else {
      B=eta*b
      MSE=eta*mse
    }

    MSE[which(MSE==0)]=exp(-700)


  }
  w <- weights[t,]
  
  object <- list(model = "NBM",coefficients = w, loss.type = loss.type)
  
  R <- R
  object$parameters <- list(eta = eta)
  object$weights <- weights
  object$prediction <- pred
  object$training <- list(R = R, w0 = w0,
                          eta=eta,b=rbind(training$b,b),
                          mse=rbind(training$mse,mse),MSE=MSE,B=B)
  
  return(object)
}