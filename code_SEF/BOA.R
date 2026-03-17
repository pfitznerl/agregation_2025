#SOCO version of BOA with sleeping experts awaken by transition object
BOA <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                w0 = NULL, training = NULL) {
  experts <- as.matrix(experts)
  names.experts=colnames(experts)
  experts_for_xgboost=experts #to train xgboost with real data (without putting missing experts or observations equal to zero)

  N <- ncol(experts)
  T <- nrow(experts)

  eV=eValue(N=N)
  threshold=0.5 # to wake up the experts
  xgboost_threshold=0.5
  transObject_threshold=probTransition(names.experts=names.experts,predType="predErr",nTrain=365*3*18*33*5,
                                        threshold=xgboost_threshold) #probTransition object to predict if the agregation will have a loss bigger than a threshold

  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
      w0 <- rep(1, N)
  }
  #dataToSave=data.frame()
  awake <- as.matrix(awake)
  is_na_spe=is.na(as.data.frame(experts) %>% select(contains(c("Q10","Q30","Q70","Q90")))) #matix true false if spe experts are na
    idx.na <- which(is.na(experts))
    awake[idx.na] <- 0
    experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    id_na_y=which(is.na(y))
    awake[id_na_y,] <- 1 #to avoid 0 division
    experts[id_na_y,] <- 0
    y[id_na_y]<- 0 #has to be done after id_na_y=which(is.na(y))
    }

  colnames(awake)=names.experts

  liste_kf=mylist <- sapply(paste("kf_",names.experts,sep=""),function(x) NULL) #list, will contain the kalman filters for each expert
  X=mylist <- sapply(paste("X_",names.experts,sep=""),function(x) NULL) #list, will contain the predictors for the kalman filter
  Y=c() #list, will contain the observations for the kalman filter

  R <- rep(0, N)
  R.reg <- R
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  w <- w0
  V <- rep(0,N)
  eta <- matrix(0, ncol = N, nrow = T + 1)
  ref_experts <- as.data.frame(experts) %>% select(-contains(c("Q10","Q30","Q70","Q90")))
  ref_experts<-as.matrix(ref_experts)
  ref_N=ncol(ref_experts)
  ref_awake <- as.data.frame(awake) %>% select(-contains(c("Q10","Q30","Q70","Q90")))
  ref_awake<-as.matrix(ref_awake)
  ref_w0=rep(1,ref_N)
  ref_w=ref_w0
    ref_eta <- matrix(0, ncol = ref_N, nrow = T + 1)
  ref_V <- rep(0,ref_N)
  ref_R <- rep(0, ref_N)
  ref_R.reg <- ref_R
  shapleyValues=data.frame() 
  NERR=0#$$ count of the large mistake of GBRT
  A1=0
  B1=0
  C1=0
  D1=0
  A2=0
  B2=0
  C2=0
  D2=0
  a1=-1
  b1=1
  a2=-1
  b2=1
  asleep1=0.5
  asleep2=0.5
  score_new_vs_old1=NULL #we only add a score if the sleeping experts are woken up
  score_new_vs_old2=NULL #we only add a score if the sleeping experts are woken up

  lnSpread<-c()#will contain the log of the expert's variance
  if (!is.null(training$R)) {
      w0 <- training$w0
      R <- training$R
      R.reg <- training$R.reg
      w <- truncate1(exp(log(training$eta)+log(w0) + training$eta * R.reg))
      eta[1, ] <- training$eta
      V <- training$V
      ref_w0 <- training$ref_w0
      ref_R <- training$ref_R
      ref_R.reg <- training$ref_R.reg
      ref_w <- truncate1(exp(log(training$ref_eta)+log(ref_w0) + training$ref_eta * ref_R.reg))
      ref_eta[1, ] <- training$ref_eta
      ref_V <- training$ref_V
      eV <- training$eV
      transObject_threshold <- training$transObject_threshold
      liste_kf <- training$liste_kf
      last_l <- training$last_l
      X <- training$X #for kalman filter
      Y <- training$Y #for kalman filter
      lnSpread <- training$lnSpread
      A1 <- training$A1
      B1 <- training$B1
      C1 <- training$C1
      D1 <- training$D1
      A2 <- training$A2
      B2 <- training$B2
      C2 <- training$C2
      D2 <- training$D2
      shapleyValues <- training$shapleyValues
      asleep1<-mean(training$score_new_vs_old1)
      if(is.na(asleep1))asleep1=0.5
      asleep2<-mean(training$score_new_vs_old2)
      if(is.na(asleep2))asleep2=0.5
      a1<-training$a1
      b1<-training$b1
      a2<-training$a2
      b2<-training$b2
      NERR<-training$NERR
      #dataToSave=training$dataToSave
  }
  date.valid<-training$date.valid#$$
  sta<-training$sta#$$
  ech<-training$ech#$$
  window<-training$window#$$
  iteration<-training$iteration#$$

  for (t in 1:T) {
  if(iteration==1) print("sleeping BOA")

  fortinSDQ_T=0 #variance of the ensemble forecast pearp, like in fortin 2014
  for (iexp in 1:N){
      fortinSDQ_T=fortinSDQ_T+ (mean(experts)-experts[iexp])^2
      }
  fortinSDQ_T=fortinSDQ_T/(N-1)
  fortinSDQ_T=sqrt(fortinSDQ_T*(N+1)/N)
  lnSpread <- c(lnSpread,log(fortinSDQ_T))
  beta <- sqrt(var(lnSpread))-sqrt(var(lnSpread[-c(length(lnSpread))])) #cf The Relationship between Ensemble Spread and Ensemble Mean Skill whitaker et loughe 1998

  runVar=computeRunVar(sta,res,iteration,names.experts) #Variance of the variance in the run, and mean variance

  liste_kf_pred_mean=mylist <- rep(NA,N) #vector, will contain the kalman filters prediction's for each expert
  liste_kf_pred_sd=mylist <- rep(NA,N) #vector, will contain the kalman filters standard deviations predictions for each expert
      for (iexp in 1:N) {
          if (!is.null(liste_kf[[iexp]])) {
      newX=matrix(c(1,experts[t,iexp],last_l[iexp]),1,3)
      liste_kf[[iexp]]=predict(liste_kf[[iexp]],newX,online=FALSE,type="model") #prediction of the kalman filter based on newX, without y[t]
      liste_kf_pred_mean[iexp]=liste_kf[[iexp]]$pred_mean
      liste_kf_pred_sd[iexp]=liste_kf[[iexp]]$pred_sd
    } else {
      newX=matrix(c(1,experts[t,iexp],0),1,3)
      d=ncol(newX)
      theta=matrix(1,d,1)
      P=diag(d)
      sig=1
      liste_kf_pred_mean[iexp]=t(theta)%*%t(newX)
      liste_kf_pred_sd[iexp]= sqrt(crossprod(t(newX), P %*% t(newX)) + sig^2)
          }
      }
  names(liste_kf_pred_mean)=paste("kf_pred_mean_",names.experts,sep="")
  names(liste_kf_pred_sd)=paste("kf_pred_sd_",names.experts,sep="")

  #biased experts put to sleep by default
  if(length(which(colnames(awake)=="Q10"))>0) awake[t,"Q10"]=0
  if(length(which(colnames(awake)=="Q30"))>0) awake[t,"Q30"]=0
  if(length(which(colnames(awake)=="Q70"))>0) awake[t,"Q70"]=0
  if(length(which(colnames(awake)=="Q90"))>0) awake[t,"Q90"]=0


  p <- awake[t, ] * w/sum(awake[t, ] * w)
  ref_p <- ref_awake[t, ] * ref_w/sum(ref_awake[t, ] *
      ref_w)
  dataFirstLeadTime = getFirstLeadTime(sta, res, date.valid, ech = 3)
  gradDepth = transObject_threshold$gradDepth
  newDataThreshold = arrangeData(experts_for_xgboost[t,], idx.na, pred = ref_experts[t, ] %*% ref_p, transObject_threshold$training_data,
      gradDepth, eV, liste_kf_pred_sd, liste_kf_pred_mean,
      dataFirstLeadTime, beta, runVar)
  shapleyValues = predict(transObject_threshold, newDataThreshold)
  predErr = sum(shapleyValues)
  if (iteration <= 100) shapleyValues = NULL
  pastErr = (experts[t, ] %*% p) - y


  if (abs(predErr)>=threshold) {#transition probabilities update
      if (predErr>0) { #sum(is_na_spe[t,c("Q10","Q30")])!=2 si au moins un des deux quantilles n'est pas na
          print("on reveille les experts endormis froid")
          print(sta)
          print(iteration)
          s1=asleep1 #if specialized expert is not missing
          s2=asleep1 #if specialized expert is not missing
          s1=1
          s2=1
          if(is_na_spe[t,c("Q10")]) s1=0 #if specialized expert is missing
          if(is_na_spe[t,c("Q30")]) s2=0 #if specialized expert is missing
          if(length(which(colnames(awake)=="Q10"))>0) awake[t,"Q10"]=s1
          if(length(which(colnames(awake)=="Q30"))>0) awake[t,"Q30"]=s2
          if(length(which(colnames(awake)=="Q70"))>0) awake[t,"Q70"]=0
          if(length(which(colnames(awake)=="Q90"))>0) awake[t,"Q90"]=0
      }else if (predErr<0) {#sum(is_na_spe[t,c("Q10","Q30")])!=2 at least one of the quantiles is not NA
          print("on reveille les experts endormis chaud")
          print(sta)
          print(iteration)
          s1=asleep2 #if specialized expert is not missing
          s2=asleep2 #if specialized expert is not missing
          s1=1
          s2=1
          if(is_na_spe[t,c("Q90")]) s1=0 #if specialized expert is missing
          if(is_na_spe[t,c("Q70")]) s2=0 #if specialized expert is missing
          if(length(which(colnames(awake)=="Q10"))>0) awake[t,"Q10"]=0
          if(length(which(colnames(awake)=="Q30"))>0) awake[t,"Q30"]=0
          if(length(which(colnames(awake)=="Q70"))>0) awake[t,"Q70"]=s2
          if(length(which(colnames(awake)=="Q90"))>0) awake[t,"Q90"]=s1
      }
    print(paste("predErr : ",predErr,sep=""))
    print(paste("err : ",pastErr,sep=""))
    # Weight update
    w <- awake[t, ] * w
    p <- w / sum(w)
    newErr=(experts[t, ] %*% p) - y[t]
    print(paste("new err : ",newErr,sep=""))

    if (predErr>0) {
      new_vs_old= abs(pastErr) - abs(newErr) #1
      a1=min(a1,new_vs_old,na.rm=TRUE) #borne inferieure de new_vs_old
      b1=max(b1,new_vs_old,na.rm=TRUE) #borne superieure de new_vs_old
      score_new_vs_old1= (new_vs_old - a1) / (b1-a1) #1/(1+exp(-new_vs_old))
    } else {
      new_vs_old= abs(pastErr) - abs(newErr) #1
      a2=min(a2,new_vs_old,na.rm=TRUE) #borne inferieure de new_vs_old
      b2=max(b2,new_vs_old,na.rm=TRUE) #borne superieure de new_vs_old
      score_new_vs_old2= (new_vs_old - a2) / (b2-a2) #1/(1+exp(-new_vs_old))
    }
  }
  

  #$$
  if(iteration>=100){ #if we start to predict the error (almost no chance that the prediction is equal to zero)
    if(pastErr>=threshold & predErr>=threshold) A1=A1+1
    if(pastErr<=-threshold & predErr<=-threshold) A2=A2+1
    if(pastErr>=threshold & predErr<threshold) B1=B1+1
    if(pastErr<=-threshold & predErr>-threshold) B2=B2+1
    if(pastErr<threshold & predErr>=threshold) C1=C1+1
    if(pastErr>-threshold & predErr<=-threshold) C2=C2+1
    if(pastErr<threshold & predErr<threshold) D1=D1+1
    if(pastErr>-threshold & predErr>-threshold) D2=D2+1
  }

  pred <- experts[t, ] %*% p
  # if (abs(predErr)>=threshold)pred=y[t] #oracle to do like if the specialized experts where perfect when they are used
  ref_pred <- ref_experts[t, ] %*% ref_p
  weights[t, ] <- p
  prediction[t] <- pred

  #loss of the prediction
  lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
  ref_lpred <- lossPred(ref_pred, y[t], ref_pred, loss.type = loss.type, loss.gradient = loss.gradient)#$$
  #vector, losses for each of the M experts
  lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
  ref_lexp <- lossPred(ref_experts[t, ], y[t], ref_pred, loss.type = loss.type, loss.gradient = loss.gradient)


  # Vector, instantaneous regret
  r <-  awake[t, ] * (lpred - lexp) # -l_{j,t} in wintenberger
  ref_r <-  ref_awake[t, ] * (ref_lpred - ref_lexp)#$$

  # Update the learning rates
  V <- V + 2.2*r^2 #vector, cumulated quadratic loss for each experts
  ref_V <- ref_V + 2.2*ref_r^2#$$
  if(length(which(V!=0))==N) {
    eta[t + 1, ] <- pmin(sqrt(1/V),exp(350)) 
  } else {
    eta[t + 1, which(V!=0)]= pmin(sqrt(1/V),exp(350)) [which(V>0)]
    eta[t + 1, which(V==0)]= 0
  }
  #$$
  if(length(which(ref_V!=0))==N) {
    ref_eta[t + 1, ] <- pmin(sqrt(1/ref_V),exp(350)) 
  } else {
    ref_eta[t + 1, which(ref_V!=0)]= pmin(sqrt(1/ref_V),exp(350)) [which(ref_V>0)]
    ref_eta[t + 1, which(ref_V==0)]= 0
  }
  r.reg <- r - eta[t + 1, ] * r^2
  ref_r.reg <- ref_r - ref_eta[t + 1, ] * ref_r^2
  R <- R + r
  R.reg <- R.reg + r.reg
  ref_R <- ref_R + ref_r#$$
  ref_R.reg <- ref_R.reg + ref_r.reg#$$

  #number of times GBRT wakes up the wrong experts
  if ((predErr >= threshold && lpred <= -threshold) || (predErr <= -threshold && lpred >= threshold)) NERR=NERR+1

  w <- truncate1(exp(log(eta[t + 1, ])+log(w0) + eta[t + 1, ] * R.reg))
  w<-w/sum(w)
  ref_w <- truncate1(ref_eta[t + 1, ]*exp(log(ref_eta[t + 1, ])+log(ref_w0) + ref_eta[t + 1, ] * ref_R.reg))
  ref_w<-ref_w/sum(ref_w)

  #kalman update
  newY=y[t]
  Y=c(Y,newY)
  for (iexp in 1:N){
    if(t>1)last_loss=last_l[iexp] else last_loss=0
    newX=matrix(c(1,experts[t,iexp],last_loss),1,3)
    X[[iexp]]=rbind(X[[iexp]],newX)
    d=ncol(X[[iexp]])
    liste_kf[[iexp]]=statespace(X[[iexp]],Y,kalman_params = list(theta=matrix(1,d,1),P=diag(d),Q=diag(d),sig=1),online=TRUE) #kalman filter update with the new observation
  }
  eV = eUpdate(eV, lexp = lossPred(experts[t, ], y[t],
      pred[t], loss.type, loss.gradient = FALSE), lpred = lossPred(pred[t],
      y[t], pred[t], loss.type, loss.gradient = FALSE),
      newexperts = experts, newy = y[t], idx.na = idx.na)

  transObject_threshold = update(transObject=transObject_threshold,
      new_training_data=newDataThreshold, loss = (ref_pred - y[t]), threshold=xgboost_threshold)
  # # $$
  # toto=cbind(transObject_threshold$training_data,sta,ech)
  # if( ! dir.exists(paste("/scratch/work/pfitznerl/training_agreg/",sta,sep=""))) dir.create(paste("/scratch/work/pfitznerl/training_agreg/",sta,sep=""))
  # if( ! dir.exists(paste("/scratch/work/pfitznerl/training_agreg/",sta,"/",ech,sep=""))) dir.create(paste("/scratch/work/pfitznerl/training_agreg/",sta,"/",ech,sep=""))
  # write.table(toto,paste("/scratch/work/pfitznerl/training_agreg/",sta,"/",ech,"/",sta,"_",ech,"_",iteration,
  #     "_1253.txt",sep=""), sep=";",col.names = TRUE, row.names = FALSE)
  # # $$
  }
  object <- list(model = "BOA", loss.type = loss.type, loss.gradient = loss.gradient,
      coefficients = w/sum(w))
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  awake_transition = rbind(training$awake, awake[t, ])
  object$awake_transition <- awake_transition
  AllshapleyValues = rbind.fill(training$shapleyValues, as.data.frame(shapleyValues))
  object$training <- list(eta = eta[T + 1, ], R = R, w0 = w0,
      R.reg = R.reg, V = V, r.reg = rbind(training$r.reg, r.reg),
      experts = rbind(training$experts, experts), y = c(training$y,
          y), prediction = c(training$prediction, prediction),
      r = rbind(training$r, r), eV = eV, transObject_threshold = transObject_threshold,
      sta = sta, ech = ech, window = window,
      awake_transition = awake_transition, liste_kf = liste_kf,
      last_l = lexp, X = X, Y = Y, lnSpread = lnSpread, A1 = A1,
      B1 = B1, C1 = C1, D1 = D1, A2 = A2, B2 = B2, C2 = C2,
      D2 = D2, ref_eta = ref_eta[T + 1, ], ref_R = ref_R, ref_w0 = ref_w0,
      ref_R.reg = ref_R.reg, ref_V = ref_V, ref_r.reg = rbind(training$ref_r.reg,
          ref_r.reg), ref_r = rbind(training$ref_r, ref_r),
      shapleyValues = AllshapleyValues, a1 = a1, b1 = b1, score_new_vs_old1 = c(training$score_new_vs_old1,
          score_new_vs_old1), score_new_vs_old2 = c(training$score_new_vs_old2,
          score_new_vs_old2), a2 = a2, b2 = b2, NERR = NERR)
  class(object) <- "mixture"
  return(object)
}



#SOCO version of BOA with sleeping experts awaken by transition object
# when training data is already in a txt file cf ££ in other scripts
# needs transition_learning_xgboost_multiple_ech_sta
totoBOA <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                w0 = NULL, training = NULL) {
  experts <- as.matrix(experts)
  names.experts=colnames(experts)
  experts_for_xgboost=experts #to train xgboost with real data (without putting missing experts or observations equal to zero)

  N <- ncol(experts)
  T <- nrow(experts)

  eV=eValue(N=N)
  threshold=0.5 # to wake up the experts
  xgboost_threshold=0.5
  transObject_threshold=probTransition(names.experts=names.experts,predType="predErr",nTrain=365*3*18*33*5,
                                        threshold=xgboost_threshold) #probTransition object to predict if the agregation will have a loss bigger than a threshold

  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
      w0 <- rep(1, N)
  }
  #dataToSave=data.frame()
  awake <- as.matrix(awake)
  is_na_spe=is.na(as.data.frame(experts) %>% select(contains(c("Q10","Q30","Q70","Q90")))) #matix true false if spe experts are na
    idx.na <- which(is.na(experts))
    awake[idx.na] <- 0
    experts[idx.na] <- 0

  if (length(which(is.na(y)))>0){ #to deal with missing observations =>no loss for every expert
    id_na_y=which(is.na(y))
    awake[id_na_y,] <- 1 #to avoid 0 division
    experts[id_na_y,] <- 0
    y[id_na_y]<- 0 #has to be done after id_na_y=which(is.na(y))
    }

  colnames(awake)=names.experts

  liste_kf=mylist <- sapply(paste("kf_",names.experts,sep=""),function(x) NULL) #list, will contain the kalman filters for each expert
  X=mylist <- sapply(paste("X_",names.experts,sep=""),function(x) NULL) #list, will contain the predictors for the kalman filter
  Y=c() #list, will contain the observations for the kalman filter

  R <- rep(0, N)
  R.reg <- R
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  w <- w0
  V <- rep(0,N)
  eta <- matrix(0, ncol = N, nrow = T + 1)
  ref_experts <- as.data.frame(experts) %>% select(-contains(c("Q10","Q30","Q70","Q90")))
  ref_experts<-as.matrix(ref_experts)
  ref_N=ncol(ref_experts)
  ref_awake <- as.data.frame(awake) %>% select(-contains(c("Q10","Q30","Q70","Q90")))
  ref_awake<-as.matrix(ref_awake)
  ref_w0=rep(1,ref_N)#$$
  ref_w=ref_w0#$$
    ref_eta <- matrix(0, ncol = ref_N, nrow = T + 1)
  ref_V <- rep(0,ref_N)#$$
  ref_R <- rep(0, ref_N)#$$
  ref_R.reg <- ref_R#$$
  shapleyValues=data.frame() #$$
  NERR=0#$$ count of the large mistake of GBRT
  A1=0
  B1=0
  C1=0
  D1=0
  A2=0
  B2=0
  C2=0
  D2=0
  a1=-1
  b1=1
  a2=-1
  b2=1
  asleep1=0.5
  asleep2=0.5
  score_new_vs_old1=NULL #we only add a score if the sleeping experts are woken up
  score_new_vs_old2=NULL #we only add a score if the sleeping experts are woken up

  lnSpread<-c()#will contain the log of the expert's variance
  if (!is.null(training$R)) {
      w0 <- training$w0
      R <- training$R
      R.reg <- training$R.reg
      w <- truncate1(exp(log(training$eta)+log(w0) + training$eta * R.reg))
      eta[1, ] <- training$eta
      V <- training$V
      ref_w0 <- training$ref_w0
      ref_R <- training$ref_R
      ref_R.reg <- training$ref_R.reg
      ref_w <- truncate1(exp(log(training$ref_eta)+log(ref_w0) + training$ref_eta * ref_R.reg))
      ref_eta[1, ] <- training$ref_eta
      ref_V <- training$ref_V
      eV <- training$eV
      transObject_threshold <- training$transObject_threshold
      liste_kf <- training$liste_kf
      last_l <- training$last_l
      X <- training$X #for kalman filter
      Y <- training$Y #for kalman filter
      lnSpread <- training$lnSpread
      A1 <- training$A1
      B1 <- training$B1
      C1 <- training$C1
      D1 <- training$D1
      A2 <- training$A2
      B2 <- training$B2
      C2 <- training$C2
      D2 <- training$D2
      shapleyValues <- training$shapleyValues
      asleep1<-mean(training$score_new_vs_old1)
      if(is.na(asleep1))asleep1=0.5
      asleep2<-mean(training$score_new_vs_old2)
      if(is.na(asleep2))asleep2=0.5
      a1<-training$a1
      b1<-training$b1
      a2<-training$a2
      b2<-training$b2
      NERR<-training$NERR
      newDataThreshold<-training$newDataThreshold
      #dataToSave=training$dataToSave
  }
  date.valid<-training$date.valid#$$
  sta<-training$sta#$$
  ech<-training$ech#$$
  window<-training$window#$$
  iteration<-training$iteration#$$

  for (t in 1:T) {
  if(iteration==1) print("sleeping BOA")



  #biased experts put to sleep by default
  if(length(which(colnames(awake)=="Q10"))>0) awake[t,"Q10"]=0
  if(length(which(colnames(awake)=="Q30"))>0) awake[t,"Q30"]=0
  if(length(which(colnames(awake)=="Q70"))>0) awake[t,"Q70"]=0
  if(length(which(colnames(awake)=="Q90"))>0) awake[t,"Q90"]=0


  p <- awake[t, ] * w/sum(awake[t, ] * w)
  ref_p <- ref_awake[t, ] * ref_w/sum(ref_awake[t, ] *
      ref_w)

  shapleyValues = predict(transObject_threshold, newDataThreshold)
  predErr = sum(shapleyValues)
  if (iteration <= 100) shapleyValues = NULL
  pastErr = (experts[t, ] %*% p) - y

  if (abs(predErr)>=threshold) {#transition probabilities update
      if (predErr>0) { #sum(is_na_spe[t,c("Q10","Q30")])!=2 si au moins un des deux quantilles n'est pas na
          print("on reveille les experts endormis froid")
          print(sta)
          print(iteration)
          s1=asleep1 #if specialized expert is not missing
          s2=asleep1 #if specialized expert is not missing
          s1=1
          s2=1
          if(is_na_spe[t,c("Q10")]) s1=0 #if specialized expert is missing
          if(is_na_spe[t,c("Q30")]) s2=0 #if specialized expert is missing
          if(length(which(colnames(awake)=="Q10"))>0) awake[t,"Q10"]=s1
          if(length(which(colnames(awake)=="Q30"))>0) awake[t,"Q30"]=s2
          if(length(which(colnames(awake)=="Q70"))>0) awake[t,"Q70"]=0
          if(length(which(colnames(awake)=="Q90"))>0) awake[t,"Q90"]=0
      }else if (predErr<0) {#sum(is_na_spe[t,c("Q10","Q30")])!=2 at least one quantile is not NA
          print("on reveille les experts endormis chaud")
          print(sta)
          print(iteration)
          s1=asleep2 #if specialized expert is not missing
          s2=asleep2 #if specialized expert is not missing
          s1=1
          s2=1
          if(is_na_spe[t,c("Q90")]) s1=0 #if specialized expert is missing
          if(is_na_spe[t,c("Q70")]) s2=0 #if specialized expert is missing
          if(length(which(colnames(awake)=="Q10"))>0) awake[t,"Q10"]=0
          if(length(which(colnames(awake)=="Q30"))>0) awake[t,"Q30"]=0
          if(length(which(colnames(awake)=="Q70"))>0) awake[t,"Q70"]=s2
          if(length(which(colnames(awake)=="Q90"))>0) awake[t,"Q90"]=s1
      }
    print(paste("predErr : ",predErr,sep=""))
    print(paste("err : ",pastErr,sep=""))
    # Weight update
    w <- awake[t, ] * w
    p <- w / sum(w)
    newErr=(experts[t, ] %*% p) - y[t]
    print(paste("new err : ",newErr,sep=""))

    if (predErr>0) {
      new_vs_old= abs(pastErr) - abs(newErr) #1

      a1=min(a1,new_vs_old,na.rm=TRUE) #borne inferieure de new_vs_old
      b1=max(b1,new_vs_old,na.rm=TRUE) #borne superieure de new_vs_old
      score_new_vs_old1= (new_vs_old - a1) / (b1-a1) #1/(1+exp(-new_vs_old))
    } else {
      new_vs_old= abs(pastErr) - abs(newErr) #1
      a2=min(a2,new_vs_old,na.rm=TRUE) #borne inferieure de new_vs_old
      b2=max(b2,new_vs_old,na.rm=TRUE) #borne superieure de new_vs_old
      score_new_vs_old2= (new_vs_old - a2) / (b2-a2) #1/(1+exp(-new_vs_old))

    }

  }

  #$$
  if(iteration>=100){ #if we start to predict the error (almost no chance that the prediction is equal to zero)
    if(pastErr>=threshold & predErr>=threshold) A1=A1+1
    if(pastErr<=-threshold & predErr<=-threshold) A2=A2+1
    if(pastErr>=threshold & predErr<threshold) B1=B1+1
    if(pastErr<=-threshold & predErr>-threshold) B2=B2+1
    if(pastErr<threshold & predErr>=threshold) C1=C1+1
    if(pastErr>-threshold & predErr<=-threshold) C2=C2+1
    if(pastErr<threshold & predErr<threshold) D1=D1+1
    if(pastErr>-threshold & predErr>-threshold) D2=D2+1
  }

  pred <- experts[t, ] %*% p
  #if (abs(predErr)>=threshold)pred=y[t] #oracle to do like if the specialized experts where perfect when they are used
  ref_pred <- ref_experts[t, ] %*% ref_p
  weights[t, ] <- p
  prediction[t] <- pred

  #loss of the prediction
  lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
  ref_lpred <- lossPred(ref_pred, y[t], ref_pred, loss.type = loss.type, loss.gradient = loss.gradient)#$$
  #vector, losses for each of the M experts
  lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
  ref_lexp <- lossPred(ref_experts[t, ], y[t], ref_pred, loss.type = loss.type, loss.gradient = loss.gradient)


  # Vector, instantaneous regret
  r <-  awake[t, ] * (lpred - lexp) # -l_{j,t} in wintenberger
  ref_r <-  ref_awake[t, ] * (ref_lpred - ref_lexp)#$$

  # Update the learning rates
  V <- V + 2.2*r^2 #vector, cumulated quadratic loss for each experts
  ref_V <- ref_V + 2.2*ref_r^2#$$
  if(length(which(V!=0))==N) {
    eta[t + 1, ] <- pmin(sqrt(1/V),exp(350)) 
  } else {
    eta[t + 1, which(V!=0)]= pmin(sqrt(1/V),exp(350)) [which(V>0)]
    eta[t + 1, which(V==0)]= 0
  }
  #$$
  if(length(which(ref_V!=0))==N) {
    ref_eta[t + 1, ] <- pmin(sqrt(1/ref_V),exp(350)) 
  } else {
    ref_eta[t + 1, which(ref_V!=0)]= pmin(sqrt(1/ref_V),exp(350)) [which(ref_V>0)]
    ref_eta[t + 1, which(ref_V==0)]= 0
  }
  r.reg <- r - eta[t + 1, ] * r^2
  ref_r.reg <- ref_r - ref_eta[t + 1, ] * ref_r^2
  R <- R + r
  R.reg <- R.reg + r.reg
  ref_R <- ref_R + ref_r#$$
  ref_R.reg <- ref_R.reg + ref_r.reg#$$

  #number of times GBRT wakes up the wrong experts
  if ((predErr >= threshold && lpred <= -threshold) || (predErr <= -threshold && lpred >= threshold)) NERR=NERR+1

  w <- truncate1(exp(log(eta[t + 1, ])+log(w0) + eta[t + 1, ] * R.reg))
  w<-w/sum(w)
  ref_w <- truncate1(ref_eta[t + 1, ]*exp(log(ref_eta[t + 1, ])+log(ref_w0) + ref_eta[t + 1, ] * ref_R.reg))
  ref_w<-ref_w/sum(ref_w)


  # ££
  df=data_training_all[which(data_training_all$iter<=iteration & data_training_all$sta==sta),] # uniquement obs avant t
  df=subset(df, select = -c(iter)) # on enleve colonne iter
  oversample=df[which(abs(df$lossAgreg)>=xgboost_threshold),]
  df=rbind(df,oversample,oversample)
  y_lossAgreg=df$lossAgreg
  df=subset(df, select = -c(lossAgreg,sta))
  X_mat <- as(as.matrix(df),"dgCMatrix")
  data_training_all_iter <- xgb.DMatrix(data=X_mat,label=y_lossAgreg)

  if(!is.null(length(which(data_training_all$iter==(iteration+1) &
                          data_training_all$sta==sta &
                          data_training_all$ech==ech)))) {
    df1=data_training_all[which(data_training_all$sta==sta &
                                data_training_all$ech==ech &
                                data_training_all$iter==(iteration+1)),] # next obs
    df1=subset(df1, select = -c(iter,lossAgreg,sta)) # remove column iter and y
    X_mat <- as(as.matrix(df1),"dgCMatrix")
    newDataThreshold <- xgb.DMatrix(data=X_mat)
  }
 
  # # ££
  df=data_training_all[which(data_training_all$iter<=iteration),] # uniquement obs avant t
  df=subset(df, select = -c(iter)) # on enleve colonne iter
  df$sta=as.factor(df$sta)
  oversample=df[which(abs(df$lossAgreg)>=xgboost_threshold),]
  df=rbind(df,oversample,oversample)
  y_lossAgreg=df$lossAgreg
  df=subset(df, select = -c(lossAgreg))
  dummies <- dummyVars(~sta,data=df)
  cat_matrix <- predict(dummies,newdata=df)
  X <- cbind(subset(df, select = -c(sta)),cat_matrix)
  X_mat <- as(as.matrix(X),"dgCMatrix")
  data_training_all_iter <- xgb.DMatrix(data=X_mat,label=y_lossAgreg)

  # if(!is.null(length(which(data_training_all$iter==(iteration+1) &
  #                         data_training_all$sta==sta &
  #                         data_training_all$ech==ech)))) {
  #   df1=data_training_all[which(data_training_all$sta==sta &
  #                               data_training_all$ech==ech &
  #                               data_training_all$iter==(iteration+1)),] # next obs
  #   df1=subset(df1, select = -c(iter,lossAgreg)) # remove column iter and y
  #   df1$sta=factor(df1$sta,levels=levels(df$sta))
  #   dummies <- dummyVars(~sta,data=df1)
  #   cat_matrix <- predict(dummies,newdata=df1)
  #   X <- cbind(subset(df1, select = -c(sta)),cat_matrix)
  #   X_mat <- as(as.matrix(X),"dgCMatrix")
  #   newDataThreshold <- xgb.DMatrix(data=X_mat)
  # }
  # # ££

  transObject_threshold = update(transObject=transObject_threshold,
      xgb_training=data_training_all_iter, loss = (ref_pred - y[t]), threshold=xgboost_threshold)

  }
  object <- list(model = "BOA", loss.type = loss.type, loss.gradient = loss.gradient,
      coefficients = w/sum(w))
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  awake_transition = rbind(training$awake, awake[t, ])
  object$awake_transition <- awake_transition
  AllshapleyValues = rbind.fill(training$shapleyValues, as.data.frame(shapleyValues))
  object$training <- list(eta = eta[T + 1, ], R = R, w0 = w0,
      R.reg = R.reg, V = V, r.reg = rbind(training$r.reg, r.reg),
      experts = rbind(training$experts, experts), y = c(training$y,
          y), prediction = c(training$prediction, prediction),
      r = rbind(training$r, r), eV = eV, transObject_threshold = transObject_threshold,
      sta = sta, ech = ech, window = window,
      awake_transition = awake_transition, liste_kf = liste_kf,
      last_l = lexp, X = X, Y = Y, lnSpread = lnSpread, A1 = A1,
      B1 = B1, C1 = C1, D1 = D1, A2 = A2, B2 = B2, C2 = C2,
      D2 = D2, ref_eta = ref_eta[T + 1, ], ref_R = ref_R, ref_w0 = ref_w0,
      ref_R.reg = ref_R.reg, ref_V = ref_V, ref_r.reg = rbind(training$ref_r.reg,
          ref_r.reg), ref_r = rbind(training$ref_r, ref_r),
      shapleyValues = AllshapleyValues, a1 = a1, b1 = b1, score_new_vs_old1 = c(training$score_new_vs_old1,
          score_new_vs_old1), score_new_vs_old2 = c(training$score_new_vs_old2,
          score_new_vs_old2), a2 = a2, b2 = b2, NERR = NERR,
      newDataThreshold=newDataThreshold)
  class(object) <- "mixture"
  return(object)
}