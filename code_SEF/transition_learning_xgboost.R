library(xgboost)
library(caret)
library(Matrix)

##########################################################################
##########################################################################
# a constructor/initialization function for the "probTransition" class, which enables to
#predict a probability vector of the experts lossews or the order of the expert losses or which will be the best expert
#training_data: data to train the model
#nexp: number of experts
#names.experts : names of the experts
#RF: random forest to predict the probability vector of which will be the best expert
#predType : what to predict :
#                           - loss, loss of all the experts,
#                           - best, probability vector for every expert of being the best,
#                           - order, order of the expert depending on there loss,
#                           - predErr, predict if the loss of the agregation will be bigger than a threshold
#nTrain: the maximum number of past iterations to train the model
#maxDist: if flag.distance=TRUE, then if the distance between the training data and the data that has to be predicted is bigger than maxDist,
#we don't use the probTransition prediction, default is infinity in order to always use the model for every prediction
#normType: is the type of norm which will be used to compute the distance between the training data and the new data,
#infinity, one and two norms
#modelType, type of the machine learning algorithm to use for learning and predicting the data, knn or rf (random forest)
#minTrain, minimum number of training data to use the models
# threshold, to predict for predType predErr and also used to raise weight of some observations via case.wt
# case.wt weight vector, Vector of non-negative weights (does not have to sum to 1) for sampling cases.
# Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples.
# gradDepth, to compute and add the gradients with respect to 1,...,gradDepth day before
# varNames names of the variables that have effectively been used for the training (some of the columns may
# have been removed (cf rm_missing_col() in modelXXXupdate functions) if to much data was missing)
probTransition <- function(names.experts,nTrain=100,training_data=data.frame(),modelType="rf",
                            predType=NULL,maxDist=exp(999),normType="infinity",minTrain=100,
                            threshold=2.5,case.wt=NULL,flag.distance=FALSE,gradDepth=3,varNames=NULL) {
    # we can add our own integrity checks
    if(predType!="loss" & predType!="best" & predType!="order" & predType!="predErr") stop("unknown type of prediction")
    
    if (predType=="predErr"  & (is.null(threshold)) ) stop("you have to give a threshold when using predErr")

    if (nrow(training_data)!=0) {
        if (modelType=="rf") {
            stop("you should train with the available data")
            # model <-rfsrc(c("l.raw.aro","l.mos.aro","l.raw.arp","l.mos.arp","l.raw.cep",
            #                         "l.mos.cep","l.Q10","l.Q30","l.Q50","l.Q70","l.Q90")~.,
            #                         data=training_data,na.action="na.omit")
        } else if (modelType=="knn") {
            model <- "toto fait quelque chose !"
        }
    } else model<- NULL

    nexp=length(names.experts)

    object <- list(training_data = training_data, nexp=nexp, model=model, modelType=modelType,predType=predType,
                    nTrain=nTrain,maxDist=maxDist,normType=normType,minTrain=minTrain,threshold=threshold,
                    names.experts=names.experts,case.wt=case.wt,flag.distance=flag.distance,gradDepth=gradDepth)
    class(object)<- "probTransition"
    return(object)
}
##########################################################################
##########################################################################



##########################################################################
##########################################################################

#updates objects of the class probTransition, which enables to predict a probability vector of which will be the best expert
#transObject, the probTransition object with the random forest model which has to be updated
#new_training_data, the new data of the current iteration to train the model
# loss, loss of the agregation if predType="predErr"
update.probTransition <- function(transObject,new_training_data,loss=NULL,threshold){
    predType=transObject$predType
    cls <- sapply(new_training_data, class) #type of the columns

    for (i in which(cls=="character")) new_training_data[,i]=as.factor(new_training_data[,i]) #character to factor

    if(predType=="loss") {
        transObject <- modelLossUpdate(transObject,new_training_data,loss)
    } else if(predType=="best") {
        transObject <- modelBestUpdate(transObject,new_training_data)
    } else if(predType=="order") {
        transObject <- modelOrderUpdate(transObject,new_training_data)
    } else if(predType=="predErr") {
        transObject <- modelpredErrUpdate(transObject,new_training_data,loss,threshold)
    } else {
        stop("problem with predType")
    }
    return(transObject)
}

#entry : probTransition object and the new data
#output : pred, the prediction of the random forest, boolean for predErr, matrix otherwise
predict.probTransition <- function(transObject,newData,weights=NULL){

    nexp=transObject$nexp
    names.experts=transObject$names.experts
    flag.distance=transObject$flag.distance
    predType=transObject$predType
    nTrain=transObject$nTrain
    maxDist=transObject$maxDist
    modelType=transObject$modelType
    minTrain=transObject$minTrain
    nbTrain=nrow(transObject$training_data) #number of training observations
    if (nbTrain>=nTrain & flag.distance) { #distance only needed when there are enough training observations to use the prediction model
        distance=dist_pred_train(transObject$training_data[,names.experts],newData,transObject$normType)#distance between data to be predicted and the training data
    } else distance=0
    #whe only use the model prediction when at least nTrain iterations trained it,
    #and when the newData is similar (distance<=maxDist) to some of the training data
    if (nbTrain>=minTrain & distance<= maxDist) {
        if (!is.null(transObject$model)) {
            if (modelType=="rf"){

                #to match the variables of training_data, because some of the columns may have been removed in the
                #training (cf rm_missing_col() in modelXXXupdate functions) if to much data was missing
                varNames <- names(newData)[(names(newData) %in% transObject$varNames)]
                newData=newData[varNames]
                newData=as.numeric(newData) #to transform into a vector (as.vector() does not work)
                newData <- as(newData, "sparseVector") # to transform into a sparseVector such that xgboost can use it

                if (predType!="predErr") {
                    newData=newData[names(newData) != "pred"]
                    #theta share to best expert
                    lossPred <- predict(transObject$model,newData) #get the the loss prediction
                    ibestPred=which(abs(lossPred)==min(abs(lossPred)))
                    #ibestWeight=which(c(weights)==max(c(weights)))
                    if (length(ibestPred)>1) print("the RF predicts that there will be several best expert")
                    transition_prob_matrix=matrix(0, nrow = nexp,ncol = nexp, byrow = T)
                    #transition_prob_matrix[ibestPred,]=0.05
                    lossRank=rank(abs(lossPred),ties.method = "min")
                    transition_prob_matrix[ibestPred,]=((lossRank-1)/100)/length(ibestPred)
                    transition_prob_matrix=transition_prob_matrix+(1-((lossRank-1)/100))*diag(nexp)

                    pred=transition_prob_matrix/apply(transition_prob_matrix,2,sum)

                } else { #if predType=="predErr"
                    # prediction.model <- predict(transObject$model,newData)
                    #pred <- prediction.model$predicted #get the prediction of the agregation's error for rfsrc
                    #pred <- predict(transObject$model,newData)
                    #to have the shapley values/contribution of each feature
                    pred <- predict(transObject$model,newData, predcontrib = TRUE)# , predinteraction =TRUE)
                    #the last column of the shapley values is not lossAgreg but the bias,
                    # cf https://cran.r-project.org/web/packages/xgboost/xgboost.pdf#Rfn.xgb.ggplot.shap.summary :
                    # "The last "+ 1" column in a matrix corresponds to bias
                    colnames(pred)=c(varNames,"bias")
                }
            } else if (modelType=="knn"){
                transition_prob_vector="toto fait quelque chose"
            }
        } else stop("the random forest is empty")
    } else if(predType=="best") {
        pred=diag(nexp)
    } else if(predType=="loss") {
        pred=diag(nexp)
        #pred=rep(1/nexp,nexp)
    } else if(predType=="order") {
        pred=diag(nexp)
    } else if(predType=="predErr") {
        pred=0
    }
    return(pred)
}
##########################################################################
##########################################################################



##########################################################################
##########################################################################

# trains the random forest with the observations of the experts losses
# loss : loss of the aggregation, if >threshold, we raise the weight probibility
#of the observation to be sampled in the bootstrap of the rf
# RF to predict the loss of the experts
modelLossUpdate <- function(transObject,new_training_data,loss) {
    names(loss)=paste(transObject$names.experts,"Loss",sep="")
    toPredict=paste(transObject$names.experts,"Loss",sep="")
    #new_training_data=as.data.frame(new_training_data)
    new_training_data=cbind(new_training_data,t(loss))
    nTrain=transObject$nTrain
    modelType=transObject$modelType

    if(!is.null(transObject$training_data)) {ndata=nrow(transObject$training_data)} else ndata=0 #number of observations

    # #boosting, more weight on all the observation of all lead times of the run where the aggregation had a loss>threshold
    # bthreshold=abs(loss)>transObject$threshold
    # if (!bthreshold){transObject$case.wt=c(transObject$case.wt,rep(1,nrow(new_training_data))) #if small loss
    # } else transObject$case.wt=c(transObject$case.wt,rep(4,nrow(new_training_data)))# if large loss, more chances to be sampled in the random forest

    transObject$training_data=rbind(transObject$training_data,new_training_data)#we add the new data to the training data

    if (ndata>=nTrain) transObject$training_data=tail(transObject$training_data,n=nTrain) #we only keep nTrain data
    training_data=rm_missing_col(transObject$training_data,0.5) #remove variables with to much missing data

    #model training
    if (modelType=="rf") {
        #when using get.mv.formula, we don't need "(y1,y2,y3)~."
        transObject$model=rfsrc(get.mv.formula(toPredict), data = training_data,
                                na.action = "na.omit")#,ntree=250)#,nodesize=1)#,case.wt=transObject$case.wt)
        # transObject$model=rfsrc.fast(get.mv.formula(colnames(transObject$training_data[,colPred])),
        #     data = transObject$training_data, na.action = "na.omit", forest=TRUE,ntree=60)
    } else if (modelType=="knn") {
        transObject$model="fait quelque chose toto"
    }

    return(transObject)
}

#trains the random forest with the observations of the order of the expert with respect to there losses
modelOrderUpdate <- function(transObject,new_training_data) {
    names.experts=transObject$names.experts
    new_training_data=as.data.frame(new_training_data)
    nTrain=transObject$nTrain
    modelType=transObject$modelType
    names.experts=transObject$names.experts
    if(!is.null(transObject$training_data)) {ndata=nrow(transObject$training_data)} else ndata=0 #number of observations
    bestExp=which_orderExp(loss,names.experts)#vector, 1 if the expert has the lowest loss, 0 otherwise
    new_training_data=cbind(new_training_data,bestExp)#we add the true outcome to new_training_data
    transObject$training_data=rbind(transObject$training_data,new_training_data)#we add the new data to the training data
    if (ndata>=nTrain) transObject$training_data=tail(transObject$training_data,n=nTrain) #we only keep nTrain data
    #training_data=rm_missing_col(transObject$training_data,0.3) #remove variables with to much missing data
    
    #model training
    if (modelType=="rf") {
        transObject$model=rfsrc(get.mv.formula(names(loss)),data = transObject$training_data, na.action = "na.omit")
    } else if (modelType=="knn") {
        transObject$model="fait quelque chose toto"
    }

    return(transObject)
}

#trains the random forest with the observation of the best expert
modelBestUpdate <- function(transObject,new_training_data) {
    names.experts=transObject$names.experts
    new_training_data=as.data.frame(new_training_data)
    nTrain=transObject$nTrain
    modelType=transObject$modelType
    names.experts=transObject$names.experts
    if(!is.null(transObject$training_data)) {ndata=nrow(transObject$training_data)} else ndata=0 #number of observations
    orderExp=which_bestExp(loss,names.experts)#vector, 1 if the expert has the lowest loss, 0 otherwise
    new_training_data=cbind(new_training_data,orderExp)#we add the true outcome to new_training_data
    transObject$training_data=rbind(transObject$training_data,new_training_data)#we add the new data to the training data
    if (ndata>=nTrain) transObject$training_data=tail(transObject$training_data,n=nTrain) #we only keep nTrain data
    
    #model training
    if (modelType=="rf") {
        transObject$model=rfsrc(get.mv.formula(names(loss)),data = transObject$training_data, na.action = "na.omit")
    } else if (modelType=="knn") {
        transObject$model="fait quelque chose toto"
    }

    return(transObject)
}

# trains the random forest with the observations of the experts losses
# loss : difference between prediction and observation
# RF to predict the loss of the agregation
modelpredErrUpdate <- function(transObject,new_training_data,loss,threshold) {
    transObject$threshold=threshold
	lossAgreg=loss
    names.experts=transObject$names.experts
    nTrain=transObject$nTrain
    modelType=transObject$modelType
    minTrain=transObject$minTrain
    if(!is.null(transObject$training_data)) {ndata=nrow(transObject$training_data)} else ndata=0 #number of observations
    new_training_data=cbind(new_training_data,lossAgreg)#we add the true outcome to new_training_data

    transObject$training_data=rbind(transObject$training_data,new_training_data) #training data without oversampling
    if (ndata>=nTrain) transObject$training_data=tail(transObject$training_data,n=nTrain) #we only keep nTrain data

    #object, indiquant coment controler l'entrainement de xgboost, ici cv pour cross validation
    xgb_trcontrol = trainControl(method = "cv", number = 5, allowParallel = TRUE, 
        verboseIter = FALSE, returnData = FALSE)
    #grille de paramètres à tester pour l'apprentissage de xgboost
    # colsample_bytree: pourcentage des colonnes pris pour construire un arbre
    # colsample_bynode is the subsample ratio of columns for each node (split)
    xgbGrid <- expand.grid(nrounds = c(2),#c(100,200),  
                        max_depth = 7, #(1, 3, 5, 10, 15),
                        colsample_bytree = 0.8, #seq(0.5, 0.9, length.out = 5),
                        ## valeurs par défaut : 
                        eta = 0.1,
                        gamma=0,
                        min_child_weight = 5,
                        subsample = 1
                        )

    training_data=transObject$training_data

    oversample=training_data[which(abs(training_data$lossAgreg)>=threshold),]
    training_data=rbind(training_data,oversample,oversample)

    # rm_missing_col has to be after the oversampling to avoid conflicts when a variable was missing and is no more missing
    training_data=rm_missing_col(training_data,0.5) #remove variables with to much missing data
    transObject$varNames=colnames(training_data)

    #model training
    if(ndata>minTrain-10){#to train the RF before we make a prediction
        cls <- sapply(training_data, class) #type of the columns
        iRemove=c()
        for (i in which(cls=="factor")) {if(length(levels(training_data[,i])) < 2 ) iRemove=c(iRemove,i)}
        if (!is.null(iRemove)) training_data=training_data[,-iRemove] #sparse.model.matrix needs at least 2 levels
        #remove rows with NAs, otherwise sparse.model.matrix will remove them and nrow(yTrain)!=nrow(xTrain)
        #so the first iterations are missing because of gradeV which is NA
        training_data=training_data[complete.cases(training_data), ]


        xTrain = sparse.model.matrix(lossAgreg ~ ., data = training_data)[, -1] #convert factors to dummy variables
        yTrain = as.numeric(data.matrix(subset(training_data, select = c(lossAgreg))))

        dtrain=xgb.DMatrix(data = xTrain, label = yTrain)#,weight=w)
        
        if (modelType=="rf") {

            # results_05_    oversample2  zSAVE_PSS_leaves_avec_kf_new_519
            transObject$model=xgboost(data=xTrain,label=yTrain,max_depth = 8,gamma=0, min_child_weight=100,
               eta = 0.4, nthread = 1, nrounds = 7, colsample_bynode = 0.85, subsample=1, objective = "reg:squarederror",verbose = 0)

        } else if (modelType=="knn") {
            stop("seulement RF avec xgboost")
        }
    }

    return(transObject)
}




#gives the best expert of the iteration
which_bestExp <-function(lexp,names.experts){
    nexp=length(lexp)
    ibestExp=which(lexp==min(lexp)) #column of the best expert
    orderExp=rep(0,nexp)
    orderExp[ibestExp]=1 #best expert=1, other experts=0
    orderExp=as.data.frame(t(orderExp))
    names(orderExp)<-paste(names.experts,"Order",sep="")
    return(orderExp)
}

#gives the best expert of the iteration
which_orderExp <-function(lexp,names.experts){
    stop("edit which_orderExp function !!!!")
    nexp=length(lexp)
    ibestExp=which(lexp==min(lexp)) #column of the best expert
    bestExp=rep(0,nexp)
    bestExp[ibestExp]=1 #best expert=1, other experts=0
    bestExp=as.data.frame(t(bestExp))
    names(bestExp)<-paste(names.experts,"Best",sep="")
    return(bestExp)
}

#arrange the data to predict and train with probTransition objects
#experts : prediction of the experts
#pred : prediction of the aggregation for the prediction ; for the learning prediction whith the updated weights  
#oldData : data of the previous iteration
arrangeData <- function(experts,idx.na=NULL,pred=NULL,oldData=NULL,gradDepth=1,eV=NULL,liste_kf_pred_sd=NULL,liste_kf_pred_mean=NULL,
                        dataFirstLeadTime=NULL,beta=NULL,runVar=NULL) {#,eV=NULL) {
    if (!is.null(idx.na)) { experts[idx.na]=NA
    } else warning("WARNING: you didn't specify if they are missing experts => the missing experts are put to zero by opera...!!!")
    names.experts=names(experts)
    nexp=length(names.experts)
    min_T=min(experts,na.rm=TRUE)
    max_T=max(experts,na.rm=TRUE)
    mean_T=round(mean(experts,na.rm=T),2)
    if(min_T==Inf) min_T=NA #too deal with cases where all experts are NA
    if(max_T==-Inf) max_T=NA # too deal with cases where all experts are NA
    if(is.nan(mean_T))mean_T=NA # too deal with cases where all experts are NA
    var_T=round(var(experts,na.rm=TRUE),2) #variance of all the experts
    SD_T= round(sqrt(var(experts,na.rm=TRUE)),2) #standard deviation of all the experts
    
    fortinRMSE=0 #rmse estimation like in fortin 2014
    for (iq in 1:nexp){
        fortinRMSE=fortinRMSE+ (mean(experts)-experts[iq])^2
    }
    fortinRMSE=fortinRMSE/(nexp-1)
    fortinRMSE=round(sqrt(fortinRMSE*(nexp+1)/nexp),2)

    idq=which(substr(names(experts),1,1)=="Q")
    nq=length(idq) #number of experts of the ensemble forecast pearp
    #varQ_T=round(var(experts[idq],na.rm=TRUE),2) #variance of the ensemble forecast pearp
    SDQ_T= round(sqrt(var(experts[idq],na.rm=TRUE)),2) #standard deviation of the esemble forecast pearp

    #difference between the mean of the experts and the aggregation, to link the agregation to the variance !?
    #because fortin 2014 says that variance and MEAN of the ensemble are linked
    diff_agreg_mean=pred-mean_T
    diff_agreg_meanQ=pred-mean(experts[idq],na.rm=TRUE)

    if(length(which(names.experts=="Q30"))>0)SDQ_T_normal=round(experts["Q70"]-experts["Q30"],2) #variance estimator if the distribution is normal
    if(length(which(names.experts=="Q20"))>0)SDQ_T_normal=round(experts["Q70"]-experts["Q20"],2) #TODO verrue !!!!!! it's the Q30 and not the Q20, it would be better to change ths in the code which creates the .txt

    fortinQRMSE=0 #RMSE estimation of the ensemble forecast pearp, like in fortin 2014
    for (iq in idq){
        fortinQRMSE=fortinQRMSE+ (mean(experts[idq])-experts[iq])^2
    }
    fortinQRMSE=fortinQRMSE/(nq-1)
    fortinQRMSE=round(sqrt(fortinQRMSE*(nq+1)/nq),2)

    amp=max_T - min_T

    if ( length(which(names.experts=="mos.aro"))==1 & length(which(names.experts=="mos.arp"))==1) {
        mos.aro_mos.arp=experts[which(names.experts=="mos.aro")]-experts[which(names.experts=="mos.arp")]
    } else mos.aro_mos.arp=NA
    if ( length(which(names.experts=="mos.aro"))==1 & length(which(names.experts=="mos.cep"))==1) {
        mos.aro_mos.cep=experts[which(names.experts=="mos.aro")]-experts[which(names.experts=="mos.cep")]
    } else mos.aro_mos.cep=NA
    if ( length(which(names.experts=="mos.arp"))==1 & length(which(names.experts=="mos.cep"))==1) {
        mos.arp_mos.cep=experts[which(names.experts=="mos.arp")]-experts[which(names.experts=="mos.cep")]
    } else mos.arp_mos.cep=NA
    if ( length(which(names.experts=="Q50"))==1 & length(which(names.experts=="mos.aro"))==1) {
        Q50_mos.aro=experts[which(names.experts=="Q50")]-experts[which(names.experts=="mos.aro")]
    } else Q50_mos.aro=NA
    if ( length(which(names.experts=="Q50"))==1 & length(which(names.experts=="mos.arp"))==1) {
        Q50_mos.arp=experts[which(names.experts=="Q50")]-experts[which(names.experts=="mos.arp")]
    } else Q50_mos.arp=NA
    if ( length(which(names.experts=="Q50"))==1 & length(which(names.experts=="mos.cep"))==1) {
        Q50_mos.cep=experts[which(names.experts=="Q50")]-experts[which(names.experts=="mos.cep")]
    } else Q50_mos.cep=NA
    mean_max=max_T - mean_T
    mean_min=mean_T - min_T
    df_experts=as.data.frame(t(experts))#numeric vector
    df_coVariables=as.data.frame(t(c(SD_T,SDQ_T,
                                    mean_min,mean_max,diff_agreg_meanQ)))#numeric vector
    names(df_coVariables)<-c("SD_T","SDQ_T",
                            "mean_min","mean_max","diff_agreg_meanQ")
    newData=df_coVariables
    
    # vecMinMax=c()#logical vector
    # for (iexp in 1:length(experts)){
    #         vecMinMax=c(vecMinMax,experts[iexp]==min_T)
    #         names(vecMinMax)[length(vecMinMax)]=paste("isMin_",names.experts[iexp],sep="")
    #         vecMinMax=c(vecMinMax,experts[iexp]==max_T)
    #         names(vecMinMax)[length(vecMinMax)]=paste("isMax_",names.experts[iexp],sep="")
    # }
    # newData=cbind(newData,t(vecMinMax))

    # if (!is.null(pred)) newData=cbind(newData,pred)

    if (!is.null(liste_kf_pred_mean)) {
        kf_mean=t(as.numeric(liste_kf_pred_mean[1]))
        colnames(kf_mean)= names(liste_kf_pred_mean)[1]
        newData=cbind(newData,kf_mean)
    }
    
    if (!is.null(liste_kf_pred_sd)) {
        kf_sd=t(as.numeric(liste_kf_pred_sd[1]))
        colnames(kf_sd)= names(liste_kf_pred_sd)[1]
        newData=cbind(newData,kf_sd)
    }

    # if (!is.null(eV)){
    #     nbcol=length(newData)
    #     lastIt=dim(eV$M)[3]
    #     if(eV$grad) {
    #         newData=cbind(newData,t(eV$gradM[2,-c(2),lastIt])) #Ev expert 2 vs j, and remove the eV comparing expert 1 to himself
    #         names(newData)[(nbcol+1):(nbcol+nexp-1)]<-paste("grad.eV_",2,"_",2:nexp,sep="")
    #         # nbcol=length(newData)
    #         # newData=cbind(newData,t(eV$gradM[4,-c(4),lastIt])) #Ev expert 4 vs j, and remove the eV comparing expert 1 to himself
    #         # names(newData)[(nbcol+1):(nbcol+nexp-1)]<-paste("grad.eV_",4,"_",2:nexp,sep="")
    #         # nbcol=length(newData)
    #         # newData=cbind(newData,t(eV$gradM[6,-c(6),lastIt])) #Ev expert 6 vs j, and remove the eV comparing expert 1 to himself
    #         # names(newData)[(nbcol+1):(nbcol+nexp-1)]<-paste("grad.eV_",6,"_",2:nexp,sep="")
    #         # nbcol=length(newData)
    #         # newData=cbind(newData,t(eV$gradM[9,-c(9),lastIt])) #Ev expert 11 vs j, and remove the eV comparing expert 1 to himself
    #         # names(newData)[(nbcol+1):(nbcol+nexp-1)]<-paste("grad.eV_",9,"_",2:nexp,sep="")
    #     } else {
    #         newData=cbind(newData,t(eV$M[1,2:nexp,lastIt])) #Ev expert 1 vs j, and remove the eV comparing expert 1 to himself
    #         names(newData)[(nbcol+1):(nbcol+nexp-1)]<-paste("eV_",1,"_",2:nexp,sep="")
    #     }
    # }

    #newData=addGrad(oldData,newData,gradDepth)

    # if (!is.null(beta)) {
    #     newData=cbind(newData,beta)
    #     names(newData)[length(newData)]="beta"
    #     newData=cbind(newData, SD_T - beta)
    #     names(newData)[length(newData)]="diff_SD"
    # }

    if (!is.null(runVar)) {
        newData=cbind(newData,runVar$runMeanVar) #mean variance during the run
        names(newData)[length(newData)]="runMeanVar"
        newData=cbind(newData,runVar$runVarVar) # variance of the variance during the run
        names(newData)[length(newData)]="runVarVar"
        newData=cbind(newData,var_T-runVar$runMeanVar)
        names(newData)[length(newData)]="diff_runMeanVar"
    }

    #add the first lead time observation and diferrences between experts and observation
    if (!is.null(dataFirstLeadTime)) {
        dataDiff=dataFirstLeadTime - dataFirstLeadTime[,"obs.t"] #difference between experts and observation
        dataFirstLeadTime=cbind(dataFirstLeadTime[,"obs.t"],dataDiff[-c(1)])
        names(dataFirstLeadTime)<-paste("firstLeadT.",c("obs.t",names(dataFirstLeadTime)[-c(1)]),sep="")
        newData=cbind(newData,dataFirstLeadTime)
    }

    return(newData)
}

#recovers training data from different lead times
# flagloss to add the losses of the experts, need to be false if used for threshold !?
recupLeadTimes <- function(ech,sta,res,iteration,names.experts,dataFirstLeadTime,oldData,gradDepth,flagLoss=FALSE){ #},window,agreg){
    names.experts[which(names.experts=="Q30")]="Q20" # TODO verrue !!!!!! it's the Q30 and not the Q20, it would be better to change ths in the code which creates the .txt
    dataLeadTimes=NULL

    iechMax=which(liste_ech_xgboost<=ech) #we don't take the lead times after ech, for simplicity (not the same experts after 48h, and we avoid problems with observation availability)

    runVar=computeRunVar(sta,res,iteration,names.experts) #Variance of the variance in the run, and mean variance

    for (iech in iechMax){ #to take the previous lead times into account
    #for (iech in iechMax){
        file_name=dir(path = dir_data_with, pattern = paste("donnee_",parametre,"_",res,"_",liste_ech_xgboost[iech],"_",sta,sep=""),full.names=TRUE)
        #il y a un ou des espaces avant les NA..., d'ou "  NA"
        donnee = read.table(file_name,head=T,sep=";",na.strings=c("       NA","      NA","     NA","    NA","   NA","  NA"," NA","NA"))
        iRowOld=nrow(oldData)-iech+1 #correspondinf lead time of the previous iteration
        if (is.null(dataLeadTimes)) {
            keep_data=as.numeric(donnee[iteration,names.experts]) #convert to a vector
            idx.na <- which(is.na(keep_data)) #missing experts
            names(keep_data)=names.experts
            keep_data=arrangeData(experts=keep_data,idx.na=idx.na,pred=NULL,oldData=oldData[iRowOld,],gradDepth=gradDepth,
                                    dataFirstLeadTime=dataFirstLeadTime,runVar=runVar)
            expLoss=(donnee[iteration,names.experts]-donnee[iteration,"obs.t"])^2 #square loss of the experts
            colnames(expLoss)=paste(names.experts,"Loss",sep="")
            if (flagLoss) { dataLeadTimes=cbind(keep_data,expLoss)
            } else dataLeadTimes = keep_data
        } else {
            keep_data=as.numeric(donnee[iteration,names.experts]) #convert to a vector
            idx.na <- which(is.na(keep_data)) #missing experts
            names(keep_data)=names.experts
            keep_data=arrangeData(experts=keep_data,idx.na=idx.na,pred=NULL,oldData=oldData[iRowOld,],gradDepth=gradDepth,
                                    dataFirstLeadTime=dataFirstLeadTime,runVar=runVar)
            expLoss=(donnee[iteration,names.experts]-donnee[iteration,"obs.t"])^2 #square loss of the experts
            colnames(expLoss)=paste(names.experts,"Loss",sep="")
            if (flagLoss) keep_data=cbind(keep_data,expLoss)
            dataLeadTimes=rbind(dataLeadTimes,keep_data)
        }
    }
    colnames(dataLeadTimes)[which(colnames(dataLeadTimes)=="Q20")]="Q30" # TODO verrue !!!!!! it's the Q30 and not the Q20, it would be better to change ths in the code which creates the .txt
    colnames(dataLeadTimes)[which(colnames(dataLeadTimes)=="Q20Loss")]="Q30Loss" # TODO verrue !!!!!! it's the Q30 and not the Q20, it would be better to change ths in the code which creates the .txt
    colnames(dataLeadTimes)[which(colnames(dataLeadTimes)=="isMax_Q20")]="isMax_Q30" # TODO verrue !!!!!! it's the Q30 and not the Q20, it would be better to change ths in the code which creates the .txt
    colnames(dataLeadTimes)[which(colnames(dataLeadTimes)=="isMin_Q20")]="isMin_Q30" # TODO verrue !!!!!! it's the Q30 and not the Q20, it would be better to change ths in the code which creates the .txt
    #$$rownames(dataLeadTimes)=lead_times[1:iechMax]
    rownames(dataLeadTimes)=lead_times[iechMax]
    return(dataLeadTimes)
}

#add the gradient of some variables to the training data
#oldData : data of the previous iteration
addGrad <- function(oldData,newData,gradDepth) {
    for (i in 1:gradDepth){#compute and add the gradients with respect to i day before
        oldD=oldData[nrow(oldData)-i+1,] #iteration t-i
        oldNames=names(newData)
        if (nrow(oldD)>0){
            if(!is.null(newData$var_T) & !is.null(oldD$var_T) ) gradVar=newData$var_T - oldD$var_T else gradVar=NA
            if(!is.null(newData$amp) & !is.null(oldD$amp)) gradAmp=newData$amp - oldD$amp else gradAmp=NA
            # if(!is.null(newData$mos.aro_mos.arp) & !is.null(oldD$mos.aro_mos.arp)) {
            #     grad.mos.aro_mos.arp=newData$mos.aro_mos.arp - oldD$mos.aro_mos.arp
            # } else grad.mos.aro_mos.arp=NA
            # if(!is.null(newData$mos.aro_mos.cep) & !is.null(oldD$mos.aro_mos.cep)) {
            #     grad.mos.aro_mos.cep=newData$mos.aro_mos.cep - oldD$mos.aro_mos.cep
            # }else grad.mos.aro_mos.cep=NA
            # if(!is.null(newData$mos.arp_mos.cep) & !is.null(oldD$mos.arp_mos.cep)) {
            #     grad.mos.arp_mos.cep=newData$mos.arp_mos.cep - oldD$mos.arp_mos.cep
            # }else grad.mos.arp_mos.cep=NA

            # if(!is.null(newData$Q50_mos.aro) & !is.null(oldD$Q50_mos.aro)) {
            #     grad.Q50_mos.aro=newData$Q50_mos.aro - oldD$Q50_mos.aro
            # }else grad.Q50_mos.aro=NA
            # if(!is.null(newData$Q50_mos.arp) & !is.null(oldD$Q50_mos.arp)) {
            #     grad.Q50_mos.arp=newData$Q50_mos.arp - oldD$Q50_mos.arp
            # }else grad.Q50_mos.arp=NA
            # if(!is.null(newData$Q50_mos.cep) & !is.null(oldD$Q50_mos.cep)) {
            #     grad.Q50_mos.cep=newData$Q50_mos.cep - oldD$Q50_mos.cep
            # }else grad.Q50_mos.cep=NA


            # if(!is.null(newData$mos.aro) & !is.null(oldD$mos.aro)) {
            #     grad.mos.aro=newData$mos.aro - oldD$mos.aro
            # }else grad.mos.aro=NA
            # if(!is.null(newData$mos.arp) & !is.null(oldD$mos.arp)) {
            #     grad.mos.arp=newData$mos.arp - oldD$mos.arp
            # }else grad.mos.arp=NA
            # if(!is.null(newData$mos.cep) & !is.null(oldD$mos.cep)) {
            #     grad.mos.cep=newData$mos.cep - oldD$mos.cep
            # }else grad.mos.cep=NA
            # if(!is.null(newData$Q10) & !is.null(oldD$Q10)) {
            #     grad.Q10=newData$Q10 - oldD$Q10
            # }else grad.Q10=NA
            # if(!is.null(newData$Q90) & !is.null(oldD$Q90)) {
            #     grad.Q90=newData$Q90 - oldD$Q90
            # }else grad.Q90=NA
        } else {
            gradVar=NA
            gradAmp=NA
            # grad.mos.aro_mos.arp=NA
            # grad.mos.aro_mos.cep=NA
            # grad.mos.arp_mos.cep=NA
            # grad.Q50_mos.aro=NA
            # grad.Q50_mos.arp=NA
            # grad.Q50_mos.cep=NA
            # grad.mos.aro=NA
            # grad.mos.arp=NA
            # grad.mos.cep=NA
            # grad.Q10=NA
            # grad.Q90=NA
        }
        newData=cbind(newData,gradVar,gradAmp
                        # ,grad.mos.aro_mos.arp,grad.mos.aro_mos.cep,grad.mos.arp_mos.cep,grad.mos.aro,grad.mos.arp,grad.mos.cep
                        # ,grad.Q50_mos.aro,grad.Q50_mos.arp,grad.Q50_mos.cep
                        # ,grad.Q10,grad.Q90
                        )
        newGradNames=paste("grad",i,".",c("Var","Amp"
                        # ,"mos.aro_mos.arp","mos.aro_mos.cep","mos.arp_mos.cep","mos.aro","mos.arp","mos.cep"
                        # ,"Q50_mos.aro","Q50_mos.arp","Q50_mos.cep"
                        # ,"Q10","Q90"
                        ),sep="")
        #grad1XXX grad of XXX between iteration t and iteration t-1
        #grad2XXX grad of XXX between iteration t and iteration t-2
        #grad3XXX grad of XXX between iteration t and iteration t-3
        names(newData)=c(oldNames,newGradNames)
    }
    return(newData)
}

#compute the distance (depending on normType) between the training data and the new data
dist_pred_train <- function(training_data,newData,normType) {
    nl=nrow(training_data)
    nc=length(newData)
    duplicateData=matrix(newData,nrow=nl,ncol=ncol(newData),byrow=T)#we duplicate the newData in a matrix
    train_data=as.matrix(training_data[,1:nc]) #[,1:nc] : only the predictors
    diff= train_data - duplicateData
    if (normType=="infinity") {
        abs_diff=abs(diff)
        max_abs_diff=apply(abs_diff,1,max)
        distance=min(max_abs_diff,na.rm=T)
    } else if (normType=="one") {
        abs_diff=abs(diff)
        mean_abs_diff=apply(abs_diff,1,mean)
        distance=min(mean_abs_diff,na.rm=T)
    } else if (normType=="two") {
        diffdiff=diff*diff
        mean_diffdiff=apply(diffdiff,1,mean)
        sqrt_mean_diffdiff=sqrt(mean_diffdiff)
        distance=min(sqrt_mean_diffdiff,na.rm=T)
    } else stop("problem with normType")
    return(distance)
}


#get the predictions of the lead time ech and the corresponding observation
getFirstLeadTime <- function(sta,res,date.valid,ech) {
    file_name=dir(path = dir_data_with, pattern = paste("donnee_t_",res,"_",ech,"_",sta,sep=""),full.names=TRUE)
    #il y a un ou des espaces avant les NA..., d'ou "  NA"
    data = read.table(file_name,head=T,sep=";",na.strings=c("       NA","      NA","     NA","    NA","   NA","  NA"," NA","NA"))
    dateObs=which(substr(data$valid,1,10)==substr(date.valid,1,10))
    if (length(dateObs)!=1) stop("dans la fonction getFirstLeadTime, il y a un problème avec la date
                                    à laquelle on veut récupérer l'obs, est-elle au bon format ?")
    data = data[dateObs,!names(data) %in% c("ech","res","insee","valid")]
}

#compute the variance of the variance during a run
computeRunVar <- function(sta,res,iteration,names.experts){
    names.experts[which(names.experts=="Q30")]="Q20" # TODO verrue !!!!!! it's the Q30 and not the Q20, it would be better to change ths in the code which creates the .txt
    varT=c()
    listeEch=c("3","6","9","12","15","18","21","24","27","30","33","36","39","42","45","48","57","72","84") #TODO pas beau en dur
    for (ech in listeEch) {
        file_name=dir(path = dir_data_with, pattern = paste("donnee_t_",res,"_",ech,"_",sta,sep=""),full.names=TRUE)
        #il y a un ou des espaces avant les NA..., d'ou "  NA"
        data = read.table(file_name,head=T,sep=";",na.strings=c("       NA","      NA","     NA","    NA","   NA","  NA"," NA","NA"))
        varT = c(varT,apply(data[iteration,names(data) %in% c(names.experts)],1,var,na.rm=T))
    }
    runVarVar=var(varT,na.rm=T)
    runMeanVar=mean(varT,na.rm=T)
    return(list(runMeanVar=runMeanVar,runVarVar=runVarVar))
}
