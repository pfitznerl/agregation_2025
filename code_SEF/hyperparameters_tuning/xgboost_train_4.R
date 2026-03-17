######################################
######################################
######################################
# tune les hyperparamètres de xgboost sur les jeux de donnés (ici, mélange de toutes les stations et échéances)
# créés par creation_donnees.R dans data
######################################
######################################
######################################

print("start date")
print(system("date"))
cat("\n")
rm(list=ls())
library(xgboost)
# library(caret)
# library(randomForestSRC)
library(plyr)
library(dplyr)
library(mlr3verse)
library(data.table)
source("transition_learning_xgboost.R")

liste_data=dir(path="data_creation", pattern = paste("RData",sep=""),full.names=TRUE)

liste_sta=substring(liste_data, 51,58)
liste_sta=unique(liste_sta)

parameters_list = list()


set.seed(777)  # Pour la 'reproductibilité


iparam=0



liste_nrounds=c(1,2,3,4,5,6,7)
liste_max_depth=c(4,5,6,7,8,9,10)
liste_eta=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
liste_colsample_bynode=c(0.65,0.7,0.75,0.8,0.85,0.9)
liste_subsample=c(1)
liste_gamma=c(0)
liste_min_child_weight=c(5,10,20,35,60,100)


print(paste("liste_nrounds : ",liste_nrounds,sep=""))
print(paste("liste_max_depth : ",liste_max_depth,sep=""))
print(paste("liste_eta : ",liste_eta,sep=""))
print(paste("liste_colsample_bynode : ",liste_colsample_bynode,sep=""))
print(paste("liste_subsample : ",liste_subsample,sep=""))
print(paste("liste_min_child_weight : ",liste_min_child_weight,sep=""))
print(paste("liste_gamma : ",liste_gamma,sep=""))
for (inrounds in liste_nrounds){
    for (imax_depth in liste_max_depth){
        for (ieta in liste_eta){
            for (igamma in liste_gamma){
                for (icolsample_bynode in liste_colsample_bynode) {
                    for(isubsample in liste_subsample){
                        for (imin_child_weight in liste_min_child_weight){
                            iparam=iparam+1
                            param <- list(booster = "gbtree",
                                    objective = "reg:squarederror",
                                    nrounds=inrounds,
                                    max_depth = imax_depth,
                                    eta = ieta,
                                    gamma = igamma,
                                    subsample = isubsample,
                                    colsample_bynode = icolsample_bynode,
                                    min_child_weight = imin_child_weight
                            )
                            parameters <- as.data.frame(param)
                            parameters_list[[iparam]] <- parameters
                        }
                    }
                }
            }
        }
    }
}

# Create object that contains all randomly created hyperparameters
parameters_df = do.call(rbind, parameters_list)

print(paste("il y a ",nrow(parameters_df)," jeux d'hyperparamètres",sep=""))

PSS_param=data.frame()

nExp=0
args <- commandArgs(TRUE)
nExp=as.numeric(args[1])


#training
for (i in nExp){
    print(paste("training ",i,sep=""))

    #######################################
    #######################################
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
    options(na.action='na.pass') #to allow NA in sparse matrices
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
    #######################################
    #######################################

    ficData=paste("data/data_sample_BOA_",i,".RData",sep="")
    load(ficData) #load data : training and validation

    threshold=3.5
    print(paste("threshold : ",threshold,sep=""))

    #define predictor and response variables in training set
    training_matrix <-model.matrix(err ~.-1, data = training)

    validation_matrix <-model.matrix(err ~.-1, data = validation)

    dtrain <- xgb.DMatrix(data = training_matrix, label = training$err)
    dvalid <- xgb.DMatrix(data = validation_matrix, label = validation$err)

    for (row in 1:nrow(parameters_df)){
        if(row==1) print(paste(i,"ème jeux de training",sep=""))
        if(row%%250==0) print(paste(row,"ème jeux de paramètres",sep=""))
        xgb_model <- xgb.train(data=dtrain,
                            booster = "gbtree",
                            objective = "reg:squarederror",
                            max_depth = parameters_df$max_depth[row],
                            eta = parameters_df$eta[row],
                            subsample = parameters_df$subsample[row],
                            colsample_bynode = parameters_df$colsample_bynode[row],
                            min_child_weight = parameters_df$min_child_weight[row],
                            gamma = parameters_df$gamma[row],
                            nrounds= parameters_df$nrounds[row],
                            eval_metric = "rmse",
                            verbose=0,
                            print_every_n = 10
        )

        predicted = predict(xgb_model, dvalid)

        #skill for e>=threshold°C and e<threshold°C
        A=length(which(predicted>=threshold & validation$err >=threshold))
        B=length(which(predicted>=threshold & validation$err <threshold))
        C=length(which(predicted<threshold & validation$err >=threshold))
        D=length(which(predicted<threshold & validation$err <threshold))

        #skill for e<=-threshold°C and e>threshold°C
        A2=length(which(predicted<=-threshold & validation$err <=-threshold))
        B2=length(which(predicted<=-threshold & validation$err >-threshold))
        C2=length(which(predicted>-threshold & validation$err <=-threshold))
        D2=length(which(predicted>-threshold & validation$err >-threshold))
        newData=cbind(i,row,parameters_df[row,],A,B,C,D,A2,B2,C2,D2)
        colnames(newData)[1]="training"
        colnames(newData)[2]="nparam"
        PSS_param=rbind(PSS_param,newData)
    }

}

# on sauvegarde les données
# PSS_param contains on each line, ntraining then number of the set of hyperparameters,
# then the hyperparameters and then the numbers A,B,C,D needed to compute the PSS
print(paste("we save : results/results_",nExp,".RData",sep=""))
nameR=paste("results/results_",nExp,".RData",sep="")
save(PSS_param,file=nameR)



