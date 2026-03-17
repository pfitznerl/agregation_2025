rm(list=ls())
library(xgboost)
# library(caret)
# library(randomForestSRC)
library(plyr)
library(dplyr)
library(mlr3verse)
library(data.table)

PSS_param_tot=data.frame()



# PSS_param contains on each line, ntraining then number of the set of hyperparameters,
# then the hyperparameters and then the numbers A,B,C,D needed to compute the PSS and then the RMSE
for (i in 1:100){
    ficData=paste("results/results_",i,".RData",sep="")
    print(ficData)
    load(ficData)
    PSS_param_tot=rbind(PSS_param_tot,PSS_param)
}


tot_PSS=data.frame()
for (iparam in 1:max(PSS_param_tot$nparam)){
    data=PSS_param_tot[which(PSS_param_tot$nparam==iparam),]#on selectionne uniquement les lignes des hyperparametres iparam
    tot_A=sum(data$A)
    tot_B=sum(data$B)
    tot_C=sum(data$C)
    tot_D=sum(data$D)
    tot_A2=sum(data$A2)
    tot_B2=sum(data$B2)
    tot_C2=sum(data$C2)
    tot_D2=sum(data$D2)
    mean_rmse=mean(data$RMSE)
    param=subset(data[1,], select = c(nrounds,max_depth,eta,colsample_bynode,subsample,gamma,min_child_weight))
    new=cbind(param,tot_A,tot_B,tot_C,tot_D,tot_A2,tot_B2,tot_C2,tot_D2)#hyperparamètres ipram, puis somme des A,B,C,D associés aux hyperparamètres iparam sur tous les jeux d'entrainements
    tot_PSS=rbind(tot_PSS,new)
}

#cf confusion matrix wikipedia
tot_PSS$hit_rate=tot_PSS$tot_A/(tot_PSS$tot_A+tot_PSS$tot_C)
tot_PSS$false_rate=tot_PSS$tot_B/(tot_PSS$tot_B+tot_PSS$tot_D)
tot_PSS$false_ratio=tot_PSS$tot_B/(tot_PSS$tot_A+tot_PSS$tot_B)
tot_PSS$true_ratio=tot_PSS$tot_A/(tot_PSS$tot_A+tot_PSS$tot_B)
tot_PSS$PSS=tot_PSS$hit_rate-tot_PSS$false_rate

tot_PSS$hit_rate2=tot_PSS$tot_A2/(tot_PSS$tot_A2+tot_PSS$tot_C2)
tot_PSS$false_rate2=tot_PSS$tot_B2/(tot_PSS$tot_B2+tot_PSS$tot_D2)
tot_PSS$false_ratio2=tot_PSS$tot_B2/(tot_PSS$tot_A2+tot_PSS$tot_B2)
tot_PSS$true_ratio2=tot_PSS$tot_A2/(tot_PSS$tot_A2+tot_PSS$tot_B2)
tot_PSS$PSS2=tot_PSS$hit_rate2-tot_PSS$false_rate2

mean_PSS=(tot_PSS$PSS+tot_PSS$PSS2)/2

cat("\n")
print("PSS max : ")
print(max(mean_PSS,na.rm=TRUE))
print("hyperparameters with the best PSS :")
print(tot_PSS[which(mean_PSS==max(mean_PSS,na.rm=TRUE)),])
cat("\n")


print("ne pas regarder rmse, c'est faux")


