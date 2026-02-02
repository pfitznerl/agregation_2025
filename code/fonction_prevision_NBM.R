prevision_NBM <-function(X,Y,window,T,sta,ech){#fonction qui fait la prevision par NBM avec la fenêtre glissante window, retourne la liste des previsions
    liste_predictions=c()#vecteur qui contiendra l'ensemble des predictions de l'agregation
    matrix_weights=matrix(1/ncol(X),ncol=ncol(X),nrow=0)#vecteur qui contiendra l'ensemble des poids de l'agregation, pas de ligne w0 pour avoir meme longueur que Y, sinon pb de plot
    data_NBM = mixture(model="NBM",loss.type="square",parameters=list(eta=eta),loss.gradient=TRUE)
    data_NBM$training$sta=sta
    data_NBM$training$ech=ech
    file_pour_date=dir(path = dir_data_with, pattern = paste("donnee_",parametre,"_",res,"_",ech,"_",sta,sep=""),full.names=TRUE)#verrue, il vaudrait mieux que ce soit en argument de prevision_NBM
    data_pour_date=read.table(file_pour_date,head=T,sep=";",na.strings=c("       NA","      NA","     NA","    NA","   NA","  NA"," NA","NA")) #il y a un ou des espaces avant les NA..., d'ou "  NA"
    date_valid=data_pour_date$valid
    for (t in seq(1,T)){
        data_NBM$training$iteration=t
        data_NBM$training$window=window
        data_NBM$training$date_valid=date_valid[t]
        if(data_NBM$T>(window -1)){
            new = mixture(model="NBM",loss.type="square",parameters=list(eta=eta),loss.gradient=TRUE)
            new$training$sta=sta
            new$training$ech=ech
            data_NBM$awake=data_NBM$awake[2:window,]
            data_NBM$b=data_NBM$b[2:window,]
            data_NBM$mse=data_NBM$MSE[2:window,]
            data_NBM$T=window-1
            data_NBM$prediction=matrix(data_NBM$prediction[2:window], ncol=1, nrow=(window -1))
            data_NBM$Y=matrix(data_NBM$Y[2:window],ncol=1,nrow=(window -1))
            data_NBM$experts=data_NBM$experts[2:window,]
            data_NBM$loss=mean(loss(c(data_NBM$prediction), c(data_NBM$Y), loss.type = data_NBM$loss.type))
            data_NBM$training$mse <- data_NBM$training$mse[2:window,]
            data_NBM$T=window-1
        }
        data_NBM = predict(data_NBM,newexpert=X[t,],newY=Y[t],online=TRUE)
        liste_predictions=c(liste_predictions,tail(data_NBM$pred,1)) #on ajoute la dernière prediction de l'agregation
    }
    NBM = mixture(model="NBM",loss.type="square") #objet mixture avec les poids qui auraient permis de faire les prevision de liste_predicitons
    NBM$weights=as.matrix(data_NBM$weights) #as.matrix(matrix_weights)
    NBM$prediction=c(liste_predictions)
    NBM$Y=Y
    NBM$experts=X
    NBM$names.experts=colnames(X)
    NBM$T=T
    NBM$d=data_NBM$d
    NBM$loss=mean(loss(c(liste_predictions), c(Y), loss.type = data_NBM$loss.type))
    NBM$awake=matrix(1,nrow=T,ncol=ncol(X))
    idx.na <- which(is.na(X)) #missing values of the experts
    NBM$awake[idx.na] <- 0 #like in predictReal of Opera
    NBM$experts<- as.matrix(NBM$experts)
    NBM$experts[idx.na] <- 0 #like in predictReal of Opera
    print(paste("NBM window : ", window, "OK"))
    return(list(pred=liste_predictions,weights=NBM$weights,mixture=NBM))
}


