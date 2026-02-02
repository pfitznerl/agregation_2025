prevision_EWA <-function(X,Y,window,T,sta,ech){#fonction qui fait la prevision par EWA avec la fenêtre glissante window, retourne la liste des previsions
    liste_predictions=c()#vecteur qui contiendra l'ensemble des predictions de l'agregation
    matrix_weights=matrix(1/ncol(X),ncol=ncol(X),nrow=0)#vecteur qui contiendra l'ensemble des poids de l'agregation, pas de ligne w0 pour avoir meme longueur que Y, sinon pb de plot
    data_EWA = mixture(model="EWA",loss.type="square",parameters=list(eta=eta),loss.gradient=TRUE)
    data_EWA$training$sta=sta
    data_EWA$training$ech=ech
    for (t in seq(1,T)){
        data_EWA$training$iteration=t
        data_EWA$training$window=window
        if(data_EWA$T>(window -1)){
            data_EWA$training$R=data_EWA$training$R-data_EWA$training$r[1,]
            data_EWA$training$cumulativeLoss <- data_EWA$training$cumulativeLoss -data_EWA$training$iLoss[1]
            data_EWA$prediction=matrix(data_EWA$prediction[2:window], ncol=1, nrow=(window -1))
            data_EWA$Y=matrix(data_EWA$Y[2:window],ncol=1,nrow=(window -1))
            data_EWA$experts=data_EWA$experts[2:window,]
            data_EWA$awake=data_EWA$awake[2:window,]
            data_EWA$loss=mean(loss(c(data_EWA$prediction), c(data_EWA$Y), loss.type = data_EWA$loss.type))
            data_EWA$training$r=data_EWA$training$r[2:window,]
            data_EWA$training$lpred=data_EWA$training$lpred[2:window]
            data_EWA$training$iLoss=data_EWA$training$iLoss[2:window]
            data_EWA$training$cumulativeLoss=data_EWA$training$cumulativeLoss[2:window]
            lmax=apply(data_EWA$training$r+data_EWA$training$lpred, 1, max) #max for each iteration
            lmin=apply(data_EWA$training$r+data_EWA$training$lpred, 1, min) #min for eac iteration
            data_EWA$training$lossMax=max(abs(lmax-lmin),na.rm=TRUE)
            data_EWA$training$V=data_EWA$training$V-data_EWA$training$v[1]
            data_EWA$training$v=data_EWA$training$v[2:window]
            data_EWA$T=window-1
        }
        data_EWA = predict(data_EWA,newexpert=X[t,],newY=Y[t],online=TRUE)
        liste_predictions=c(liste_predictions,tail(data_EWA$pred,1)) #on ajoute la dernière prediction de l'agregation
    }
    EWA = mixture(model="EWA",loss.type="square") #objet mixture avec les poids qui auraient permis de faire les prevision de liste_predicitons
    EWA$weights=as.matrix(data_EWA$weights)
    EWA$prediction=c(liste_predictions)
    EWA$Y=Y
    EWA$experts=X
    EWA$names.experts=colnames(X)
    EWA$T=T
    EWA$d=data_EWA$d
    EWA$loss=mean(loss(c(liste_predictions), c(Y), loss.type = data_EWA$loss.type))
    EWA$awake=matrix(1,nrow=T,ncol=ncol(X))
    idx.na <- which(is.na(X)) #missing values of the experts
    EWA$awake[idx.na] <- 0 #like in predictReal of Opera
    EWA$experts<- as.matrix(EWA$experts)
    EWA$experts[idx.na] <- 0 #like in predictReal of Opera
    print(paste("EWA window : ", window, "OK"))
    return(list(pred=liste_predictions,weights=EWA$weights,mixture=EWA))
}


