prevision_MLpol <-function(X,Y,window,T,sta,ech){#fonction qui fait la prevision par MLpol avec la fenêtre glissante window, retourne la liste des previsions
    liste_predictions=c()#vecteur qui contiendra l'ensemble des predictions de l'agregation
    liste_eta=array(0,dim=c(0,ncol(X))) #array of the eta for each expert, each line is one iteration
    colnames(liste_eta) <- colnames(X)
    matrix_weights=matrix(1/ncol(X),ncol=ncol(X),nrow=0)#vecteur qui contiendra l'ensemble des poids de l'agregation, pas de ligne w0 pour avoir meme longueur que Y, sinon pb de plot
    data_MLpol = mixture(model="MLpol",loss.type="square",loss.gradient=TRUE)
    data_MLpol$training$sta=sta
    data_MLpol$training$ech=ech

    for (t in seq(1,T)){
        data_MLpol$training$iteration=t
        if(data_MLpol$T>(window -1)){
            data_MLpol$training$R=data_MLpol$training$R-data_MLpol$training$r[1,]
            data_MLpol$parameters$eta=data_MLpol$parameters$eta[2:window,]  	
            data_MLpol$prediction=matrix(data_MLpol$prediction[2:window], ncol=1, nrow=(window -1))
            data_MLpol$Y=matrix(data_MLpol$Y[2:window],ncol=1,nrow=(window -1))
            data_MLpol$experts=data_MLpol$experts[2:window,]
            data_MLpol$awake=data_MLpol$awake[2:window,]
            data_MLpol$training$r=data_MLpol$training$r[2:window,]
            data_MLpol$training$B=max(data_MLpol$training$r^2)
            data_MLpol$T=window-1
            data_MLpol$loss=mean(loss(c(data_MLpol$prediction), c(data_MLpol$Y), loss.type = data_MLpol$loss.type))
        }
        data_MLpol = predict(data_MLpol,newexpert=X[t,],newY=Y[t],online=TRUE)
        liste_predictions=c(liste_predictions,tail(data_MLpol$pred,1)) #on ajoute la dernière prediction de l'agregation
        liste_eta=rbind(liste_eta,data_MLpol$parameters$eta[nrow(data_MLpol$parameters$eta),])
    }
    MLpol=mixture(model="MLpol",loss.type="square") #objet mixture avec les poids qui auraient permis de faire les prevision de liste_predicitons
    MLpol$weights=as.matrix(data_MLpol$weights)
    MLpol$prediction=c(liste_predictions)
    MLpol$Y=Y
    MLpol$experts=X
    MLpol$names.experts=colnames(X)
    MLpol$T=T
    MLpol$d=data_MLpol$d
    MLpol$loss=mean(loss(c(liste_predictions), c(Y), loss.type = data_MLpol$loss.type))
    MLpol$awake_transition=data_MLpol$awake_transition
    MLpol$awake=matrix(1,nrow=T,ncol=ncol(X))
    idx.na <- which(is.na(X)) #missing values of the experts
    MLpol$awake[idx.na] <- 0 #like in predictReal of Opera
    MLpol$experts<- as.matrix(MLpol$experts)
    MLpol$experts[idx.na]=0 #like in predictReal of Opera
    print(paste("MLpol window : ", window, "OK"))
    return(list(pred=liste_predictions,weights=MLpol$weights,mixture=MLpol,liste_eta=liste_eta))
}