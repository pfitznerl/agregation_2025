prevision_MLprod <-function(X,Y,window,T,sta,ech){#fonction qui fait la prevision par MLprod avec la fenêtre glissante window, retourne la liste des previsions
    liste_predictions=c()#vecteur qui contiendra l'ensemble des predictions de l'agregation
    liste_eta=array(0,dim=c(0,ncol(X))) #array of the eta for each expert, each line is one iteration
    colnames(liste_eta) <- colnames(X)
    matrix_weights=matrix(1/ncol(X),ncol=ncol(X),nrow=0)#vecteur qui contiendra l'ensemble des poids de l'agregation, pas de ligne w0 pour avoir meme longueur que Y, sinon pb de plot
    data_MLprod = mixture(model="MLprod",loss.type="square",loss.gradient=TRUE)
    data_MLprod$training$sta=sta
    data_MLprod$training$ech=ech
    for (t in seq(1,T)){
        data_MLprod$training$iteration=t
        if(data_MLprod$T>(window -1)){
            data_MLprod$training$R=data_MLprod$training$R- log(1 + data_MLprod$awake[window,] * data_MLprod$parameters$eta[window,] * data_MLprod$training$r[1,])
            data_MLprod$parameters$eta=data_MLprod$parameters$eta[2:window,]  	
            data_MLprod$prediction=matrix(data_MLprod$prediction[2:window], ncol=1, nrow=(window -1))
            data_MLprod$Y=matrix(data_MLprod$Y[2:window],ncol=1,nrow=(window -1))
            data_MLprod$experts=data_MLprod$experts[2:window,]
            data_MLprod$awake=data_MLprod$awake[2:window,]
            data_MLprod$training$r=data_MLprod$training$r[2:window,]
            data_MLprod$training$L=apply(data_MLprod$training$r^2,2,max)
            data_MLprod$training$maxLoss=apply(data_MLprod$training$r,2,max)
            data_MLprod$T=window-1
            data_MLprod$loss=mean(loss(c(data_MLprod$prediction), c(data_MLprod$Y), loss.type = data_MLprod$loss.type))
        }
        data_MLprod = predict(data_MLprod,newexpert=X[t,],newY=Y[t],online=TRUE)
        liste_predictions=c(liste_predictions,tail(data_MLprod$pred,1)) #on ajoute la dernière prediction de l'agregation
        matrix_weights=rbind(matrix_weights,c(data_MLprod$coefficients)) #on ajoute les derniers poids de l'agregation
        liste_eta=rbind(liste_eta,data_MLprod$parameters$eta[nrow(data_MLprod$parameters$eta),])
    }
    MLprod=mixture(model="MLprod",loss.type="square") #objet mixture avec les poids qui auraient permis de faire les prevision de liste_predicitons
    MLprod$weights=as.matrix(data_MLprod$weights)
    MLprod$prediction=c(liste_predictions)
    MLprod$Y=Y
    MLprod$experts=X
    MLprod$names.experts=colnames(X)
    MLprod$T=T
    MLprod$d=data_MLprod$d
    MLprod$loss=mean(loss(c(liste_predictions), c(Y), loss.type = data_MLprod$loss.type))
    MLprod$awake_transition=data_MLprod$awake_transition
    MLprod$awake=matrix(1,nrow=T,ncol=ncol(X))
    idx.na <- which(is.na(X)) #missing values of the experts
    MLprod$awake[idx.na] <- 0 #like in predictReal of Opera
    MLprod$experts<- as.matrix(MLprod$experts)
    MLprod$experts[idx.na]=0 #like in predictReal of Opera
    print(paste("MLprod window : ", window, "OK"))
    return(list(pred=liste_predictions,weights=MLprod$weights,mixture=MLprod,liste_eta=liste_eta))
}