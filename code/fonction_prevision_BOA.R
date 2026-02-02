prevision_BOA <-function(X,Y,window,T,sta,ech,date.valid){#fonction qui fait la prevision par BOA avec la fenêtre glissante window, retourne la liste des previsions
    liste_predictions=c()#vecteur qui contiendra l'ensemble des predictions de l'agregation
    liste_eta=array(0,dim=c(0,ncol(X))) #array of the eta for each expert, each line is one iteration
    colnames(liste_eta) <- colnames(X)
    matrix_weights=matrix(1/ncol(X),ncol=ncol(X),nrow=0)#vecteur qui contiendra l'ensemble des poids de l'agregation, pas de ligne w0 pour avoir meme longueur que Y, sinon pb de plot
    data_BOA = mixture(model="BOA",loss.type="square",loss.gradient=TRUE)
    data_BOA$training$sta=sta
    data_BOA$training$ech=ech
    for (t in seq(1,T)){
        data_BOA$training$date.valid=date.valid[t]
        data_BOA$training$iteration=t
        if(data_BOA$T>(window -1)){
            data_BOA$training$R=data_BOA$training$R-data_BOA$training$r[1,]
            data_BOA$training$R.reg=data_BOA$training$R.reg-data_BOA$training$r.reg[1,]
            data_BOA$training$L=data_BOA$training$L-data_BOA$training$l
            data_BOA$training$ref_L=data_BOA$training$ref_L-data_BOA$training$ref_l
            data_BOA$training$all_L=data_BOA$training$all_L-data_BOA$training$all_l
            data_BOA$training$V=data_BOA$training$V-2.2*data_BOA$training$r[1,]^2
            data_BOA$training$V[which(data_BOA$training$V<0)]=0
            data_BOA$parameters$eta=data_BOA$parameters$eta[2:window,]  	
            data_BOA$prediction=matrix(data_BOA$prediction[2:window], ncol=1, nrow=(window -1))
            data_BOA$Y=matrix(data_BOA$Y[2:window],ncol=1,nrow=(window -1))
            data_BOA$experts=data_BOA$experts[2:window,]
            data_BOA$awake=data_BOA$awake[2:window,]
            data_BOA$training$r=data_BOA$training$r[2:window,]
            data_BOA$training$r.reg=data_BOA$training$r.reg[2:window,]

            data_BOA$T=window-1
            data_BOA$loss=mean(loss(c(data_BOA$prediction), c(data_BOA$Y), loss.type = data_BOA$loss.type))
        }

        data_BOA = predict(data_BOA,newexpert=X[t,],newY=Y[t],online=TRUE)
        liste_predictions=c(liste_predictions,tail(data_BOA$pred,1)) #on ajoute la dernière prediction de l'agregation
        liste_eta=rbind(liste_eta,data_BOA$parameters$eta[nrow(data_BOA$parameters$eta),])
    }
    BOA = mixture(model="BOA",loss.type="square") #objet mixture avec les poids qui auraient permis de faire les prevision de liste_predicitons
    BOA$weights=as.matrix(data_BOA$weights)
    BOA$prediction=c(liste_predictions)
    BOA$Y=Y
    BOA$experts=X
    BOA$names.experts=colnames(X)
    BOA$T=T
    BOA$d=data_BOA$d
    BOA$loss=mean(loss(c(liste_predictions), c(Y), loss.type = data_BOA$loss.type))
    BOA$awake_transition=data_BOA$awake_transition
    BOA$awake=matrix(1,nrow=T,ncol=ncol(X))
    idx.na <- which(is.na(X)) #missing values of the experts
    BOA$awake[idx.na] <- 0 #like in predictReal of Opera
    BOA$experts<- as.matrix(BOA$experts)
    BOA$experts[idx.na]=0 #like in predictReal of Opera
    NERR=data_BOA$training$NERR #number of times GBRT wakes up the wrong experts
    if(is.null(data_BOA$training$A1)) A1= NA else A1=data_BOA$training$A1
    if(is.null(data_BOA$training$B1)) B1= NA else B1=data_BOA$training$B1
    if(is.null(data_BOA$training$C1)) C1= NA else C1=data_BOA$training$C1
    if(is.null(data_BOA$training$D1)) D1= NA else D1=data_BOA$training$D1
    if(is.null(data_BOA$training$A2)) A2= NA else A2=data_BOA$training$A2
    if(is.null(data_BOA$training$B2)) B2= NA else B2=data_BOA$training$B2
    if(is.null(data_BOA$training$C2)) C2= NA else C2=data_BOA$training$C2
    if(is.null(data_BOA$training$D2)) D2= NA else D2=data_BOA$training$D2
    hit_rate1=A1/(A1+C1)
    false_rate1=B1/(B1+D1)
    PSS1=hit_rate1-false_rate1
    false_ratio1=B1/(A1+B1)
    hit_rate2=A2/(A2+C2)
    false_rate2=B2/(B2+D2)
    PSS2=hit_rate2-false_rate2
    false_ratio2=B2/(A2+B2)
    GSS=(PSS1+PSS2)/2
    print(paste("BOA window : ", window, "OK"))
    return(list(pred=liste_predictions,weights=BOA$weights,mixture=BOA,liste_eta=liste_eta,A1=A1,B1=B1,C1=C1,D1=D1,A2=A2,B2=B2,C2=C2,D2=D2,NERR=NERR))
}