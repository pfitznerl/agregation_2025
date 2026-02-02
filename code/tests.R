# To prevent running on a login node !!!!
hostName = Sys.info()[4]
if (substr(hostName, 1, 12) == "belenoslogin") stop("!!!!!!!!!!!!!!!! You are on a login node !!!!!!!!!!!!!!!!!!!!!")
if (substr(hostName, 1, 12) == "taranislogin") stop("!!!!!!!!!!!!!!!! You are on a login node !!!!!!!!!!!!!!!!!!!!!")

graphics.off()
setwd("./.")

library(forecast)
library(foreach)
library(EnvStats)
suppressMessages(library(doMC))
library(parallel)
library(dplyr)

registerDoMC(100)

dir_out="sorties/"
dir_data_with="data/donnee_avec_pearp"
parametre="t"
res="00"
liste_ech=c("6","9","12","15","18","21","24","27","30","33","36","39","42","45","48","57","72","84") #without lead time 3 because we use the observation of lead time 3 in the random forests

liste_exp_without=c("sd.aro","as.aro","sd.arp","as.arp","sd.cep","as.cep")
liste_exp_with=c("raw.aro","mos.aro","raw.arp","mos.arp","raw.cep","mos.cep","Q10","Q30","Q50","Q70","Q90") #experts to use if available in the used data

liste_sta=dir(path = dir_data_with, pattern = paste("donnee_",parametre,"_00_48_",sep=""),full.names=TRUE)

liste_window=c(1253)

liste_agreg<- c("bestConvex","BOA","EWA","MLpol","MLprod","NBM","UNIFORM","bestExpert","raw.arp","mos.arp","raw.cep","mos.cep")
liste_agreg<- c("UNIFORM","bestExpert")

E1=c()
E2=c()


# Array of dimensions window*exp*exp*ech*sta
DM <- array(NA, dim = c(length(liste_window),length(liste_agreg),length(liste_agreg),length(liste_ech),length(liste_sta)))
QT <- array(NA, dim = c(length(liste_window),length(liste_agreg),length(liste_agreg),length(liste_ech),length(liste_sta)))
for (iwindow in 1:length(liste_window)) {
    window=liste_window[iwindow]
    for (ista in 1:length(liste_sta)) {
        sta=liste_sta[ista]
        sta=strsplit(sta, "_")[[1]][7]
        sta=strsplit(sta, ".txt")[[1]][1]
        for(iech in 1:length(liste_ech)) {
            ech=liste_ech[iech]
            for (iagreg1 in 1:length(liste_agreg)){
                agreg1=liste_agreg[iagreg1]
                AGREG1=agreg1
                if(AGREG1=="raw.arp"|AGREG1=="mos.arp"|AGREG1=="raw.cep"|AGREG1=="mos.cep"){#|AGREG1=="UNIFORM") {
                    AGREG1="BOA"
                }
                data1=read.table(paste(dir_out,AGREG1,"_",sta,"_",window,"_",ech,".txt",sep=""),header=TRUE,sep=";")
                icol1=which(colnames(data1)==agreg1)
                e1=data1[,icol1]-data1$Y
                if(agreg1=="UNIFORM")E1=c(E1,e1)
                if(agreg1=="bestExpert")E2=c(E2,e1)
                idNA1=which(!is.na(e1))
                if(length(icol1)>0) e1=e1[idNA1] else e1=rep(0,1253)
                for (iagreg2 in 1:length(liste_agreg)) {
                    agreg2=liste_agreg[iagreg2]
                    AGREG2=agreg2
                    if(AGREG2=="raw.arp"|AGREG2=="mos.arp"|AGREG2=="raw.cep"|AGREG2=="mos.cep"){#|AGREG2=="UNIFORM") {
                        AGREG2="BOA"
                    }
                    print(paste(agreg1," vs ",agreg2,"      ech : ",ech,"    sta :",sta))
                    data2=read.table(paste(dir_out,AGREG2,"_",sta,"_",window,"_",ech,".txt",sep=""),header=TRUE,sep=";")
                    icol2=which(colnames(data2)==agreg2)
                    e2=data2[,icol2]-data2$Y
                    if(length(icol1)>0 & length(icol2)>0) e2=e2[idNA1] else e2=rep(0,1253)#on supprime les lignes correspondant aux NA de e1, même si on devrait avoir idNA1=idNA2 normalement
                    idNA2=which(!is.na(e2))
                    if(length(icol1)>0 & length(icol2)>0) e2=e2[idNA2] else e2=rep(0,1253)
                    if(length(icol2)>0 & length(icol2)>0) e1=e1[idNA2] else e1=rep(0,1253)

                    if (iagreg1!=iagreg2 & sum(e1-e2,na.rm=T)!=0) {
                        dm=dm.test(e1^2,e2^2,h=ech,alternative="less",varestimator = "acf",power=2)$p.value
                        data_qt=as.data.frame(rbind(cbind(e1,agreg1),cbind(e2,agreg2)))
                        data_qt[,1]=as.numeric(data_qt[,1])
                        colnames(data_qt)=c("err","agreg")
                        qt=quantileTest(abs(e1),abs(e2),alternative="less",target.quantile=0.95)$p.value
                    } else if (sum(e1-e2,na.rm=T)==0) {
                       dm=1
                       qt=1
                    } else {
                        dm=NA
                        qt=NA
                    }
                    DM[iwindow,iagreg1,iagreg2,iech,ista]=dm
                    QT[iwindow,iagreg1,iagreg2,iech,ista]=qt
                }
            }
        }
    }
}

# Array of dimensions window*exp*exp
BH_DM <- array(NA, dim = c(length(liste_window),length(liste_agreg),length(liste_agreg)),
            dimnames = list(liste_window,liste_agreg,liste_agreg))

for (iwindow in 1:length(liste_window)) {
    window=liste_window[iwindow]
    for (iagreg1 in 1:length(liste_agreg)){
        for (iagreg2 in 1:length(liste_agreg)) {
            p_vals=c(DM[iwindow,iagreg1,iagreg2,,])
            sta_x_ech=length(p_vals)
            bh=p.adjust(p_vals, method = 'BH', n = sta_x_ech)
            bh=length(which(bh <= 0.05))/sta_x_ech
            percentage_bh=100*bh/length(bh) #percentage of significative DM tests after BH
            BH_DM[iwindow,iagreg1,iagreg2]=percentage_bh
        }
    }
}

# Array of dimensions window*exp*exp
BH_QT <- array(NA, dim = c(length(liste_window),length(liste_agreg),length(liste_agreg)),
            dimnames = list(liste_window,liste_agreg,liste_agreg))

for (iwindow in 1:length(liste_window)) {
    window=liste_window[iwindow]
    for (iagreg1 in 1:length(liste_agreg)){
        for (iagreg2 in 1:length(liste_agreg)) {
            p_vals=c(QT[iwindow,iagreg1,iagreg2,,])
            sta_x_ech=length(p_vals)
            bh=p.adjust(p_vals, method = 'BH', n = sta_x_ech)
            bh=length(which(bh <= 0.05))/sta_x_ech
            percentage_bh=100*bh/length(bh) #percentage of significative DM tests after BH
            BH_QT[iwindow,iagreg1,iagreg2]=percentage_bh
        }
    }
}

nameRdata="out.Rdata"
save.image(nameRdata)


