library(maps)
library(ggrepel)

############################

#plot of a map of FRance with the stations on it
plot_map <- function(liste_station){
  coordinates=data.frame()
  data_sta=read.table("liste_station_obs_asmixtes",header=T)
  for (file_station in liste_station) {
    sta=substr(strsplit(file_station, "_")[[1]][7],1,8) #insee of the station
    if(substr(sta,1,1)==0) sta=substr(sta,2,substr(sta,2,nchar(sta))) #remove the first 0 if the insse begins with 0
    coordinates=rbind(coordinates,data_sta[which(data_sta$indic==sta),])
  }

  france <- map_data("france")
  plot.map <- ggplot()+
              geom_polygon(data = france, aes(x=long, y = lat, group = group),fill = "white",colour="grey") +
              #coord_fixed(1.3) +
              theme(axis.line = element_line(color='black'), #to remove the gridlines and put the background in white
                panel.background = element_rect(fill='transparent'),
                plot.background = element_rect(fill='transparent'),
                #plot.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.title.y = element_text(size=25),axis.title.x = element_text(size=25),
                axis.text.x = element_text(size = 25),axis.text.y = element_text(size = 25))+
              geom_point(data = coordinates, mapping = aes(x = lon, y = lat), size=5, color = "black") +
              geom_text_repel(data = coordinates, aes(x = lon, y = lat, label = indic), size = 7,max.overlaps=nrow(coordinates))
                # scaleBar(lon=-4.5, lat=42, distanceLon=200, distanceLat=20, distanceLegend=50,
                #           dist.unit = "km", orientation=FALSE)
                #  coord_sf(crs = 4326)
              
           
  ggsave(file=paste(dir_plot,"map_sta.pdf",sep=""), width=10, height=10, dpi=600)
  print(paste(dir_plot,"map_sta.pdf",sep=""))
}



#plot, boxplot and histogram of the parameters alpha eta, for FSBOA and FSCALIB
plot_a_e <- function(data,param,sta,window,start.date,end.date,names.experts){
  if(param=="alpha") {
    #data=log10(data)
    dir=dir_alpha
  }
  if (param=="eta"){
    #data= log2(data)
    dir=dir_eta
  }
  cnames=names.experts
  if (is.null(ncol(data))){#if only one parameter
    print(paste(param,"_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    pdf(paste(dir_plot,dir,param,"_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    plot(data)
    dev.off()
    print(paste(param,"_boxplot_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    pdf(paste(dir_plot,dir,param,"_boxplot_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    boxplot(data,outline=FALSE,main="boxplot of the learning rate, without the outliners")
    dev.off()
    print(paste(param,"_histogram_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    pdf(paste(dir_plot,dir,param,"_histogram_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    hist(data)
    dev.off()
  }
  else {#if multiple parameter (one for each expert for example)
    print(paste(param,"_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    pdf(paste(dir_plot,dir,param,"_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    matplot(data, type = c("b"),pch=1,col = 1:ncol(data),ylim=c(0,1))
    legend("topleft", legend = cnames, col=1:ncol(data), pch=1)
    dev.off()
    print(paste(param,"_boxplot_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    pdf(paste(dir_plot,dir,param,"_boxplot_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    boxplot(data,outline=FALSE,main="boxplot of the learning rate for each expert, without the outliners",col=1:ncol(data))
    legend("topleft", legend = cnames, col=1:ncol(data), pch=1)
    dev.off()
    print(paste(param,"_histogram_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    pdf(paste(dir_plot,dir,param,"_histogram_",agreg,"_window",window,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
    hist(data[,1],col=1,main="histogram")
    for (iexp in 2:ncol(data)){
      hist(data[,iexp],col=iexp,add=TRUE)
    }
    dev.off()
  }
}

################################
#plot of the rmse for the lead time ech
plot_rmse_ech <-function(liste_window,rmse.ech,rmse.lim){
  png(paste(dir_plot,dir_rmse,"rmse_resume","_ech",ech,"_tot.png",sep=""))
  plot(liste_window, unlist(rmse.ech[[1]]),ylim=rmse.lim,type="b",col=liste_couleur[1])
  legend("topright",legend=liste_agreg,fill=liste_couleur[1:length(liste_agreg)])
  abline(h=min(unlist(rmse.ech[[1]])), lty = "dotted",col=liste_couleur[1])
  #text(x=400,y=min(liste_recap[[1]]$rmse_mean),labels=min(liste_recap[[1]]$rmse_mean),col=liste_couleur[1])
  if (length(liste_recap)>1){#if more than one aggregation strategy
    for (iagreg in 2:length(liste_recap)){
      points(liste_window, unlist(rmse.ech[[iagreg]]),type="b",col=liste_couleur[iagreg])
      abline(h=min(unlist(rmse.ech[[iagreg]])), lty = "dotted",col=liste_couleur[iagreg])
      #text(x=400,y=min(liste_recap[[l]]$rmse_mean),labels=min(liste_recap[[l]]$rmse_mean),col=liste_couleur[l])
    }
  }
  dev.off()
}
################################

################################
#plot of the parameters for the lead time ech,
#(ech can be "AllEch", this means that all the parameters of all the lead times have been concatenated)
plot_parameters <-function(data_param,param,window,ech,n.sta){
  
  if(param=="alpha") {
    #data_param=log10(data_param)
    dir_param=dir_alpha
  }
  else if (param=="eta"){
    #data_param= log2(data_param)
    dir_param=dir_eta
  }

  print(paste(dir_plot,dir_param,param,"_resume_ech",ech,"_",window,"_tot.png",sep=""))
  png(paste(dir_plot,dir_param,param,"_resume_ech",ech,"_",window,"_tot.png",sep=""))
  if (is.null(ncol(data_param))){#if only one parameter
    plot(data_param, type = c("b"),pch=1,ylim=c(0,1))
  }
  else {#if multiple parameter (one for each expert for example)
    #only works if all the stations have the same number of iterations
    #!!!!!!!!!! 84087001 a une iteration de moins ? !!!!!!!!!!!!
    n.iteration=nrow(data_param)/n.sta
    data=matrix(nrow=n.iteration,ncol=0)
    if ( as.integer(n.iteration)!=n.iteration) warning("all the stations don't have the same number of iteration,
                                                      it's a problem for the next plot")
    for (i in 1:n.sta) {#we transform the columns in n*leadtimes columns in order to plot all the parameters of one iteration on the same axis
      start=1+(i-1)*n.iteration
      end=i*n.iteration
      data=cbind(data,data_param[start:end,])#data_param[start:end,] are all the parameters for one lead time
    }
    matplot(data, type = c("b"),pch=1,col = 1:ncol(data_param),ylim=c(0,1))
    legend("topleft", legend = colnames(data_param), col=1:ncol(data_param), pch=1)
  }
  dev.off()
}
################################

################################
#boxplot of the parameters for the lead time ech,
#(ech can be "AllEch", this means that all the parameters of all the lead times have been concatenated)
boxplot_parameters <-function(data_param,param,window,ech){
  if(param=="alpha") {
    #data_param=log10(data_param)
    dir_param=dir_alpha
  }
  else if (param=="eta"){
    #data_param= log2(data_param)
    dir_param=dir_eta
  }
  print(paste(dir_plot,dir_param,param,"_resume_boxplot_ech",ech,"_",window,"_tot.pdf",sep=""))
  pdf(paste(dir_plot,dir_param,param,"_resume_boxplot_ech",ech,"_",window,"_tot.pdf",sep=""))
  if (is.null(ncol(data_param))){#if only one parameter
    boxplot((data_param))
  }
  else {#if multiple parameter (one for each expert for example)
    boxplot(data_param,outline=FALSE,main="boxplot of the learning rate for each expert, without the outliners",col=1:ncol(data_param))
    legend("topleft", legend = colnames(data_param), col=1:ncol(data_param), pch=1)
  }
  dev.off()
}
################################



#plot des obs, des previsions des experts et de l'agregation en fonction du temps/des itérations de l'agregation
#pour une fenêtre glissante une echeance et une station
plot_obs_pred_window <-function(X,pred.agreg,Y,date.valid,window,sta,ech,agreg){
  X=as.data.frame(X)
  liste_expert=names(X)
  pred.agreg=as.data.frame(pred.agreg)
  df<-cbind(date.valid,pred.agreg,Y,X)
  if(length(which(substr(df$date.valid,1,10)==plot.periode[1]))==0 || length(which(substr(df$date.valid,1,10)==plot.periode[2]))==0)stop("plot period not in the data")
  df=df[which(substr(df$date.valid,1,10)==plot.periode[1]):which(substr(df$date.valid,1,10)==plot.periode[2]),] #we only keep the data between plot.periode[1] and plot.periode[2]
  names(df)[names(df) == "Y"] <- "obs"
  names(df)[names(df) == "pred.agreg" | names(df) == "V1"] <- agreg #V1 for bestConvex which has no name i don't know why
  names(df)[names(df) == "mos.aro" ] <- "ppm.aro"
  names(df)[names(df) == "mos.arp" ] <- "ppm.arp"
  names(df)[names(df) == "mos.cep" ] <- "ppm.ifs"
  names(df)[names(df) == "raw.cep" ] <- "raw.ifs"
  liste_expert[liste_expert == "mos.aro" ] <- "ppm.aro"
  liste_expert[liste_expert == "mos.arp" ] <- "ppm.arp"
  liste_expert[liste_expert == "mos.cep" ] <- "ppm.ifs"
  liste_expert[liste_expert == "raw.cep" ] <- "raw.ifs"
  df$date.valid=as.POSIXct(df$date.valid)
  df$date.valid=as.Date(df$date.valid) #type date from lubridate
  datebreaks <- seq(df$date.valid[1], df$date.valid[length(df$date.valid)], by="1 month")# by "1 month" or by "1 day"
  df<-df %>%
      dplyr::select(date.valid, c(all_of(liste_expert),obs)) %>%#,all_of(agreg)
      gather(key="modele", value="temperature",-date.valid)

  df <- transform(df, modele=factor(modele, levels=unique(modele))) #without this, ggplot don't plot the lines in the expected order (the order of the dataframe's columns)
  #liste de couleur pour les experts
  # Dynamically generate default color values
  #color.values = hue_pal()(length(liste_expert))
  color.values = c("springgreen2","springgreen4","gold2","goldenrod4","hotpink","deeppink","blue","deepskyblue","orange","brown","red","purple")
  color.values=color.values[1:length(liste_expert)]
  names(color.values) = liste_expert
  s.colors=eval(parse(text=paste('c(obs="black",',agreg,'="yellow")',sep=""))) #special colors for the observations and the agregation
  color.values = c(color.values, s.colors)
  width.values=rep(2,length(liste_expert))
  names(width.values)=c(liste_expert)
  s.width=eval(parse(text=paste('c(obs=4,',agreg,'=4)',sep=""))) #special width for the observations and the agregation
  width.values=c(width.values,s.width)
  line.type.values=c("dashed","solid","dashed","solid","dashed","solid","dashed","dashed","dashed","dashed","dashed") # dashed for raw outputs and solid for postprocessed outputs
  line.type.values=line.type.values[1:length(liste_expert)]
  names(line.type.values) = liste_expert
  s.type.values=eval(parse(text=paste('c(obs="solid",',agreg,'="solid")',sep=""))) #solid lines for the observations and the agregation, dashed lines for the others
  line.type.values=c(line.type.values,s.type.values)
  #shape.values=c(1:14)

  plot.df<-ggplot(df, aes(x=date.valid, y=temperature))+
          #theme_gray() +
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) + #to have the lines
          geom_point(aes(color = modele, size = modele ,linetype=modele),show.legend = FALSE) + #to have the points/shapes
          geom_vline(xintercept=as.numeric(datebreaks[1])+12, linetype='dashed', color='grey40',size=2)+#start of the event in chamonix
          geom_vline(xintercept=as.numeric(datebreaks[1])+22, linetype='dashed', color='grey40',size=2)+#end of the event in chamonix
          scale_linetype_manual(values=line.type.values) + #solid or dashed lines, and no legend
          #scale_shape_manual(values = shape.values) + #type of shape for the points
          scale_colour_manual(values=color.values) + #colours
          scale_size_manual(values = width.values) + #sizes
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) + #x axes with datebreaks, %d for days number, %b for month
          labs(x="Month (of the winter 2021-2022)", y= "Temperature (°C)") + #title, x and y axes description
          theme(plot.title = element_text(colour="black", size=25), axis.title.x = element_text(colour="black",size=25),
            axis.title.y = element_text(colour="black",size=25),legend.text=element_text(size=25),
            axis.text.x = element_text(colour="black",size = 25),axis.text.y = element_text(colour="black",size = 25)) +
          theme(legend.title=element_blank())

  ggsave(file=paste(dir_plot,dir_obs_pred,"obs_pred_",agreg,"_window",window,"_",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf",sep=""), width=18, height=7, dpi=600)
  print(paste(dir_plot,dir_obs_pred,"obs_pred_",agreg,"_window",window,"_",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf has been created",sep=""))

  #we derive the error of the modeles
  df[which(df$modele!="obs"),]$temperature=abs(df[which(df$modele!="obs"),]$temperature-df[which(df$modele=="obs"),]$temperature)
  plot.experts<-ggplot(df[which(df$modele!="obs"),], aes(x=date.valid, y=temperature))+
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) + #to have the lines
          geom_point(aes(shape = modele, color = modele, size = modele ,linetype=modele)) + #to have the points/shapes
          scale_linetype_manual(values=line.type.values) + #solid or dashed lines, and no legend
          #scale_shape_manual(values = shape.values) + #type of shape for the points
          scale_colour_manual(values=color.values) + #colours
          scale_size_manual(values = width.values) + #sizes
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) + #x axes with datebreaks
          labs(x="Date", y= "Absolute error (°C)",title="Absolute error of the temperature forecast, Chamonix, lead time 48h") + #title, x and y axes description
          theme(plot.title = element_text(colour="black", size=23,face="bold"), axis.title.x = element_text(colour="black",size=20,face="bold"),
            axis.title.y = element_text(colour="black",size=20,face="bold"),legend.text=element_text(size=20))

  plot.obs<-ggplot(df[which(df$modele=="obs"),], aes(x=date.valid, y=temperature))+
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) +
          geom_point(aes(shape = modele, color = modele, size = modele ,linetype=modele)) +
          scale_linetype_manual(values=line.type.values) +
          #scale_shape_manual(values = shape.values) +
          scale_colour_manual(values=color.values) +
          scale_size_manual(values = width.values) +
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) +
          labs(x="Date", y= "Temperature (°C)") +
          theme(plot.title = element_text(colour="black", size=20,face="bold"), axis.title.x = element_text(colour="black",size=17,face="bold"),
            axis.title.y = element_text(colour="black",size=17,face="bold"),legend.text=element_text(size=17))

  #we put the two plots together
  plot.df <- ggarrange(plot.experts, plot.obs,
                    labels = c("Absolute temperature error of the experts and the agregation  (°C)", "Observations  (°C)"),
                    label.x=0,
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend="right",
                    heights=c(2,1)) #different plot heights

  ggsave(file=paste(dir_plot,dir_obs_pred,"obs_err_pred_",agreg,"_window",window,"_",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf",sep=""), width=15, height=5, dpi=600)
  print(paste(dir_plot,dir_obs_pred,"obs_err_pred_",agreg,"_window",window,"_",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf has been created",sep=""))
}

#plot des obs, des previsions de l'agregation pour plusieurs fenêtre glissantes en fonction du temps/des itérations de l'agregation
#pour une echeance et une station
plot_obs_pred <-function(predictions.algo,Y,date.valid,liste_window,sta,ech,agreg){
  pred.agreg=NULL
  for (iliste in 1:length(predictions.algo)){
    pred.agreg=cbind(pred.agreg,predictions.algo[[iliste]])
  }
  list_agreg<-paste(agreg,liste_window,sep="")
  colnames(pred.agreg)<-list_agreg
  df<-cbind(as.data.frame(date.valid),Y,pred.agreg)
  df<-as.data.frame(df)
  df=df[which(substr(df$date.valid,1,10)==plot.periode[1]):which(substr(df$date.valid,1,10)==plot.periode[2]),] #we only keep the data between plot.periode[1] and plot.periode[2]
  names(df)[names(df) == "Y"] <- "obs"
  df$date.valid=as.POSIXct(df$date.valid)
  df$date.valid=as.Date(df$date.valid) #type date from lubridate
  datebreaks <- seq(df$date.valid[1], df$date.valid[length(df$date.valid)], by="1 month")

  df<-df %>%
      dplyr::select(date.valid, c(all_of(list_agreg),obs)) %>%
      gather(key="modele", value="temperature",-date.valid)

  df <- transform(df, modele=factor(modele, levels=unique(modele))) #without this, ggplot don't plot the lines in the expected order (the order of the dataframe's columns)

  # Dynamically generate default color values
  color.values = hue_pal()(length(list_agreg))
  names(color.values) = list_agreg
  s.colors=eval(parse(text=paste('c(obs="black")',sep=""))) #special colors for the observations
  color.values = c(color.values, s.colors)
  width.values=rep(2,length(list_agreg))
  names(width.values)=c(list_agreg)
  s.width=eval(parse(text=paste('c(obs=2)',sep=""))) #special width for the observations
  width.values=c(width.values,s.width)
  line.type.values=rep("dashed",length(list_agreg))
  names(line.type.values) = list_agreg
  s.type.values=eval(parse(text=paste('c(obs="solid")',sep=""))) #solid lines for the observations
  line.type.values=c(line.type.values,s.type.values)
  #shape.values=c(1:14)

  plot.df<-ggplot(df, aes(x=date.valid, y=temperature))+
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) +
          geom_point(aes(shape = modele, color = modele, size = modele ,linetype=modele)) +
          scale_linetype_manual(values=line.type.values,guide = "none") +
          #scale_shape_manual(values = shape.values) +
          scale_colour_manual(values=color.values) +
          scale_size_manual(values = width.values) +
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) +
          labs(x="Date", y= "Temperature (C°)", title="Temperature forecast of the agregation for different sliding windows and temperature observations") +
          theme(plot.title = element_text(colour="black", size=20,face="bold"), axis.title.x = element_text(colour="black",size=17,face="bold"),
            axis.title.y = element_text(colour="black",size=17,face="bold"),legend.text=element_text(size=17))

  ggsave(file=paste(dir_plot,dir_obs_pred,"obs_pred_",agreg,"_",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf",sep=""), width=15, height=5, dpi=600)
  print(paste(dir_plot,dir_obs_pred,"obs_pred_",agreg,"_",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf has been created",sep=""))
  
  #we derive the error of the modeles
  df[which(df$modele!="obs"),]$temperature=df[which(df$modele!="obs"),]$temperature-df[which(df$modele=="obs"),]$temperature
  plot.agreg<-ggplot(df[which(df$modele!="obs"),], aes(x=date.valid, y=temperature))+
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) +
          geom_point(aes(shape = modele, color = modele, size = modele ,linetype=modele)) +
          scale_linetype_manual(values=line.type.values,guide = "none") +
          #scale_shape_manual(values = shape.values) +
          scale_colour_manual(values=color.values) +
          scale_size_manual(values = width.values) +
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) +
          labs(x="Date", y= "Temperature error (°C)") +
          theme(plot.title = element_text(colour="black", size=20,face="bold"), axis.title.x = element_text(colour="black",size=17,face="bold"),
            axis.title.y = element_text(colour="black",size=17,face="bold"),legend.text=element_text(size=17))

  plot.obs<-ggplot(df[which(df$modele=="obs"),], aes(x=date.valid, y=temperature))+
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) +
          geom_point(aes(shape = modele, color = modele, size = modele ,linetype=modele)) +
          scale_linetype_manual(values=line.type.values,guide = "none") +
          #scale_shape_manual(values = shape.values) +
          scale_colour_manual(values=color.values) +
          scale_size_manual(values = width.values) +
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) +
          labs(x="Date", y= "Obsevations (°C)") +
          theme(plot.title = element_text(colour="black", size=20,face="bold"), axis.title.x = element_text(colour="black",size=17,face="bold"),
            axis.title.y = element_text(colour="black",size=17,face="bold"),legend.text=element_text(size=17))

  #we put the two plots together
  plot.df <- ggarrange(plot.agreg, plot.obs,
                    labels = c("Temperature error of the agregation for different sliding windows (°C)", "Observations (°C)"),
                    label.x=0,
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend="right",
                    heights=c(2,1)) #different plot heights
  
  ggsave(file=paste(dir_plot,dir_obs_pred,"obs_err_pred_",agreg,"_",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf",sep=""), width=15, height=5, dpi=600)
  print(paste(dir_plot,dir_obs_pred,"obs_err_pred_",agreg,"_",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf has been created",sep=""))
}

#plot des obs, des previsions de toutes les agregation pour une fenêtre glissantes en fonction du temps/des itérations de l'agregation
#pour une echeance et une station
plot_obs_agreg <-function(pred.agreg,Y,date.valid,window,sta,ech){
  df<-cbind(as.data.frame(date.valid),Y,pred.agreg)
  df<-as.data.frame(df)
  df=df[which(substr(df$date.valid,1,10)==plot.periode[1]):which(substr(df$date.valid,1,10)==plot.periode[2]),] #we only keep the data between plot.periode[1] and plot.periode[2]
  names(df)[names(df) == "Y"] <- "obs"
  df$date.valid=as.POSIXct(df$date.valid)
  df$date.valid=as.Date(df$date.valid) #type date from lubridate
  datebreaks <- seq(df$date.valid[1], df$date.valid[length(df$date.valid)], by="1 month")

  df<-df %>%
      dplyr::select(date.valid, c(all_of(liste_agreg),obs)) %>%
      gather(key="modele", value="temperature",-date.valid)

  df <- transform(df, modele=factor(modele, levels=unique(modele))) #without this, ggplot don't plot the lines in the expected order (the order of the dataframe's columns)

  # Dynamically generate default color values
  color.values = hue_pal()(length(liste_agreg))
  names(color.values) = liste_agreg
  s.colors=eval(parse(text=paste('c(obs="black")',sep=""))) #special colors for the observations
  color.values = c(color.values, s.colors)
  width.values=rep(2,length(liste_agreg))
  names(width.values)=c(liste_agreg)
  s.width=eval(parse(text=paste('c(obs=2)',sep=""))) #special width for the observations
  width.values=c(width.values,s.width)
  line.type.values=rep("dashed",length(liste_agreg))
  names(line.type.values) = liste_agreg
  s.type.values=eval(parse(text=paste('c(obs="solid")',sep=""))) #solid lines for the observations
  line.type.values=c(line.type.values,s.type.values)
  #shape.values=c(1:14)

  plot.df<-ggplot(df, aes(x=date.valid, y=temperature))+
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) +
          geom_point(aes(shape = modele, color = modele, size = modele ,linetype=modele)) +
          scale_linetype_manual(values=line.type.values,guide = "none") +
          #scale_shape_manual(values = shape.values) +
          scale_colour_manual(values=color.values) +
          scale_size_manual(values = width.values) +
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) +
          labs(x="Date", y= "Temperature (C°)", title="Temperature forecast of the agregations and temperature observations") +
          theme(plot.title = element_text(colour="black", size=20,face="bold"), axis.title.x = element_text(colour="black",size=17,face="bold"),
            axis.title.y = element_text(colour="black",size=17,face="bold"),legend.text=element_text(size=17))

  ggsave(file=paste(dir_plot,dir_obs_pred,"obs_pred_",ech,"_window",window,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf",sep=""), width=15, height=5, dpi=600)
  print(paste(dir_plot,dir_obs_pred,"obs_pred_",ech,"_window",window,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf has been created",sep=""))
  
  #we derive the error of the modeles
  df[which(df$modele!="obs"),]$temperature=df[which(df$modele!="obs"),]$temperature-df[which(df$modele=="obs"),]$temperature
  plot.agreg<-ggplot(df[which(df$modele!="obs"),], aes(x=date.valid, y=temperature))+
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) +
          geom_point(aes(shape = modele, color = modele, size = modele ,linetype=modele)) +
          scale_linetype_manual(values=line.type.values,guide = "none") +
          #scale_shape_manual(values = shape.values) +
          scale_colour_manual(values=color.values) +
          scale_size_manual(values = width.values) +
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) +
          labs(x="Date", y= "Temperature error (°C)") +
          theme(plot.title = element_text(colour="black", size=20,face="bold"), axis.title.x = element_text(colour="black",size=17,face="bold"),
            axis.title.y = element_text(colour="black",size=17,face="bold"),legend.text=element_text(size=17))

  plot.obs<-ggplot(df[which(df$modele=="obs"),], aes(x=date.valid, y=temperature))+
          theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
          geom_line(aes(color = modele, size = modele ,linetype=modele)) +
          geom_point(aes(shape = modele, color = modele, size = modele ,linetype=modele)) +
          scale_linetype_manual(values=line.type.values,guide = "none") +
          #scale_shape_manual(values = shape.values) +
          scale_colour_manual(values=color.values) +
          scale_size_manual(values = width.values) +
          scale_x_date(breaks=datebreaks,labels=date_format("%b")) +
          labs(x="Date", y= "Temperature (°C)") +
          theme(plot.title = element_text(colour="black", size=20,face="bold"), axis.title.x = element_text(colour="black",size=17,face="bold"),
            axis.title.y = element_text(colour="black",size=17,face="bold"),legend.text=element_text(size=17))

  #we put the two plots together
  plot.df <- ggarrange(plot.agreg, plot.obs,
                    labels = c("Temperature error of the agregations (°C)", "Observations (°C)"),
                    label.x=0,
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend="right",
                    heights=c(2,1)) #different plot heights

  ggsave(file=paste(dir_plot,dir_obs_pred,"obs_err_pred_",ech,"_window",window,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf",sep=""), width=15, height=5, dpi=600)
  print(paste(dir_plot,dir_obs_pred,"obs_err_pred_",ech,"_window",window,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf has been created",sep=""))
}





#cf https://github.com/Dralliag/opera/blob/master/R/plot-mixture.R
lance_plot_weights <-function(weights,names.experts,plot_name,date.valid,plot.periode,T,max_experts = 50){
  print(plot_name)
  weights=as.data.frame(weights)
  K <- ncol(weights)
  w.order <- order(colMeans(weights),decreasing = FALSE) # order the experts by increasing mean weight
  col = c("springgreen2","springgreen4","gold2","goldenrod4","hotpink","deeppink","blue","deepskyblue","orange","brown","red","purple")
  col = col[1:length(names.experts)]
  col=col[w.order]
  names(weights) <- names.experts
  names(weights)[names(weights) == "mos.aro" ] <- "ppm.aro"
  names(weights)[names(weights) == "mos.arp" ] <- "ppm.arp"
  names(weights)[names(weights) == "mos.cep" ] <- "ppm.ifs"
  names(weights)[names(weights) == "raw.cep" ] <- "raw.ifs"
  weights <- weights[, w.order]

  if ( length(which(substr(date.valid,1,10)==plot.periode[1])) !=0 & length(which(substr(date.valid,1,10)==plot.periode[2])) != 0 ) {
    id.date=which(substr(date.valid,1,10)==plot.periode[1]):which(substr(date.valid,1,10)==plot.periode[2])
  } else id.date = 1:length(date.valid)
  
  weights=weights[id.date,] #we only keep the data between plot.periode[1] and plot.periode[2]
  date.valid=date.valid[id.date] #we only keep the data between plot.periode[1] and plot.periode[2]
  date.valid=as.POSIXct(date.valid)
  date.valid=as.Date(date.valid) #type date from lubridate
  if (ncol(weights) > max_experts) {
    tmp_weights <- weights[]
    tmp_weights <- cbind(rowSums(tmp_weights[1:(ncol(tmp_weights) - max_experts)]), 
                          tmp_weights[, (ncol(tmp_weights) - max_experts + 1):ncol(tmp_weights)])
    names(tmp_weights)[1] <- "others"
    tmp_K <- min(K, max_experts + 1)
    tmp_cols <- c(rev(col)[1:(tmp_K-1)], "grey")
    
  } else {
    tmp_weights <- weights
    tmp_cols <- rev(col)
    tmp_K <- K
  }
  
  w.summed = apply(tmp_weights,1,sum)
  pdf(paste(dir_plot,dir_weights,plot_name,sep=""))
  #$$plot(w.summed, type = "l", col = col, lwd = 2, axes = F, xlim = c(1, T), 
  plot(w.summed, type = "l", col = col, lwd = 2, xaxt = "n",
        ylim = c(0,max(w.summed)), ylab = "", xlab = "", main = "Weights associated with the experts")
  mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
  abline(h=400, col="black",lwd = 5)
  #$$x.idx <- c(1, 1:T, T:1)
  x.idx <- c(1, 1:length(date.valid),length(date.valid):1) # ok !
  #x.idx <- c(plot.periode[1], as.factor(date.valid),rev(as.factor(date.valid)))
  i.remaining = rep(TRUE, tmp_K)
  for (i in 1:tmp_K) {
    y.idx <- c(0, w.summed, rep(0, T))
    y.idx <- c(0, w.summed, rep(0, length(date.valid)))
    polygon(x = x.idx, y = y.idx, col = tmp_cols[i], border=NA)
    w.summed.old <- w.summed
    w.summed <- w.summed - tmp_weights[, rev(names(tmp_weights))][, i]
    i.remaining[i] <- FALSE
    writeLegend(f = w.summed.old,w.summed,name = rev(colnames(tmp_weights))[i])
  }
  dateBreaks <- seq(date.valid[1], date.valid[length(date.valid)], by="1 month")#$$
  dateBreaks.id=c()
  for (i in 1:length(dateBreaks)){
    dateBreaks.id=c(dateBreaks.id,which(date.valid==dateBreaks[i])) #ajout de l'indice dans date.valid de la valeur dateBreaks[i]
  }
  #dateBreaks.id=coordonnées dans date.valid des dates de dateBreaks#$$
  #axis(side=1, at = c(dateBreaks.id), labels=dateBreaks)#axis(side=1, at = 1:length(date.valid), labels=date.valid)
  axis(side=1, at = dateBreaks.id, labels=dateBreaks,las=1)
  #axis(1)
  #axis(2)
  box()
  names.toWrite <- rev(colnames(tmp_weights))
  names.toWrite[-(1:min(tmp_K,15))] <- ""
  mtext(side = 4, text = names.toWrite,
        #$$at = (1-cumsum(c(tmp_weights[, rev(names(tmp_weights))][T,])))  + 
        #$$  tmp_weights[, rev(names(tmp_weights))][T,]/2, las = 2, col = tmp_cols, cex= 0.5, line = 0.3)
        at = (1-cumsum(c(tmp_weights[, rev(names(tmp_weights))][length(date.valid),])))  + 
          tmp_weights[, rev(names(tmp_weights))][length(date.valid),]/2, las = 2, col = tmp_cols, cex= 0.5, line = 0.3)

  dev.off()
}
