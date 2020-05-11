

# CREATION DES LIST
lst<-split(mr6,var$Treatment)
vec<-split(var,var$Treatment)


# PUT CONTROL AT THE END OF EACH TREATMENT OF THE DATA.FRAME
tmp <- vec$Control
tmp2 <- lst$Control
for(i in 1:length(lst)){
  lst[[i]] <- rbind(lst[[i]],tmp2)
  vec[[i]] <- as.data.frame(rbind(vec[[i]],tmp))
}
lst <- lst[-which(names(lst) == "Control")]
vec <- vec[-which(names(vec) == "Control")]


# SCALING ALL TREATMENT (WITH CONTROL FOR EACH TREATMENTS)
for(i in 1:length(vec)){
  vec[[i]][,4:11] <- scale(vec[[i]][,4:11])
}


# FIT WIHT RDA AND ORDISTEP
fit <- list()
for(i in 1:length(vec)){
  fit[[i]] <- ordistep(rda(lst[[i]] ~ Moisture_percent+	TC_percent+	pH+	NH4_ug_g+	NO3_ug_g	+TN_percent+	TOC_percent	+CN,vec[[i]][,4:11]))
  names(fit)[i] <- names(vec)[i]
}


# BIG LIST WITH ALL THE INFORMATIONS NEEDED FOR THE GGPLOT
summ <- list()
for(i in 1:length(vec)){
  summ[[i]] <- list()
  summ[[i]]$summary <- summary(fit[[i]])
  
  
  # extract site of the summary
  summ[[i]]$site <- cbind(summ[[i]]$summary$sites[,1:2],
                          Treatment2 = rep(names(vec)[i],length(summ[[i]]$summary$sites[,1])),
                          Treatment = vec[[i]]$Treatment,
                          Time = vec[[i]]$Time,
                          Namegroup = vec[[i]]$name_groups
  )
  colnames(summ[[i]]$site) <- c("RDA1","RDA2","Treatment2","Treatment","Time","name_group")
  
  # extract species of the summary
  summ[[i]]$species <- cbind(summ[[i]]$summary$species[,1:2],
                             Treatment = rep(names(vec)[i],length(summ[[i]]$summary$species[,1])))
  colnames(summ[[i]]$species) <- c("RDA1","RDA2","Treatment")
  
  # extract biplot of the summary
  summ[[i]]$biplot <-  cbind(summ[[i]]$summary$biplot[,1:2],
                             Treatment = rep(names(vec)[i],length(summ[[i]]$summary$biplot[,1])))
  colnames(summ[[i]]$biplot) <- c("RDA1","RDA2","Treatment2")
  
  # extract pourc of the summary
  summ[[i]]$pourc <- data.frame(matrix(rep(0,3),ncol = 3,nrow = 1))
  colnames(summ[[i]]$pourc) <- c("RDA1","RDA2","Treatment2")
  summ[[i]]$pourc$RDA1 <- round(summ[[i]]$summary$cont$importance[2,1]*100,3)
  summ[[i]]$pourc$RDA2 <- round(summ[[i]]$summary$cont$importance[2,2]*100,3)
  summ[[i]]$pourc$Treatment2 <- names(vec)[i]
  
  
  names(summ)[i] <- names(vec)[i]
}



# rbind THE DATA.FRAME AND ADD SOME STUFF TO THE PLOT
datSite <- datSpecies <- datBiplot <- datPourc <-  data.frame()
for(i in 1:length(summ)){
  datSite <- rbind(datSite,summ[[i]]$site)
  datSpecies <- rbind(datSpecies,summ[[i]]$species)
  datBiplot <- rbind(datBiplot,summ[[i]]$biplot)
  datPourc <- rbind(datPourc,summ[[i]]$pourc)
}

datBiplot$RDA1 <- as.numeric(as.character(datBiplot$RDA1))/4
datBiplot$RDA2 <- as.numeric(as.character(datBiplot$RDA2))/4
datBiplot$x0 <- rep(0,length(datBiplot[,1]))
datBiplot$y0 <-  rep(0,length(datBiplot[,1]))
datBiplot$name <- rownames(datBiplot)



datSite$RDA1 <- as.numeric(as.character(datSite$RDA1))
datSite$RDA2 <- as.numeric(as.character(datSite$RDA2))

datPourc$RDA1 <- paste("F1 : ", datPourc$RDA1,"",sep = "")
datPourc$RDA2 <- paste("F2 : ", datPourc$RDA2,"",sep = "")
datPourc$x <- rep(-0.2,4)
datPourc$y <- rep(0.3,4)




## EXTTRACT THE CENTROID AND THE CONVEX HULL
# the names must be NULL for the ddply function ...
names(datSite$Treatment2) <- NULL
names(datSite$Treatment) <- NULL
names(datSite$Time) <- NULL
names(datSite$name_group) <- NULL
find_hull <- function(df) df[chull(df$RDA1, df$RDA2), ]
hulls <- ddply(datSite, c("name_group"), find_hull)

# FUNCTION THAT CALCULATES THE CENTROID
centroid <- function(data,varIndex){
  n <- nlevels(as.factor(data[,varIndex]))
  lvl <- levels(as.factor(data[,varIndex]))
  final <- data.frame(RDA1 = rep(0,n),
                      RDA2 = rep(0,n),
                      Treatment = rep(0,n),
                      name_group = rep(0,n))
  for(i in 1:n){
    tmp <-data[which(data[,varIndex] == lvl[i]),]
    final$RDA1[i] <- sum(tmp$RDA1)/length(tmp$RDA1)
    final$RDA2[i] <- sum(tmp$RDA2)/length(tmp$RDA2)
    final$name_group[i] <- lvl[i]
    final$Treatment[i] <- as.character(unique(tmp$Treatment))
  }
  return(final)
}

# CALCUL THE CENTROID
center <- centroid(hulls,6)
center <- center[-which(center$Treatment == "Control"),] # not desired for the plot

colnames(center) <- c("RDA1","RDA2","Treatment2","name_group")
center$Time <- as.numeric(do.call(rbind,strsplit(center$name_group,split = "_D"))[,2]) # add time to CENTER


# colour_fill <- c("Blood" = "lightskyblue4",
#                  "Control" = "coral4",
#                  "Faeces"=  "royalblue4",
#                  "Pig" = "ivory4",
#                  "Urine"="olivedrab4")


datSite$name_group <- gsub("_","-",datSite$name_group)
center$name_group <- gsub("_","-",center$name_group)
datBiplot$name <- gsub("_","-",datBiplot$name)

# tikz(file = "rda.tex", width =  5.39749, height = 6,standAlone = FALSE)
p_rda <- ggplot() +
  geom_polygon(data = datSite,aes(x = RDA1,y = RDA2,fill = Treatment,group = name_group),alpha = 0.3)+
  geom_point(data=datSite,aes(x=RDA1,y=RDA2,color=Treatment,shape = Treatment))+
  geom_path(data =center,aes(x = RDA1,y = RDA2,color = Treatment2))+
  geom_text(data=center,aes(x=RDA1,y=RDA2,label = Time,colour = Treatment2,group = as.factor(Time)),size=3.0,vjust=-0.2)+
  geom_segment(data= datBiplot,aes(x=x0, y=y0, xend=RDA1, yend=RDA2,color = Treatment2), arrow=arrow(length = unit(0.1, "inches")),colour = "black", size=0.3)+
  geom_text(data=datBiplot,aes(x=RDA1+0.01,y=RDA2+0.01,label = name),size=2.5,hjust = 0.2)+
  geom_text(data=datPourc,aes(x=x,y=y,label = RDA1),size=2.5,hjust = 0.2)+
  geom_text(data=datPourc,aes(x=x,y=y-0.04,label = RDA2),size=2.5,hjust = 0.2)+
  facet_wrap(~Treatment2)+
  scale_color_manual("Treatments",values = colour_fill)+
  scale_shape_manual("Treatments",values = c(15,16,17,18,19))+
  scale_fill_manual("Treatments",values = colour_fill, aesthetics = "fill")+
  xlab("RDA1") + ylab("RDA2")+
  geom_vline(xintercept = 0,linetype="dashed",size = 0.1)+
  geom_hline(yintercept = 0,linetype="dashed",size = 0.1)+
  theme_biol() +
  # ggtitle("18S") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="bottom",
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
# print(p_rda)