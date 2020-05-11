

# NMDS WITH THE MR6 ( double decostand on the 2 margin )
nmds<-metaMDS(mr6,k=2)

# PREPARE THE DATA.FRAME FOR THE PLOT
data.scores <- as.data.frame(scores(nmds)) # EXCTRACT THE NMDS OUTPUT
data.scores$site <- rownames(data.scores) # RIGHT NAME
var$P <- do.call(rbind,strsplit(var$X1,split = "_"))[,3] ## GET THE PLACET 
data.scores$grp <- as.factor(var$Treatment) # THE TREATMENT GROUP
data.scores$namegrp <- as.factor(var$name_groups) # THE TREATMENT PER TIME

varEnv <- var[,3:11] # ENVIRONEMENTAL VARIABLE
fit_env <-envfit(nmds,varEnv) # ENVFIT


## PREPARE THE DATA TO GGPLOT
envi <- as.data.frame(fit_env$vectors$arrows)
envi$x0 <- rep(0,length(envi$NMDS1))
envi$y0 <- rep(0,length(envi$NMDS1))
envi$NMDS1 <- envi$NMDS1/2
envi$NMDS2 <- envi$NMDS2/2
envi$NMDS1_1 <- envi$NMDS1 + 0.01 
envi$NMDS2_2 <- envi$NMDS2 + 0.01 
envi$name <- rownames(envi)



## EXTTRACT THE CENTROID AND THE CONVEX HULL
df <- data.scores
df$Time <- var$Time
find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls <- ddply(df, "namegrp", find_hull)


# REORDERING THE DATA
data.scores <- data.scores[order(data.scores$namegrp),]

# FUNCTION THAT CALCULATES THE CENTROID
centroid <- function(data,varIndex){
  n <- nlevels(as.factor(data[,varIndex]))
  lvl <- levels(as.factor(data[,varIndex]))
  final <- data.frame(NMDS1 = rep(0,n),
                      NMDS2 = rep(0,n),
                      grp = rep(0,n),
                      fac = rep(0,n))
  for(i in 1:n){
    tmp <-data[which(data[,varIndex] == lvl[i]),]
    final$NMDS1[i] <- sum(tmp$NMDS1)/length(tmp$NMDS1)
    final$NMDS2[i] <- sum(tmp$NMDS2)/length(tmp$NMDS2)
    final$fac[i] <- lvl[i]
    final$grp[i] <- as.character(unique(tmp$grp))
  }
  return(final)
}

# CALCUL THE CENTROID
center <- centroid(data.scores,5)
center$Time <- as.factor(do.call(rbind,strsplit(center$fac,"_"))[,2])

hulls$namegrp <- gsub("_","-",hulls$namegrp)
data.scores$namegrp <- gsub("_","-",data.scores$namegrp)
center$fac <- gsub("_","-",center$fac)
envi$name <- gsub("_","-",envi$name)


p_nmds <- ggplot() +
  geom_polygon(data=hulls,aes(x=NMDS1,y=NMDS2,fill=grp,group=namegrp),alpha=0.30)+
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=grp,shape = grp)) +
  geom_path(data =center,aes(x = NMDS1,y = NMDS2,color = grp))+
  geom_text(data=center,aes(x=NMDS1,y=NMDS2,label = Time,colour = grp,group = as.factor(Time)),size=2.5,vjust=-0.2)+
  geom_segment(data= envi, mapping=aes(x=x0, y=y0, xend=NMDS1, yend=NMDS2), arrow=arrow(length = unit(0.08, "inches")), size=0.4, color="black")+
  geom_text(data=envi,aes(x=NMDS1_1,y=NMDS2_2,label = name),size=2.5,hjust = 0.2)+
  scale_color_manual("Treatments",values = colour_fill)+
  # scale_color_viridis_d("Treatments")+
  scale_shape_manual("Treatments",values = c(15,16,17,18,19))+
  # scale_fill_viridis_d("Treatments")+
  scale_fill_manual("Treatments",values = colour_fill, aesthetics = "fill")+
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="bottom",
        axis.ticks = element_blank(), 
        panel.background = element_blank(),
        plot.background = element_blank())
# print(p_nmds)