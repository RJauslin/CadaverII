
var$place <-  paste(var$name_groups,do.call(rbind,strsplit(as.character(var$X1),split = "[_]"))[,3],sep = "_")

f <- do.call(rbind,strsplit(as.character(ass3$V4),split = "[|]"))[,2]
f <- gsub("_","-",f)

mrtest <- mr6

taxo <- c()
for(i in 1:nrow(mrtest)){
  tmp <- aggregate(as.data.frame(as.numeric(mrtest[i,])),by = list(f),FUN = "sum")
  tmp[,2] <- tmp[,2]/sum(tmp[,2])*100
  if(any(tmp[,2] < 1)){
    index <- which(tmp[,2] <1)
  }
  tmp[index,1] <- "Other"
  tmp[,3] <- rep(var$place[i],nrow(tmp))
  tmp <- aggregate(tmp[,2],by = list(factor(tmp[,1]),tmp[,3]),FUN = "sum")
  colnames(tmp) <- c("group","treat_time","relab")
  taxo <- rbind(taxo,tmp)
}


taxo$group <- factor(taxo$group,levels = c("Other",levels(taxo$group)[-which(levels(taxo$group) == "Other")]))
pal <- brewer.pal(nlevels(taxo$group), "Paired")
pal <- c("grey",pal)
taxo$treat_time <- gsub("_","-",taxo$treat_time)

# tikz(file = "taxoprofile.tex", width =  5.39749, height = 6,standAlone = FALSE)
p_taxo <- ggplot()+
  geom_bar(data = taxo,aes(x = treat_time,y = relab,fill = group),color = "black",stat="identity",size = 0.01)+
  # scale_fill_brewer(palette = "Set2")+
  scale_fill_manual(values= pal)+
  labs(x = "Treatment per time",
       y = "Relative abundance",
       # title = "Meuse river",
       fill = "",
       caption = NULL)+
  theme_biol()+
  theme(
    panel.border=element_blank()
  )
# print(p_taxo)
