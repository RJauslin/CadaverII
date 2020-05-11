
# Be sure that input.R have been loaded via source()
#
# source("C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/input.R")
 

mr_pr_ab <- ifelse(mr4==0,0,1) # CHANGE THE COMM MATRIX INTO OCCURENCE MATRIX  1 - 0
resp <-rowSums(mr_pr_ab)
names(resp) <- var$X1
var$resp <- resp


var$Time2 <- as.numeric(do.call(rbind,strsplit(var$Time,"D"))[,2])


# tikz(file = "loess.tex", width =  5.39749, height = 6,standAlone = FALSE)
p_loess <-ggplot(data = var,aes(x = Time2,y = resp,color = Treatment,shape = Treatment)) +
  geom_point()+
  geom_smooth(aes(group = Treatment,color = Treatment),se = F,alpha = 0.2)+
  # geom_smooth(aes(group = Treatment),method = lm,se = F)+
  scale_color_manual("Treatment",values = colour_fill)+
  # ggtitle("OTU's number evolution over time")+
  xlab("Time (days)")+ ylab("Number of OTUs")+
  scale_shape_manual("Treatment",values = c(15,16,17,18,19))+
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
