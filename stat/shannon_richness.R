# Be sure that input.R have been loaded via source()
#
# source("C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/input.R")

var$shannon <- vegan::diversity(mr6, index = "shannon", MARGIN = 1)
var$richness <- rowSums(mr4)

gg_sha_rich <- reshape2::melt(var,id = c("name_groups","Time","Treatment"),measure.vars = c("shannon","richness"))

gg_sha_rich$name_groups <- gsub("_","-",gg_sha_rich$name_groups)

# tikz(file = "shannon.tex", width =  5.39749, height = 6,standAlone = FALSE)
p_shan <- ggplot()+
  geom_point(data = gg_sha_rich,
             aes(x = name_groups,y = value,color = Treatment,shape = Treatment)) +
  facet_wrap(~variable,scale = "free")+
  labs(x = "Treatment per time",
       y = "Observed values",
       # title = "Meuse river",
       fill = "",
       caption = NULL)+
  # scale_color_viridis_d()+
  scale_color_manual("Treatments",values = colour_fill)+
  scale_shape_manual("Treatments",values = c(15,16,17,18,19))+
  theme_biol()
p_shan

