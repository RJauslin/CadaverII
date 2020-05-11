
# -----------------------------------------------
rm(list = ls())
# -----------------------------------------------
# load packages
library(gplots)
library(edgeR)
library(indicspecies)
library(igraph)
library(knitr)
library(vegan)
library(ggplot2)
library(readr)
library(knitr)
library(ggrepel)
library(plyr)
library(cluster)
library(viridis)
library(Rtsne)
library(tikzDevice)
library(gridExtra)
library(RColorBrewer)

# -----------------------------------------------
# theme plot of the paper
theme_biol <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family="sans",color = "black",size = 8),
      panel.spacing = unit(2, "lines"),
      # title
      plot.title = element_text(hjust = 0.5,size = 9),
      # axes
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      # legend
      legend.position="bottom",
      legend.title = element_text(size = 9,vjust = +1.0),
      legend.key.size = unit(0.3, "cm"),
      legend.key.width = unit(0.7,"cm") ,
      # background colors
      panel.background=element_blank(),
      panel.border=element_rect(colour = "black",fill = "transparent"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1,size = 5),
      # keep edge black facet_wrap
      # strip.background = element_rect(fill="white"),
      strip.text =element_text(color = "black",size = 8)
    )
}


# workdir <- "C:/Users/jauslinr/switchdrive/LaboBiolSol"
# workdir <- "C:/Users/jauslinr/switchdrive/LaboBiolSol/eukaryote/18S"
workdir <- "C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/data"

# mr <- read.table(file.path(workdir,"RC.mr"),header=T)
# REMOVE SOME SAMPLE
# supp <- c()
# supp <- c(supp,which(rownames(mr) == "E+Faeces+D21+P10+S132+L001"))
# supp <- c(supp,which(rownames(mr) == "E+Feces+deer+dry+MIX+S156+L001"))
# supp <- c(supp,which(rownames(mr) == "E+Feces+deer+wet+MIX+S144+L001"))
# supp <- c(supp,which(rownames(mr) == "E+Feces+goat+2+S168+L001" ))
# supp <- c(supp,which(rownames(mr) == "E+Feces+sheep+MIX+S180+L001"))
# mr <- mr[-supp,]
# rm(supp)

# ass<-read.table(file.path(workdir,"RC.ass"),sep="\t")
# var <-  read_delim("meta3.csv",";", escape_double = FALSE, trim_ws = TRUE)

mr <- read.table(file.path(workdir,"MRfinal.csv"),header=T,sep = ";")
ass<-read.table(file.path(workdir,"ASSfinal.csv"),sep="\t")
var <-  read_delim(file.path(workdir,"meta3.csv"),";", escape_double = FALSE, trim_ws = TRUE)

# Check right number
check <- function(mr,ass){
  val1 <- colSums(mr)
  val2 <- as.numeric(do.call(rbind,strsplit(as.character(ass$V1),split = "_"))[,3])
  return(any(val1 != val2))
}


#---- REMOVE BASCTERIA
ind_MEP <- which(grepl("Metazoa|Embryophyceae|Archaea|Bacteria",ass$V4))
# ind_MEP <- grepl('Metazoa|Embryophyceae|Archaea|Bacteria', ass$V4)
mr2 <- mr[,-ind_MEP]
ass2 <- ass[-ind_MEP,]
rm(ind_MEP)
check(mr2,ass2)

#---- REMOVE FUNGI
ind_FUN <- which(grepl("Opisthokonta|Fungi",ass2$V4))
mr2 <- mr2[,-ind_FUN]
ass2 <- ass2[-ind_FUN,]
rm(ind_FUN)

#---- REMOVE RARE OTU
tmp <- as.numeric(colSums(mr2))
tot <- sum(tmp)
sort(tmp/tot*100,decreasing = T)
index <- which(tmp/tot*100 > 0.01)


mr3 <- mr2[,index]
ass3 <- ass2[index,]
check(mr3,ass3)


#---- REMOVE THE SAMPLE THAT HAVE FAILED (124 occurence of OTUs)
failed <- which(rownames(mr3) == "E+C2+D15+P08")
# failed <- which(grepl("E+C2+D15+P08",as.character(rownames(mr3))))
# mr4 <- mr3
mr4 <- mr3[-failed,]
var <- var[-failed,]
var <- as.data.frame(var)
check(mr4,ass3) ## TRUE, OK because we remove a line of mr3.

#---- RAREFYING BIOLOGICAL COUNT DATA IS STATISTICALLY INADMISSIBLE
#  Despite its current popularity in microbiome analyses rarefying biological count data is statistically inadmissible
# because it requires the omission of available valid data.
#
# [ McMurdie, P.J. & Holmes, S. (2014). Waste not, want not: Why rarefying microbiome data is inadmissible.
# PLoS Comput Biol 10(4): e1003531. doi: 10.1371/journal.pcbi.1003531 ]

mr5 <- decostand(mr4,method = "hellinger",MARGIN = 2)
mr6 <- decostand(mr5,method = "tot",MARGIN = 1)
ass6 <- ass3


colour_fill <- c("Blood" = "red2",
                   "Control" = "grey50",
                   "Faeces"=  "tan4",
                   "Pig" = "hotpink2",
                   "Urine"="gold2")
