bioinc_sign <- read.table("C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/results/bioinc_sign.txt", header=T, sep="\t")

#---- find name for the color
# NAME <- do.call(rbind,strsplit(as.character(ass3$V4),split = "[|]"))[,2]
index_bioindic <- c()
for(i in 1:nrow(bioinc_sign)){
  index_bioindic <- c(index_bioindic,which(rownames(bioinc_sign)[i] == colnames(mr4)))
}

write.table(ass3[index_bioindic,],"C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/results/bioincTable.csv",sep=";",quote=F)
