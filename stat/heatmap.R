
## CALCULATE RELATIVE ABUNDANCE OF 
mr_edge_indval <- mr_indic
# mr_edge_indval <- mr_edge_indval[index,]
mr_edge_indval <- decostand(mr_edge_indval,method = "hellinger",MARGIN = 2)
mr_edge_indval <- decostand(mr_edge_indval,method = "tot",MARGIN = 1)

mr_edge_indval <- aggregate(mr_edge_indval,by = list(groups),FUN = "sum")
rownames(mr_edge_indval) <- mr_edge_indval[,1]
mr_edge_indval <- mr_edge_indval[,-1]
mr_edge_indval <- t(mr_edge_indval/rowSums(mr_edge_indval) * 100)


phylum <- mr_edge_indval
phylum <- aggregate(phylum,by = list(dat$name_otu),FUN = "sum")
rownames(phylum) <- phylum[,1]
phylum <- phylum[,-1] 


# colSums(mr_edge_indval)
# mr_edge_indval <- t(mr_edge_indval[index,])
# rownames(mr_edge_indval) <- groups[index]

rownames(mr_edge_indval)   <- paste(rownames(mr_edge_indval),dat$name_otu)
colors <- colorRampPalette(c("grey90", "darkblue", "darkred"))

# mr_edge_indval<-  t(mr_edge_indval)
heatmap.2(as.matrix(log2(mr_edge_indval+1)),
          col=colors,
          margins=c(5,12),
          trace="none",
          dendrogram = "none",
          labRow = NA,
          cexCol=1,
          srtCol=45,
          key.xlab=expression(paste("Relative OTU abundance"," (log"[2]," CPM)")))


heatmap.2(as.matrix(log2(phylum+1)), col=colors, margins=c(5,12),scale="none", trace="none"
          ,density.info="none",key.title=NA,
          srtCol=45,
          cexCol = 1,
          labRow=rownames(phylum),
          labCol = colnames(phylum),
          dendrogram="none",
          key.xlab=expression(paste("Relative OTU abundance"," (log"[2]," CPM)")))

