
# source("C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/input.R")
# source("C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/taxo_profile.R")

BiGraph <- function(mr4,var,ass3,char){
    
  # ---- Define big group of control and time 64
  
  groups <- var$name_groups
  groups[which(var$Time == "D0")]  <- "Control"
  groups[which(var$Treatment == "Control")]  <- "Control"
  index <- which(groups == "Control")
  
  # index <- c(index,which(grepl("98",groups)))
  # index <-  c(index,which(grepl("64",groups)))
  index <-  c(index,which(grepl(char,groups)))
  
  groups <- groups[index]
  groups <- factor(groups)
  groups <- factor(groups,levels = c("Control",levels(groups)[-which(levels(groups) == "Control")]))
  var_indic <- var[index,] # extract var matrix
  mr_indic <- mr4[index,] # extract mr_indic matrix
  
  #---- Define design matrix for GLM analyses
  
  design <- model.matrix(~groups)
  dge_indic <- DGEList(counts=t(mr_indic),group = groups)
  
  # More precisely, the filtering keeps OTUs that have count above min.count in n samples,
  # where n is determined by the design matrix.
  keep <- filterByExpr(dge_indic, min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)
  length(which(keep))
  dge_indic <- dge_indic[keep, , keep.lib.sizes=FALSE]
  min(colSums(t(dge_indic$counts)))
  dge_indic$counts
  mr_indic <- t(dge_indic$counts)
  
  #---- edgeR analyses
  
  disp_indic <- estimateGLMRobustDisp(dge_indic,design = design)
  glmfit_indic <- glmFit(disp_indic, design = design)
  lrt_indic <- glmLRT(glmfit_indic,coef = 2:nlevels(groups))
  indc_edge <- topTags(lrt_indic,n = Inf, p.value = 0.05)
  indc_edge <- indc_edge$table
  
  
  
  #---- Indicspecies analyses
  
  bioinc <- multipatt(mr_indic,groups,func = "r.g" ,control = how(nperm=9999))
  bioinc_sign <- bioinc$sign
  bioinc_sign <- bioinc_sign[which(bioinc_sign$p.value < 0.05),]
  
  #---- Intersect the two analyses and redefine the bioinc_sign
  
  index_edge_indval <- intersect(rownames(bioinc_sign),rownames(indc_edge))
  bioinc_sign <- bioinc_sign[index_edge_indval,]
  mr_indic <- mr_indic[,rownames(bioinc_sign)]
  
  #---- save the table
  
 
  
  
  #---- remove index 1 of the table
  bioinc_sign <- bioinc_sign[which(bioinc_sign$index != 1),]
  write.table(bioinc_sign,paste0("C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/results/bioinc_sign",char,".txt"),sep="\t",quote=F)

  # ---------------------------
  # extract each significant OTUs for each group
  l <- list()
  for(i in 1:nlevels(groups)){
    l[[i]] <- as.matrix(bioinc_sign[which(bioinc_sign[,i] == 1),])
  }
  inc_values <- do.call(rbind,l) # recreate data.frame
  colnames(inc_values)[1:nlevels(groups)] <- levels(groups)
  
  
  
  #---- find name for the color
  NAME <- do.call(rbind,strsplit(as.character(ass3$V4),split = "[|]"))[,2]
  name_otu <- c()
  for(i in 1:nrow(bioinc_sign)){
    name_otu <- c(name_otu,NAME[which(rownames(bioinc_sign)[i] == colnames(mr4))])
  }
  
  for(i in 1:length(name_otu)){
    if(!any(name_otu[i] == levels(taxo$group)[2:nlevels(taxo$group)])){
      name_otu[i] <- "Other"    
    }
  }
  
  cols <- pal
  lev <- levels(taxo$group)
  
  colorvertex <- c()
  for(i in 1:nlevels(factor(taxo$group))){
    colorvertex[which(name_otu == lev[i])] <- cols[i]
  }
  dat <- data.frame(name_otu = name_otu,rownames(bioinc_sign),col = colorvertex)
  rownames(dat) <- as.character(rownames(bioinc_sign))
  
  
  
  #---- GRAPH
  
  from <- to <- indval <- c()
  for(i in 1:nlevels(groups)){
    from <- c(from,rep(colnames(inc_values)[i],length(which(inc_values[,colnames(inc_values)[i]] == 1))))
    to <- c(to,rownames(inc_values)[which(inc_values[,colnames(inc_values)[i]]==1)])
    indval <- c(indval,inc_values[which(inc_values[,colnames(inc_values)[i]] == 1),"stat"])
  }
  
  bipart <- data.frame(from = from,to = to, indval = indval)
  bi <- graph.data.frame(bipart,directed=F)
  bi <- simplify(bi, remove.multiple=T, remove.loops=T)
  # V(bi)$color[1:4] <- "grey"
  
  ## Set node labels
  V(bi)$label <- V(bi)$name
  V(bi)$label[grepl("X",V(bi)$label)] = NA
  
  V(bi)$color <- as.character(dat[V(bi)$name, ]$col)
  V(bi)$color[1:nlevels(groups)] <- rep("white",nlevels(groups))
  V(bi)$size <- c(rep(2,nlevels(groups)), rep(5,length(V(bi)$name) - nlevels(groups) ))
  V(bi)$shape <- c(rep("circle",nlevels(groups)),rep("circle",length(V(bi)$name) - nlevels(groups) ))
  V(bi)$frame.color <- V(bi)$color
  bi_fr <- layout_with_fr(bi, niter=9999)
  
  return(list(table = bioinc_sign, bi= bi, bi_fr = bi_fr,lev = lev,cols = cols))
}
  
  # bioinc_08 <- BiGraph(mr4,var,ass3,char = "08")
  # bioinc_15 <- BiGraph(mr4,var,ass3,char = "15")
  # bioinc_21 <- BiGraph(mr4,var,ass3,char = "21")
  # bioinc_64 <- BiGraph(mr4,var,ass3,char = "64")
  # bioinc_98 <- BiGraph(mr4,var,ass3,char = "64")
  # 
  # # tikz(file = "network.tex", width =  5.39749, height = 6,standAlone = FALSE)
  # layout(matrix(c(1,2), nrow=2, byrow=T), c(1,1), c(7,3))
  # # layout.show(n = 2)
  # par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  # plot(test$bi, vertex.label.cex=2, layout=test$bi_fr, asp=0)
  # plot(1, type="n", ann=F, axes=F)
  # legend("center", pch=19, bty="n", horiz=F, x.intersp=0.4,cex = 0.75,
  #        legend= test$lev, 
  #        col=test$cols,ncol = 2)
  # 
