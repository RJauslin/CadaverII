# bioinc_sign <- read.table("bioinc_sign.txt", header=T, sep="\t")



# ---------------------------
# extract each significant OTUs for each group
l <- list()
for(i in 1:nlevels(groups)){
  l[[i]] <- as.matrix(bioinc_sign[which(bioinc_sign[,i] == 1),])
}
inc_values <- do.call(rbind,l) # recreate data.frame
colnames(inc_values)[1:nlevels(groups)] <- levels(groups)



#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

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
# cols <- viridisLite::magma(nlevels(factor(name_otu)))
# cols <- brewer.pal(nlevels(taxo$group), "Paired")
# lev <- levels(factor(name_otu))


colorvertex <- c()
# for(i in 1:nlevels(factor(name_otu))){
for(i in 1:nlevels(factor(taxo$group))){
  colorvertex[which(name_otu == lev[i])] <- cols[i]
}
dat <- data.frame(name_otu = name_otu,rownames(bioinc_sign),col = colorvertex)
rownames(dat) <- as.character(rownames(bioinc_sign))



#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


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


# tikz(file = "network.tex", width =  5.39749, height = 6,standAlone = FALSE)
layout(matrix(c(1,2), nrow=2, byrow=T), c(1,1), c(7,3))
# layout.show(n = 2)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(bi, vertex.label.cex=2, layout=bi_fr, asp=0)
plot(1, type="n", ann=F, axes=F)
legend("center", pch=19, bty="n", horiz=F, x.intersp=0.4,cex = 0.75,
       legend=lev, 
       col=cols,ncol = 2)

