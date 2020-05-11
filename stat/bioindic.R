

# ---- Define big group of control and time 64

groups <- var$name_groups
groups[which(var$Time == "D0")]  <- "Control"
groups[which(var$Treatment == "Control")]  <- "Control"
index <- which(groups == "Control")

# index_98 <- c(index,which(grepl("98",groups)))


index <-  c(index,which(grepl("64",groups)))
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


# dge_indic <- calcNormFactors(dge_indic,method = "TMM")
mr_indic <- t(dge_indic$counts)

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


disp_indic <- estimateGLMRobustDisp(dge_indic,design = design)
glmfit_indic <- glmFit(disp_indic, design = design)
lrt_indic <- glmLRT(glmfit_indic,coef = 2:nlevels(groups))
indc_edge <- topTags(lrt_indic,n = Inf, p.value = 0.05)
indc_edge <- indc_edge$table



#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

bioinc <- multipatt(mr_indic,groups,func = "r.g" ,control = how(nperm=9999))
bioinc_sign <- bioinc$sign
bioinc_sign <- bioinc_sign[which(bioinc_sign$p.value < 0.05),]

#---- Intersect the two analyses and redefine the bioinc_sign

index_edge_indval <- intersect(rownames(bioinc_sign),rownames(indc_edge))
bioinc_sign <- bioinc_sign[index_edge_indval,]
mr_indic <- mr_indic[,rownames(bioinc_sign)]

write.table(bioinc_sign,"C:/Users/Raph/switchdrive/LaboBiolSol/Rnw/stat/results/bioinc_sign.txt",sep="\t",quote=F)
