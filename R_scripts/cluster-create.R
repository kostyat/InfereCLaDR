pc = 16 #the number of principal components dimensions to use for dimensionality reduction
args <- commandArgs(trailingOnly = TRUE)
nclust = as.numeric(args[1]) #number of expression data condition clusters
nclust.genes = as.numeric(args[2]) #number of expression data gene clusters
out.dir = args[3]
library(parallel)
mc.cores = detectCores()/2

expr = read.csv('input/yeast/expression.tsv',sep='\t')

### GENE CLUSTERING ###
cat(paste("Creating the", nclust.genes, "gene clusters...\n"))

#First, normalize (or is it already normalized)
expr.scaled = t(scale(t(expr)))

#Second, reduce the number of genes down to only the genes in the GS (for more even clustering, in order to get better RNA half-life predictions for every cluster):
gs = read.csv("input/yeast/gold-standard.tsv",sep='\t')

expr.gs.scaled = expr.scaled[intersect(rownames(expr),rownames(gs)),]

d.genes <- dist(as.matrix(expr.gs.scaled))
hc.genes=hclust(d.genes)
#This creates a table of mappings from genes to cluster numbers:
ct.genes = cutree(hc.genes, k=nclust.genes)

#Write lists of genes for each gene cluster into a file:
gene.dir = paste("output", out.dir, "cluster-genes", sep='/')
if(!dir.exists(gene.dir)) { dir.create(gene.dir, recursive=TRUE) }
cat(paste("Writing the gene clusters to", gene.dir, "\n"))
for(i in 1:nclust.genes) {
  filename = paste(gene.dir, "/genes_clust", i, ".tsv", sep="")
  write.table(names(which(ct.genes==i)), filename, sep='\n', row.names=FALSE, col.names=FALSE)
}


### CONDITION CLUSTERING ###
cat(paste("Creating the", nclust, "condition clusters...\n"))

#Remove outliers before condition clustering (this is important to get more even clusters):
d.expr = dist(t(expr))
fit = hclust(d.expr)
#if there is one cluster that really dominates the expression dataset (i.e. if more than 95% of the samples are in one cluster), we keep only that cluster
dominates=TRUE
k=1
cutoff=0.95
while(dominates) {
  groups = cutree(fit, k)
  group.proportions = table(groups)/sum(table(groups))
  if(sum(group.proportions>cutoff)==1)  {
    k=k+1
  } else {
    dominates=FALSE
  }
}
groups = cutree(fit, k-1) #k-1=5 was used in the InfereCLaDR paper.
group.proportions = table(groups)/sum(table(groups))
expr.no.col.outliers = expr[,which(groups == which(group.proportions>cutoff))]
  
### Cluster conditions

expr.scaled = t(scale(t(expr.no.col.outliers)))
#Reduce dimensionality via PCA (Principal Components Analysis):
pc.nocol = prcomp(t(expr.scaled))
comp <- data.frame(pc.nocol$x[,1:pc])

set.seed(6) #always use the same random seed to keep the order of clusters the same. This is the random seed that corresponds to the results in the InfereCLaDR paper (Tchourine et al 2018).
k <- kmeans(comp, nclust, nstart=25, iter.max=1000)
#Create nclust expression matrices with the expression levels of all genes for each condition cluster:
expr.list <- mclapply(1:nclust, function(i) { return(expr.no.col.outliers[,which(k$cluster == i)]) }, mc.cores = mc.cores)

#Create the corresponding meta data files:
source("R_scripts/design_and_response.R")
meta=read.csv('input/yeast/meta_data.tsv',sep='\t')
meta.list <- mclapply(1:nclust, function(i) { return(meta.reduce(meta, colnames(expr.list[[i]]))) }, mc.cores = mc.cores)

write.expr.meta.files <- function(expr.list, meta.list, folder.name) {
  dir.create(file.path(folder.name), showWarnings = FALSE)
  for(i in 1:length(expr.list)) {
    expr.file.name = paste(folder.name, paste("expr_clust", i, ".tsv", sep=""), sep = "/")
    meta.file.name = paste(folder.name, paste("meta_clust", i, ".tsv", sep=""), sep = "/")
    write.table(expr.list[[i]], file = expr.file.name, sep="\t")
    write.table(meta.list[[i]], file = meta.file.name, sep="\t")
  }
}
#Write the expression and meta data files for each cluster:
cond.dir = paste("output", out.dir, "clusters-conditions", sep='/')
cat(paste("Writing the condition clusters to", cond.dir, "\n"))
write.expr.meta.files(expr.list, meta.list, cond.dir)



