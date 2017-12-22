#Here I will create files that contain predicted taus for all genes for every condition

args <- commandArgs(trailingOnly = TRUE)
folder = args[1]
nclust.conds = as.numeric(args[2])
nclust.genes = as.numeric(args[3])
clust.names.paper = as.logical(args[4])

seed.dirs = paste("prior_seed_", 42:61, sep='')
clusts = paste('clust', 1:nclust.conds, sep='')
taus=c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 200, 250)

library("Matrix")
source("R_scripts/evaluate.R")

gene.clusts = lapply(1:nclust.genes, function(i) {
  t(read.csv(paste(folder, paste("cluster-genes/genes_clust",i,".tsv",sep=""), sep='/'), sep='\t', header=FALSE))
})

expr = read.csv("input/yeast/expression.tsv", sep='\t')

#First we need to assign cluster numbers to the genes that were not in the gold standard (initially, only those genes got assigned clusters).
dist.all = dist(expr) # distances between each pair of genes
dist.mat = as.matrix(dist.all)

gn.clust.full = sapply(rownames(expr), function(gn) {
  which.min(sapply(gene.clusts, function(g.cl) { min(dist.mat[gn,g.cl]) } )) } )

write.table(gn.clust.full, paste(folder,"cluster-genes","gn_clust_full.tsv",sep='/'), sep='\t')

#Now I will need to create files that will contain the predicted taus for every gene depending on its condition and gene: 4 files (one for each condition cluster).

load(paste(folder, "half-life_fitting/prs_by_seed-clust-tau-gclust.RData", sep='/'))

aupr.seed.list <- lapply(pr.seed.list, function(pr.clust.list) {
  aupr.clust.mat.gnlist <- lapply(1:length(gene.clusts), function(k) {
    aupr.clust.mat <- matrix(NA, nrow=length(taus), ncol=nclust.conds, dimnames=list(taus, clusts))
    for(i in 1:length(pr.clust.list)) {
      for(j in 1:length(pr.clust.list[[i]])) {
        aupr.clust.mat[j,i] = pr.clust.list[[i]][[j]][[k]]$auc
      }
    }
    return(aupr.clust.mat)
  } )
  return(aupr.clust.mat.gnlist)
} )

aupr.maxtau.list <- lapply(aupr.seed.list, function(aupr.clust.mat.gnlist) {
  Reduce(cbind, lapply(aupr.clust.mat.gnlist, function(aupr.clust.mat) {
    return(taus[apply(aupr.clust.mat, 2, which.max)])
  } ) )
} )

aupr.maxtau.array3D <- array(unlist(aupr.maxtau.list), c(nclust.conds,nclust.genes,20))

if(clust.names.paper) {
  dimnames(aupr.maxtau.array3D) = list(c("chemostat", "transcr. inhibition", "log-phase growth", "fermentation"), c("nucleic acid metabolism", "translation", "cell wall biogenesis", "no enrichment", "protein catabolism"), c())
}

aupr.median.tab <- apply(aupr.maxtau.array3D, 1:2, median)

clust.taus.list <- lapply(1:nclust.conds, function(i) {
  sapply(gn.clust.full, function(gn.clust) { aupr.median.tab[i,gn.clust] })
})

#write the files that contain optimal taus for each gene for each cluster:
for(i in 1:nclust.conds) {
  write.table(clust.taus.list[[i]], paste(folder,"/cluster-genes/taus-predicted_clust",i,".tsv",sep=''), sep='\t')
}






