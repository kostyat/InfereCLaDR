#April 7, 2016

#The point of this code is to make a plot that consists of N_C by N_G panels (N_C rows - the condition clusters, by N_G columns - the gene clusters), such that each panel contains 20 line plots (from each prior sample) of AUPR on the 50% leave-out set of those genes on that cluster.

#Furthermore, it creates an R data file "prs_by_seed-clust-tau-gclust.RData" that contains all the precision-recall vectors (as well as vectors of true positives (TPs) and false positives (FPs), which are used for ROC-curve calculations) for every bi-cluster

#Then it also make a heatmap of where the peaks occur.

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
  t(read.csv(paste(folder, paste("../cluster-genes/genes_clust",i,".tsv",sep=""), sep='/'), sep='\t', header=FALSE))
})

cat("extracting the optimal tau (RNA half-life * ln(2)) for each bicluster\n")

pr.seed.list = list()
tpfp.seed.list = list()
for(seed.dir in seed.dirs) {
  clust.dirs = list.dirs(paste(folder, seed.dir, sep='/'), recursive=FALSE)
  load(paste(dir(clust.dirs[1], pattern="affy*", full.names=TRUE, recursive=FALSE)[1], "params_and_input.RData", sep='/'))
  cl.gns = lapply(gene.clusts, intersect, rownames(IN$gs.mat))
  pr.clust.list = list()
  tpfp.clust.list = list()
  for(clust.dir in clust.dirs) {
    tau.dirs = dir(clust.dir, pattern="affy*", full.names=TRUE, recursive=FALSE)
    pr.taus.list = list()
    tpfp.taus.list = list()
    for(tau.dir in tau.dirs) {
      load(paste(tau.dir, "combinedconf_frac_tp_50_perm_1--frac_fp_0_perm_1_1.1.RData", sep='/'))
      pr.taus.list[[tau.dir]] = lapply(cl.gns, function(gns) {
        make.ChristophsPR(gs=IN$gs.mat[gns,], priors=IN$priors[[1]][gns,], comb.confs[gns,], TRUE)
      } )
      tpfp.taus.list[[tau.dir]] = lapply(cl.gns, function(gns) {
        make.TP.FP(gs=IN$gs.mat[gns,], priors=IN$priors[[1]][gns,], comb.confs[gns,], TRUE)
      } )
    }  
    pr.clust.list[[clust.dir]] = pr.taus.list
    tpfp.clust.list[[clust.dir]] = tpfp.taus.list
  }
  pr.seed.list[[seed.dir]] = pr.clust.list
  tpfp.seed.list[[seed.dir]] = tpfp.clust.list
}

save(pr.seed.list, file=paste(folder, "prs_by_seed-clust-tau-gclust.RData", sep='/'))
save(tpfp.seed.list, file=paste(folder, "tpfps_by_seed-clust-tau-gclust.RData", sep='/'))

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

#3D cond_clust-by-gene_clust-by-prior_resample array of optimal Taus
aupr.maxtau.array3D <- array(unlist(aupr.maxtau.list), c(nclust.conds,nclust.genes,20))

if(clust.names.paper) {
  dimnames(aupr.maxtau.array3D) = list(c("chemostat", "transcr. inhibition", "log-phase growth", "fermentation"), c("nucleic acid metabolism", "translation", "cell wall biogenesis", "no enrichment", "protein catabolism"), c())
}

aupr.median.tab <- apply(aupr.maxtau.array3D, 1:2, median)
write.table(aupr.median.tab, paste(folder,'median-aupr-table.tsv',sep='/'),sep='\t')

#Making a heatmap
library(gplots)
cols = colorRampPalette(c("white", "red", "dark red"))
png(paste(folder,"aupr-heatmap.png", sep='/'), height=480, width=720)
heatmap.2(log(2)*t(aupr.median.tab), Rowv = rev(as.dendrogram(hclust(dist(t(aupr.median.tab))))), Colv = rev(as.dendrogram(hclust(dist(aupr.median.tab)))), dendrogram="none", col=cols(20), srtCol = -20, offsetRow=0.5, offsetCol=2, adjCol=c(0, 0), cexRow=2, cexCol=2, trace="column", tracecol = "black", vline=NA, key.title = "", key.xlab = "half-lives (minutes)", key.ylab='', keysize=2, margins=c(7,20), key.par=list(mar=c(2,4,0.5,1), cex=1.5))
dev.off()

#Making a big line plot
cols = c('dark blue', 'dark red', 'dark green', 'black', 'orange', 'blue', 'red', 'green', 'purple', 'grey',
'gray45', 'gold', 'goldenrod', 'yellow3', 'yellowgreen', 'saddlebrown', 'thistle3', 'khaki', 'seagreen4', 'mediumpurple')
ymin = min(sapply(aupr.seed.list, function(l) { min(sapply(l, function(l1) {min(l1) })) }))
ymax = max(sapply(aupr.seed.list, function(l) { max(sapply(l, function(l1) {max(l1) })) }))
biclust.num = nclust.conds*nclust.genes
if(clust.names.paper) {
  genes.ord = c(3,4,2,1,5)
  conds.ord = c(3,1,2,4)
} else {
  genes.ord = 1:nclust.genes
  conds.ord = 1:nclust.conds
}

tiff(paste(folder, "auprs-by-cond-and-gene-clust.tif", sep='/'), height=12, width=10, units='in', res=300)
layout(matrix(1:biclust.num, nclust.genes, nclust.conds, byrow = TRUE), widths=rep(1.2,biclust.num), heights=rep(1,biclust.num))
for(gc in genes.ord) {
  for(cn in conds.ord) {
  plot(c(),c(), xlim=c(0,log(2)*max(taus)), ylim=c(ymin,ymax), xlab="RNA half-lifes (minutes)", ylab="AUPR on 50% leave-out", mgp = c(1.6,0.6,0))
  abline(v=seq(0,175,25), h=seq(0.1, 0.4, 0.1), lty=3, col="grey")
  abline(v=log(2)*median(sapply(1:20,function(pr) { taus[which.max(aupr.seed.list[[pr]][[gc]][,cn])] })), lty=2)
    for(pr in 1:20) {
      lines(log(2)*taus, aupr.seed.list[[pr]][[gc]][,cn], col=cols[pr])
      points(log(2)*taus[which.max(aupr.seed.list[[pr]][[gc]][,cn])],max(aupr.seed.list[[pr]][[gc]][,cn]), pch=19, col = cols[pr])
    }
  }
}
dev.off()


