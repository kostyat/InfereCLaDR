args <- commandArgs(trailingOnly = TRUE)
folder = args[1]
nclust.conds = as.numeric(args[2])

cat(paste("Combining predicted interaction confidence scores across the", nclust.conds, "condition clusters\n"))

clusts = dir(folder, patter='affy', full.names=T)
source("R_scripts/evaluate.R")

Calc_Prec_Rec <- function(comb.confs, gs) {
  ord.idx <- order(comb.confs, decreasing=TRUE)

  prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx)))
  rec <- cumsum(gs[ord.idx]) / sum(gs)

  prec <- c(prec[1], prec)
  rec <- c(0, rec)

  auc <- ChristophsAUC(rec, prec)
  return(list(prec=prec, rec=rec, auc=auc))
}


confs.file = 'combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
confs.list <- lapply(clusts, function(clust) {
  load(paste(clust, confs.file, sep='/'))
  return(comb.confs)
} )
load(paste(clusts[1], 'params_and_input.RData', sep='/'))

confs.sum = Reduce('+',confs.list)

write.table(as.matrix(confs.sum), paste(folder,'comb_confs_sum.tsv',sep='/'),sep='\t')

pr = Calc_Prec_Rec(confs.sum, IN$gs.mat)

#Create matrices of 0's and 1's (0=no predicted interaction, 1=predicted interaction) at 0.5, 0.25, and 0.1 precision

# At 50% precision:
max.conf = sort(confs.sum,decreasing=TRUE)[sum(pr$prec>0.5)]
preds = data.frame(which(confs.sum >= max.conf, arr.ind=TRUE))
preds$TF = colnames(confs.sum)[preds[,2]]
preds$gene = rownames(confs.sum)[preds[,1]]
preds$ind =1
pred.mat.pre <- dcast(preds, gene~TF, value.var="ind")
pred.mat = pred.mat.pre[,-1]
rownames(pred.mat) = pred.mat.pre$gene
pred.mat[is.na(pred.mat)] = 0
write.table(pred.mat, paste(folder,"predicted-network-prec_0.5.tsv",sep='/'), sep='\t')

# At 25% precision:
max.conf.25 = sort(confs.sum,decreasing=TRUE)[sum(pr$prec>0.25)]
preds.25 = data.frame(which(confs.sum >= max.conf.25, arr.ind=TRUE))
preds.25$TF = colnames(confs.sum)[preds.25[,2]]
preds.25$gene = rownames(confs.sum)[preds.25[,1]]
preds.25$ind =1
pred.mat.25.pre <- dcast(preds.25, gene~TF, value.var="ind")
pred.mat.25 = pred.mat.25.pre[,-1]
rownames(pred.mat.25) = pred.mat.25.pre$gene
pred.mat.25[is.na(pred.mat.25)] = 0
write.table(pred.mat.25, paste(folder,"predicted-network-prec_0.25.tsv",sep='/'), sep='\t')

# At 10% precision:
max.conf.10 = sort(confs.sum,decreasing=TRUE)[sum(pr$prec>0.10)]
preds.10 = data.frame(which(confs.sum >= max.conf.10, arr.ind=TRUE))
preds.10$TF = colnames(confs.sum)[preds.10[,2]]
preds.10$gene = rownames(confs.sum)[preds.10[,1]]
preds.10$ind =1
pred.mat.10.pre <- dcast(preds.10, gene~TF, value.var="ind")
pred.mat.10 = pred.mat.10.pre[,-1]
rownames(pred.mat.10) = pred.mat.10.pre$gene
pred.mat.10[is.na(pred.mat.10)] = 0
write.table(pred.mat.10, paste(folder,"predicted-network-prec_0.1.tsv",sep='/'), sep='\t')





