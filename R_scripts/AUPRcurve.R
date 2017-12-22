#Plotting the AUPR curve:
#full.dir = "output/yeast/fullrun/affy_TFA_BBSR_1.1_2015-03-16_15-13-01"
#full.dir = "output/yeast/fullrun/affy_TFA_BBSR_1.1_2015-03-16_16-04-18"
#full.dir = "output/yeast/fullrun/affy_TFA_BBSR_1.1_2015-03-16_17-19-08"
full.dir = "output/yeast/fullrun/affy_TFA_BBSR_1.1_2015-03-16_18-45-51"

#load(paste(full.dir, "combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData", sep="/"))
load(paste(full.dir, "combinedconf_frac_tp_50_perm_1--frac_fp_0_perm_1_1.1.RData", sep="/"))
params.and.input <- paste(full.dir, 'params_and_input.RData', sep='/')
load(params.and.input)
#in case of leaveout sets like here, we need to modify the gold standard.
priors.rest = abs(IN$gs.mat) - abs(IN$priors[[1]])
#rows <- apply(abs(abs(IN$gs.mat)), 1, sum) > 0#that's if eval.on.subset=T
#cols <- apply(abs(abs(IN$gs.mat)), 2, sum) > 0#that's if eval.on.subset=T
rows <- apply(abs(priors.rest), 1, sum) > 0
cols <- apply(abs(priors.rest), 2, sum) > 0
#PR = ChristophsPR(order(comb.confs[rows,cols], decreasing=T), abs(IN$gs.mat)[rows,cols])
PR = ChristophsPR(order(comb.confs[rows,cols], decreasing=T), abs(priors.rest)[rows,cols])
#PR.rand = ChristophsPR(sample(order(comb.confs[rows,cols], decreasing=T)), abs(IN$gs.mat)[rows,cols])
PR.rand = ChristophsPR(sample(order(comb.confs[rows,cols], decreasing=T)), abs(priors.rest)[rows,cols])

prec = PR$prec
rec = PR$rec
#AUPR curve for the Leave Out set:
#pdf(file=paste(full.dir, "AUPRcurve.pdf", sep="/"), paper="a4")
#pdf(file=paste(full.dir, "AUPRcurveSGD.pdf", sep="/"), paper="a4")
#pdf(file=paste(full.dir, "AUPRcurveLeaveOut.pdf", sep="/"), paper="a4")
pdf(file=paste(full.dir, "AUPRcurveLeaveOutSGD.pdf", sep="/"), paper="a4")
plot(rec[seq(1, length(rec), length.out=100000)], prec[seq(1, length(rec), length.out=100000)], pch=".", main = c("AUPR=", PR$auc))#and same input with "lines"
lines(rec[seq(1, length(rec), length.out=100000)], prec[seq(1, length(rec), length.out=100000)], pch=".")
dev.off()

#Testing on the entire Platinum Standard:
rows.full <- apply(abs(IN$gs.mat), 1, sum) > 0#that's if eval.on.subset=T
cols.full <- apply(abs(IN$gs.mat), 2, sum) > 0
PR.full = ChristophsPR(order(comb.confs[rows.full,cols.full], decreasing=T), IN$gs.mat[rows.full,cols.full])
prec = PR.full$prec
rec = PR.full$rec
pdf(file=paste(full.dir, "AUPRcurve_fullPS.pdf", sep="/"), paper="a4")
plot(rec[seq(1, length(rec), length.out=100000)], prec[seq(1, length(rec), length.out=100000)], pch=".", main = c("AUPR=", PR.full$auc))#and same input with "lines"
lines(rec[seq(1, length(rec), length.out=100000)], prec[seq(1, length(rec), length.out=100000)], pch=".")
dev.off()

#Testing random AUPR distribution:
rand.auprs = sapply(1:100, function(x) { return(ChristophsPR(sample(order(comb.confs[rows,cols], decreasing=T)), priors.rest[rows,cols])$auc) } )
(PR$auc-mean(rand.auprs))/sd(rand.auprs)# = 36.767 standard deviations away from the mean.

