require('Matrix')
require('ggplot2')
require('multicore')
require('corpcor')

# if noself is TRUE, all self-regulatory interactions are removed
# however, if dup.self is TRUE, self interactions for TFs that have other TFs 
# with the exact same set of interactions in the prior are kept
# the motivation for dup.self=TRUE is that TFs with identical prior should have 
# identical TFA
tfa <- function(prior, exp.mat, exp.mat.halftau, noself=TRUE, dup.self=TRUE) {
  duplicates <- c()
  if (dup.self) {
    has.p <- apply(prior != 0, 2, sum) > 0
    duplicates <- duplicated(t(prior[, has.p])) |
                  duplicated(t(prior[, has.p]), fromLast = TRUE)
    duplicates <- colnames(prior)[duplicates]
  }
  tfs <- setdiff(colnames(prior), duplicates)
  if (noself) {
    diag(prior[tfs, tfs]) <- 0
  }
  
  activities <- pseudoinverse(prior) %*% exp.mat.halftau
  dimnames(activities) <- list(colnames(prior), colnames(exp.mat.halftau))
  
  has.no.act <- names(which(apply(prior != 0, 2, sum) == 0))
  
  activities[has.no.act, ] <- exp.mat[has.no.act, ]
  
  return(activities)
}

# sub-sample entries from the priors matrix; either a perc percent, or with replacement
subsample <- function(priors, perc=NULL) {
  if (is.null(perc)) {
    prior.order <- sample(which(priors != 0), replace=TRUE)
    prior.order <- prior.order[!duplicated(prior.order)]
    n.priors <- length(prior.order)
  } else {
    prior.order <- sample(which(priors != 0))
    n.priors <- floor(sum(priors != 0) * perc / 100)
  }
  p.mat <- priors * 0
  p.mat[prior.order[1:n.priors]] <- priors[prior.order[1:n.priors]]
  return(p.mat)
}


tfa.stability <- function(design.mat, response.mat, prior.mat, N=32, seed=42, cores=1) {
  set.seed(seed)
  
  # to speed things up later on, we only care about genes and TFs that are part
  # of the prior network
  p.genes <- intersect(rownames(prior.mat)[apply(prior.mat!=0, 1, sum) > 0], rownames(design.mat))
  p.tfs <- intersect(colnames(prior.mat)[apply(prior.mat!=0, 2, sum) > 0], rownames(design.mat))
  goi <- union(p.genes, p.tfs)
  
  # distance between TFA profiles is correlation based
  # if no signed priors are used, use absolute correlation
  my.sum <- mean
  if (sum(prior.mat < 0) == 0) {
    cat('looks like this is an unsigned prior - will use absolute correlation\n')
    my.sum <- function(x) mean(abs(x))
  }
  
  cat('sub-sample prior matrix: ')
  prior.list <- list()
  for (i in 1:N) {
    cat(i, '')
    prior.list[[i]] <- subsample(prior.mat[goi, p.tfs])
  }
  
  cat('\ncompute TFA for all sub-sampled prior matrices\n')
  tfa.ss <- mclapply(prior.list, function(x) tfa(x, design.mat[goi, ], response.mat[goi, ]), mc.cores=cores)
  
  # get distance between the different tfa matrices
  all.cors <- c()
  cat('get distance between tfa matrices\n')
  combs <- subset(expand.grid(1:N, 1:N), Var2 > Var1)
  cor.vecs <- mclapply(1:nrow(combs), function(i) diag(cor(t(tfa.ss[[combs[i, 1]]]), t(tfa.ss[[combs[i, 2]]]))), mc.cores=cores)
  for (i in 1:nrow(combs)) {
    all.cors <- cbind(all.cors, cor.vecs[[i]])
  }
  rownames(all.cors) <- rownames(tfa.ss[[1]])
  
  # summarize the results
  avg.cor <- apply(all.cors, 1, my.sum)
  res <- data.frame(rownames(tfa.ss[[1]]), apply(prior.mat[, rownames(tfa.ss[[1]])]!=0, 2, sum), avg.cor)
  colnames(res) <- c('ID', 'targets.in.prior', 'avg.cor')
  
  return(res)
}

plot.stabilities <- function(df, fname, is.sig=NA) {
  df$targets.in.prior <- as.factor(df$targets.in.prior)
  df$higher.than.random <- is.sig
  g <- ggplot(df, aes(x=targets.in.prior, y=avg.cor)) + 
    #geom_boxplot(color='grey50') + 
    geom_point(size=3, aes(fill=higher.than.random), position = position_jitter(width = .25), shape=21, color='black') + 
    xlab('Number of targets in prior') + 
    ylab('Stability') + 
    theme(axis.title.y = element_text(angle=0)) +
    theme(axis.text.x = element_text(angle=270, hjust=0, size=rel(0.7)))
  #dev.new()
  pdf(file=fname, width=14, height=7, useDingbats=FALSE)
  print(g)
  dev.off()
}


# set how many cores to use; don't set this too high, if you don't have enough memory
#CORES <- 12 
CORES <- 1

# here I am loading design matrix, response matrix and prior matrix from a previous 
# Inferelator run
#load('output/yeast/GScrossCompare1/affy_TFA_BBSR_1.1_2015-03-22_13-12-02/params_and_input.RData')
#load('output/yeast/affy_TFA_BBSR_1.1_2015-03-26_14-00-44/params_and_input.RData')
load('output/yeast/Platinum/affy_TFA_BBSR_1.1_2015-04-02_16-02-51/params_and_input.RData')
design.matrix <- IN$final_design_matrix
response.matrix <- IN$final_response_matrix_halftau
prior.mat <- IN$priors.mat

# get the TFA stabilities
stability.result <- tfa.stability(design.matrix, 
                                  response.matrix, 
                                  prior.mat,
                                  cores=CORES)
#plot.stabilities(stability.result, 'tmp_tfa_stability_plot.pdf')






# it is a good idea to compare the results to those of randomized versions
# of the prior matrix
N <- 12  # number of random prior networks
r.avg.cors <- c()
for (i in 1:N) {
  # randomize prior net
  r.prior.mat <- prior.mat[sample(nrow(prior.mat)), ]
  rownames(r.prior.mat) <- rownames(prior.mat)
  
  cat('run tfa.stability for random prior', i, '\n')
  r.stab.res <- tfa.stability(design.matrix, 
                              response.matrix, 
                              r.prior.mat,
                              cores=CORES)
                              
  r.avg.cors <- cbind(r.avg.cors, r.stab.res[stability.result$ID, 'avg.cor'])
}
r.larger <- sapply(1:nrow(r.avg.cors), function(j) sum(r.avg.cors[j, ] > stability.result$avg.cor[j]))
# let's call any tfa stability significant, if no random prior matrix resulted in a higher stability
plot.stabilities(stability.result, 'tmp_tfa_stability_plot_significance_platinum.pdf', is.sig=r.larger==0)



