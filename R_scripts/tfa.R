
# fix prior known interactions so the pseudoinverse does not fail
fix.pki <- function(pki) {
  diag(pki)[apply(pki != 0, 2, sum) == 0] <- 1
  # check for identical columns in pki
  m <- ncol(pki)
  fix.me <- T
  while (fix.me) {
    fix.me <- F
    for (i in which(duplicated(pki, MARGIN=2))) {
      for (j in 1:(i-1)) {
        if(all(pki[, i] == pki[, j]) == TRUE) {
          pki[i, i] <- 1
          pki[j, j] <- 1
          cat('Identical columns in rf (', i, j, ') fixed via self-interactions.\n')
          fix.me <- T
        }
      }
    }
  }
  return(pki)
}


# if noself is TRUE, all self-regulatory interactions are removed
# however, if dup.self is TRUE, self interactions for TFs that other TFs with the 
# exact same set of interactions in the prior are kept
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
  tfs <- intersect(tfs, rownames(prior))
  if (noself) {
    diag(prior[tfs, tfs]) <- 0
  }
  
  require('corpcor')
  activities <- pseudoinverse(prior) %*% exp.mat.halftau
  dimnames(activities) <- list(colnames(prior), colnames(exp.mat.halftau))
  
  has.no.act <- names(which(apply(prior != 0, 2, sum) == 0))
  has.no.act <- intersect(has.no.act, rownames(exp.mat))
  
  activities[has.no.act, ] <- exp.mat[has.no.act, ]
  
  return(activities)
}


tfa.bs <- function(prior, res.mat, des.mat) {
  K <- 50
  act.bs <- array(NA, dim=c(ncol(prior), ncol(res.mat), K))
  cond <- ncol(res.mat)
  for (k in 1:K) {
    selected <- sample(cond, replace=T)
    activities <- tfa(prior, res.mat[, selected], des.mat[, selected])
    act.bs[,selected,k] <- activities
  }
  ret <- apply(act.bs, 1:2, median, na.rm=TRUE)
  dimnames(ret) <- dimnames(activities)
  return(ret)
}

tfa.nn <- function(prior, res.mat, des.mat, cores) {
  require('nnls')
  diag(prior) <- 0
  #cor.mat <- cor(t(res.mat), t(res.mat[colnames(prior), ]))
  #prior <- prior * sign(cor.mat)
  #prior.fixed <- fix.pki(prior)
  
  #dn <- dimnames(prior)
  #prior <- refit.betas.mc(des.mat[colnames(prior), ], res.mat, prior)[, -1]
  #dimnames(prior) <- dn
  prior.fixed <- fix.pki(prior)
  
  #abs.col.avg <- apply(abs(prior.fixed), 2, sum) / apply(prior.fixed != 0, 2, sum)
  #abs.col.avg <- apply(prior.fixed, 2, sum) / apply(prior.fixed != 0, 2, sum)
  #prior.fixed <- t(t(prior.fixed) / abs.col.avg)
  
  des.mat <- matrix(0, ncol(prior), ncol(res.mat), dimnames=list(colnames(prior), colnames(res.mat)))
  for (i in 1:ncol(res.mat)) {
    des.mat[, i] <- nnls(prior.fixed, res.mat[, i])$x
  }
  return(des.mat)
}

collapse.activities <- function(activities, prior) {
  tfs <- colnames(prior)
  prior.noself <- prior
  diag(prior.noself[tfs, tfs]) <- 0
  noself.targets <- apply(prior.noself != 0, 2, sum)
  
  dup <- which(duplicated(t(prior)) & apply(prior != 0, 2, sum) > 0)
  while (length(dup) > 0) {
    dup.cols <- which(apply(prior, 2, function(x) all(x==prior[, dup[1]])))
    exemplar <- dup.cols[order(noself.targets[dup.cols], decreasing=TRUE)[1]]
    print('collapse activities')
    print(dup)
    print(dup.cols)
    print(exemplar)
    for (dup.col in dup.cols) {
      #activities[dup.cols, ] <- apply(activities[dup.cols, ], 2, mean)
      activities[dup.col, ] <- activities[exemplar, ]
    }
    dup <- setdiff(dup, dup.cols)
  }
  return(activities)
}

# given a design and response matrix, refit the betas
# useful if we don't know whether the old betas came from scaled design and 
# response matrices
refit.betas.mc <- function(X, Y, betas.old, cores) {
  X <- rbind(1, X)
  beta.rows <- mclapply(1:nrow(Y), refit.one, X, Y, betas.old, mc.cores=cores)
  beta <- matrix(unlist(beta.rows), nrow(Y), byrow=TRUE)
  return(beta=beta)
}

refit.one <- function(i, X, Y, betas.old) {
  K <- nrow(X)
  beta <- rep(0, K)
  selected <- c(TRUE, betas.old[i, ] != 0)
  x <- t(matrix(X[selected, ], sum(selected)))
  coefs <- as.numeric(solve(crossprod(x), crossprod(x, Y[i, ])))
  beta[selected] <- coefs
  return(beta)
}

