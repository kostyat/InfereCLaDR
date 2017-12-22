source("R_scripts/evaluate.R")

plot.auprs <- function(in.dir, out.dir, length.out = 100000) {
  params.and.input <- paste(in.dir, 'params_and_input.RData', sep='/')
  if (!file.exists(params.and.input)) {
    cat('No params_and_input.RData - skipping', full.dir, '\n')
    return()
  }
  load(params.and.input)
  
  files <- list.files(in.dir, "combinedconf_.+\\.RData$")
  print(file.path(in.dir, files))
  
  gs <- IN$gs.mat
  file.name = paste(paste(file.path(out.dir),gsub("/", "-",in.dir), sep='/'),".pdf",sep="")
  pdf(file.name, height=15, width=15, paper="a4")
  for (res.file in files) {
    load(file.path(in.dir, res.file))
    if(PARS$deg.rates) {
      comb.confs = comb.confs[,-1]
    }
    rows <- rep(TRUE, nrow(gs))
    cols <- rep(TRUE, ncol(gs))
    if(PARS$eval.on.subset) {
      rows <- apply(gs, 1, sum) > 0
      cols <- apply(gs, 2, sum) > 0
    }
    PR = ChristophsPR(order(comb.confs[rows,cols], decreasing=T), gs[rows,cols])
    prec = PR$prec
    rec = PR$rec
    plot(rec[seq(1, length(rec), length.out=length.out)], prec[seq(1, length(rec), length.out=length.out)], pch=".", main = c("AUPR=", PR$auc))#and same input with "lines"
    lines(rec[seq(1, length(rec), length.out=length.out)], prec[seq(1, length(rec), length.out=length.out)], pch=".")
    if(PARS$perc.tp < 100) {
      gs.li = abs(IN$priors[[which(res.file %in% files)]])
      gs.lo = gs - gs.li
      if(PARS$eval.on.subset) {
        rows <- apply(gs.lo, 1, sum) > 0
        cols <- apply(gs.lo, 2, sum) > 0
      }
      #Now we need to exclude all predictions that are in the prior:
      gs.li.inv = abs(1 - gs.li)
      comb.confs = comb.confs * gs.li.inv
      ord.idx = order(comb.confs[rows,cols], decreasing=T)
      PR = ChristophsPR(ord.idx, gs.lo[rows,cols])
      prec = PR$prec
      rec = PR$rec
      plot(rec[seq(1, length(rec), length.out=length.out)], prec[seq(1, length(rec), length.out=length.out)], pch=".", main = c("AUPR=", PR$auc))#and same input with "lines"
      lines(rec[seq(1, length(rec), length.out=length.out)], prec[seq(1, length(rec), length.out=length.out)])
    }
  } 
}
      
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  in.dir <- args[1]
}
if (length(args) == 2) {
  in.dir <- args[1]
  out.dir <- args[2]
}
cat(in.dir, out.dir)

plot.auprs(in.dir, out.dir)

#Usage:
#for d in output/bsubtilis/*/; do Rscript R_scripts/AUPRcurves.R $d AUPRs/; done
