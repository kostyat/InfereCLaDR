PARS$input.dir <- 'input/yeast'
date.time.str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

PARS$exp.mat.file <- NULL
PARS$meta.data.file <- NULL
PARS$priors.file <- 'gold-standard.tsv'
PARS$gold.standard.file <- 'gold-standard.tsv'
PARS$tf.names.file <- 'tf_names.tsv'

PARS$job.seed <- NULL
PARS$job.seed.prior <- NULL

PARS$num.boots <- 50 
PARS$cores <- detectCores()/2 

PARS$delT.max <- 60000
PARS$delT.min <- 1
PARS$tau <- NULL #RNA half-life is NULL because it's set in the inferelator_DT.pbs file

PARS$scale.desres <- TRUE
PARS$avg.diffs <- FALSE
PARS$reduce.expr <- TRUE
PARS$load.expr.meta.degrates.from.output <- TRUE

PARS$perc.tp <- 50  # percent of true priors that will be used; can be vector
PARS$perm.tp <- 1  # number of permutations of true priors
PARS$perc.fp <- 0  # percent of false priors (100 = as many false priors as 
                   # there are true priors); can be vector
PARS$perm.fp <- 1  # number of permutations of false priors
PARS$pr.sel.mode <- 'random'  # prior selection mode: 'random' or 'tf'
                              # if 'random', the true priors are randomly chosen
                              # from all priors edges, if 'tf', 
                              # PARS$perc.tp is interpreted as the percent of
                              # TFs to use for true priors and all interactions
                              # for the chosen TFs will be used


PARS$eval.priors.on.subset <- FALSE
PARS$eval.on.subset <- TRUE

PARS$method <- 'BBSR'
PARS$prior.weight <- 1.1
PARS$dr.weight <- 0.1
PARS$use.tfa <- TRUE
PARS$deg.rates <- FALSE
PARS$set.deg.rates <- FALSE

#PARS$save.to.dir <- paste('output/yeast/gpl90', PARS$method, PARS$prior.weight, sep='_')
PARS$save.to.dir <- paste(c('affy_TFA', PARS$method, PARS$prior.weight, date.time.str), collapse='_')

