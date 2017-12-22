PARS$input.dir <- 'input/yeast'
date.time.str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

PARS$exp.mat.file <- NULL
PARS$meta.data.file <- NULL
PARS$priors.file <- 'gold-standard.tsv'
PARS$gold.standard.file <- 'gold-standard.tsv'
PARS$deg.rates.file <- NULL

PARS$job.seed <- 42
PARS$job.seed.prior <- 42

PARS$num.boots <- 50 
PARS$cores <- 16 

PARS$delT.max <- 60000
PARS$delT.min <- 1
PARS$tau <- NULL

PARS$scale.desres <- TRUE
PARS$avg.diffs <- FALSE
PARS$reduce.expr <- FALSE
PARS$load.expr.meta.degrates.from.output <- TRUE

PARS$perc.tp <- 100  # percent of true priors that will be used; can be vector
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
PARS$set.deg.rates <- 'taus' #Taus for every gene are set according to the file that is fed in as an argument

PARS$save.to.dir <- paste(c('affy_TFA', PARS$method, PARS$prior.weight, date.time.str), collapse='_')

