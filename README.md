# Read this first

This GitHub repository contains all of the code necessary to run InfereCLaDR (Tchourine et al., 2016: http://www.biorxiv.org/content/early/2017/01/31/104885)

All necessary files are inside the input directory.
For yeast (*S. cerevisiae*), the expression data file is too large and cannot be uploaded onto GitHub.
We ask the user to download the expression data file from https://drive.google.com/file/d/1fL6AnF4fJ0JVPebUgRpF1D_bp3oqunrC/view?usp=sharing
 - please name it `"expression.tsv"`, and put it in the `input/yeast/` folder.

# Requirements

The software is run on a Unix system via bash shell scripts that call R scripts and schedule jobs using PBS (SLURM option will be added soon)

The following software is required to run the InfereCLaDR:
 - R version 3.1.2 or higher
   - R packages: `corpcor`, `Matrix`, `inline`, `methods`, `elasticnet`, `lars`, `parallel`, `reshape2`, and `gplots`
 - The PBS job scheduler
   - soon a version will be added where SLURM can be used instead of PBS

# Running the InfereCLaDR

To run the InfereCLaDR, simply go into the highest-level directory (i.e. where you saved the InfereCLaDR code), and run

```
bash InfereCLaDR $namedir $nCondClusts $nGeneClusts $clustNamesPaper > errors-infereCLaDR.txt
```
where
 - `$namedir` is an argument that represents the name of the output directory;
 - `$nCondClusts` and `$nGeneClusts` represent the number of gene clusters and number of condition clusters, respectively; and
 - `$clustNamesPaper` is a boolean variable. If `TRUE`, the same naming scheme for condition and gene clusters will be used as in the InfereCLaDR paper, if `FALSE`, no names will be used. `TRUE` will only work if `$nCondClusts=4` and `$nGeneClusts=5`.

For example, to recapitulate the results of the InfereCLaDR paper, use the following command:

```
bash InfereCLaDR yeast/paper 4 5 TRUE > errors-infereCLaDR.txt
```

Note that if your machine or cluster that has PBS (or SLURM - coming soon) cannot perform computations not submitted via PBS (i.e. if it's a login node that only exists for the purpose of submitting jobs to PBS), you can run:
 - lines 3-17 (from `InfereCLaDR.sh`) on any machine (i.e. without PBS)
   - setting `$namedir`, `$nCondClusts`, `$nGeneClusts`, and `$clustNamesPaper` inside the file
 - copy and paste the code and the new output folder onto the machine/cluster that has PBS
 - run line 21 (i.e. `bash half-life_fitting/run_pbs.sh output/$namedir $nCondClusts > errors-DR_fitting.txt`) on that machine/system
 - copy the output from the machine/cluster that has PBS back onto the first machine, and
 - run the rest of the lines in `InfereCLaDR.sh` on the first machine.
