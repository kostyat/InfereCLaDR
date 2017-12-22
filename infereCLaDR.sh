#!/bin/bash

namedir=yeast/example
jobname=jobs/yeast_affy_tauNULL.R
nCondClusts=4 #number of condition clusters
nGeneClusts=5 #number of gene clusters
clustNamesPaper=TRUE #if TRUE, the same names are given to clusters as in the InfereCLaDR paper, if FALSE, no gene and condition cluster names are given

mkdir -p output/$namedir/

cp inferelator.R output/$namedir/
mkdir output/$namedir/jobs
cp $jobname output/$namedir/jobs
cp -r R_scripts output/$namedir/
cp infer-final-network.sh output/$namedir/

#cluster the expression matrix and split the metadata table accordingly:
Rscript R_scripts/cluster-create.R $nCondClusts $nGeneClusts $namedir

#run the Inferelator for the set of RNA half-lives (specified in half-life_fitting/inferelator_DT.pbs)
### THIS is the only line that requires PBS!
bash half-life_fitting/run_pbs.sh output/$namedir $nCondClusts > errors-DR_fitting.txt
###

#calculate the optimal RNA half-life for every bi-cluster:
Rscript half-life_fitting/auprs-by-cond-and-gene-clust.R output/$namedir/half-life_fitting $nCondClusts $nGeneClusts $clustNamesPaper
Rscript half-life_fitting/predicted-taus-all.R output/$namedir $nCondClusts $nGeneClusts $clustNamesPaper

#Now run the Inferelator using the full prior (Gold Standard in this case) and the fitted RNA half-lives to produce the best possible network
bash infer-final-network.sh output/$namedir $nCondClusts > errors-network_final.txt

#Combine the interaction confidence scores across the condition clusters
Rscript combine-cluster_networks.R output/$namedir $nCondClusts

cp errors-infereCLaDR.txt output/$namedir/

exit 0;
