#!/bin/bash

namedir_final=$1
nCondClusts=$2
jobname=jobs/yeast_make-final-Network.R

mkdir -p $namedir_final

cp $jobname $namedir_final/jobs


for i in `seq 1 $nCondClusts`;
do
Rscript inferelator.R $jobname $namedir_final/clusters-conditions/expr_clust$i.tsv $namedir_final/clusters-conditions/meta_clust$i.tsv $i $namedir_final/cluster-genes/taus-predicted_clust$i.tsv $namedir_final
done

cp errors-network_final.txt $namedir_final

exit 0;


