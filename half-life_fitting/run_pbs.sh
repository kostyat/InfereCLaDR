#!/bin/bash

namedir_DR_fitting=$1/half-life_fitting
clust_folder=$1/cluster-conditions

mkdir -p $namedir_DR_fitting/
cur_dir=$(pwd)

jobname=jobs/yeast_fit-DegRates.R
nCondClusts=$2

cp half-life_fitting/inferelator.R $namedir_DR_fitting/
cp half-life_fitting/inferelator_DT.pbs $namedir_DR_fitting/inferelator_DT.pbs
cp half-life_fitting/run_pbs.sh $namedir_DR_fitting

clustArray=(`seq 1 $nCondClusts`)

exprs_prepend=(${clustArray[@]/#/$clust_folder/expr_clust})
exprs=(${exprs_prepend[@]/%/.tsv})

metas_prepend=(${clustArray[@]/#/$clust_folder/meta_clust})
metas=(${metas_prepend[@]/%/.tsv})

kMax=`expr ${#exprs[@]} - 1`

mkdir $namedir_DR_fitting/errors-DR_fitting/

for j in `seq 42 61`;
do
  for k in `seq 0 $kMax`;
  do
    qsub -v JOB=$jobname,EXPR=${exprs[$k]},META=${metas[$k]},CLUSTN=$[$k+1],PRSEED=$j,DIR=$namedir_DR_fitting,CUR_DIR=$cur_dir -N inferelator_prior${j[0]}_clust4_$[$k+1] half-life_fitting/inferelator_DT.pbs
  done
done


cp errors-DR_fitting.txt $namedir_DR_fitting/

exit 0;
