# MainReadme: AncestralClust Alg
	Author: Yu Sun
    Last update 2021.04.30
    Github: https://github.com/YuSunwisc/Phylo563_FinalProject
    Contact: ysun258@wisc.edu
## 1.Introduction

This project is going to do a simulation about ancestral sequence reconstruction with less randomly selected samples. It's well known that clustering is a fundamental task in the analysis of nucleotide sequences, and traditional clustering methods have mostly focused on optimizing high speed clustering of highly similar sequences, but state-of-art clustering method starting with less samples may also led to a good result. There we're going to run an "homemade" version of AncestralClust alg with ssRNA data.

## 2. Content
This repository contains 3 folders [script](script), [data](data) and [graphs](graphs), with two report files [report.md](report.md) and [report.html](report.html). Please direct to [script](script) and run all the scripts by step by step instruction in [report.md](report.md). The final reslut of cluster members are in [data/final_clust](data/final_clust) folder. All detailed questions please check readme file in each folder.

## 3. About AncestralClust(from original paper)
- Step 1: With total number of N nucleotide sequences, choose a relatively small fraction of it(i.e. k=3%N), and use wavefront alforithm to do pairwise alignment.

- Step 2: Use distance method to construct the tree structure: here we use Jukes-Cantor model with neighbor-joining algorithm for our instance.

- Step 3: If we assume there is C clusters in N, then we prune C-1 longest branch to have C initial clusters.

- Step 4: For each initial cluster, use kalign3 algotithm to do multi sequence alignemnt.

- Step 5: Use MAP estimator to construct root of each initial cluster. 

- Step 6: For each unsampled nucleotide sequence, use wavefront alforithm again to get pairwise alignment with all C roots, and assign it to the cluster with shortest JC distance.

- Step 7: If the shortest distance to any of the C ancestral sequences is larger than the average distance between clusters, the sequence is saved for the next iteration.

- Step 8: In each iteration after the first iteration, a cut of a branch in the phylogenetic tree is chosen if the the branch is longer that the average length of branches cut in the first iteration.










