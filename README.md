# Phylo563_FinalProject: Clustering with less samples(Yu Sun, 04.09.2021)
## 1.Introduction

This project is going to do a simulation about ancestral sequence reconstruction with less randomly selected samples. It's well known that clustering is a fundamental task in the analysis of nucleotide sequences, and traditional clustering methods have mostly focused on optimizing high speed clustering of highly similar sequences, but state-of-art clustering method starting with less samples may also led to a good result. But our recent analysis shows that, on some simplified models, there is a chance that this method could be wrong--some nucleotide sequences will be assigned into wrong groups. 
I'll focus on two aspects: for a fixed dataset, I will try to show when this method is working and when it gives a unreliable clusters.

## 2. AncestralClust
- Step 1: With total number of N nucleotide sequences, choose a relatively small fraction of it(i.e. k=3%N), and use wavefront alforithm to do pairwise alignment.

- Step 2: Use distance method to construct the tree structure: here we use Jukes-Cantor model with neighbor-joining algorithm for our instance.

- Step 3: If we assume there is C clusters in N, then we prune C-1 longest branch to have C initial clusters.

- Step 4: For each initial cluster, use kalign3 algotithm to do multi sequence alignemnt.

- Step 5: Use MAP estimator to construct root of each initial cluster. 

- Step 6: For each unsampled nucleotide sequence, use wavefront alforithm again to get pairwise alignment with all C roots, and assign it to the cluster with shortest JC distance.

- Step 7: If the shortest distance to any of the C ancestral sequences is larger than the average distance between clusters, the sequence is saved for the next iteration.

- Step 8: In each iteration after the first iteration, a cut of a branch in the phylogenetic tree is chosen if the the branch is longer that the average length of branches cut in the first iteration.

## 3. MAP estimator result
In our unpubilshed proof, we find some extreme cases that the unsampled sequence tends to have stronger correlation with the "wrong root", which leads to the wrong clustering. Here we're going to choose a fixed set of samples that gives, both good and bad clusters.

## 4. Dataset

Current dataset we choose comes from : 6, 18S, and Cytochrome Oxidase I (COI)) from the CALeDNA project Meyer.