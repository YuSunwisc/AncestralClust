### 1_3_muscle.sh --- muscle MSA 
### alignments drawn from consecutive blocks of a reference genome and SNP data
#
# Authors: Yu Sun 
# (Last updated 2021-04-26)
#
# DESCRIPTION: This script is running for MSA using muscle algorithm. The input file
# and the output FASTA data are in ../data folder.
#
# INSTRUCTIONS: First, download correct version of Muscle from 
# https://www.drive5.com/muscle/downloads.htm. Unzip and hange your download file 
# with name "muscle"(Please do not add any extension to it. If you're using MacOSX 
# 64 bits, then you can just run this script without any change). Move that muscle 
# file into current folder and replace the old one. 

# To run this script, navigate to current folder and run 

# "bash 1_3_muscle.sh $1"
# $1 is the FULL path of the input file 
# If you syetem require to move muscle to trash since it's unverified. Please change your privacy.
# If you see error like "no such file", change the variables "muscle" below 
# to be correct FULL paths for the binary file.

# Example:
# ```
# bash muscle.sh /Users/shuqi/Desktop/563/Phylo563_FinalProject/data/example_input.fasta 

# MUSCLE v3.8.31 by Robert C. Edgar

# http://www.drive5.com/muscle
# This software is donated to the public domain.
# Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

# example_MSA_input 2 seqs, max length 29903, avg  length 29877
# 00:00:00      2 MB(0%)  Iter   1  100.00%  K-mer dist pass 1
# 00:00:00      2 MB(0%)  Iter   1  100.00%  K-mer dist pass 2
# 00:00:40   970 MB(11%)  Iter   1  100.00%  Align node       
# 00:00:40   970 MB(11%)  Iter   1  100.00%  Root alignment
# ```


################################################################################################################################################



### ---Begin script---

## Step 1. Run muscle-MSA

muscle="/Users/shuqi/Desktop/563/Phylo563_FinalProject/scripts/muscle"
input="$1"
input_name="${1%.*}"  #${foo##*/} is greedy back trim until last /, and ${foo%/*} is greedy front trim until last /
output="$input_name"_MSA.fasta
touch $output


chmod +x $muscle
$muscle -in $input -out $output  # muscle command



### ---End script---


################################################################################################################################################
