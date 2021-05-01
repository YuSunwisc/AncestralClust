### 2_1_AncestralClust.r --- This is the 2_1_AncestralClust.r alg for R. 

# Authors: Yu Sun 
# (Last updated 2021-04-28)
#
# DESCRIPTION: This script is running for AncestralClust algorithm in R. The input file
# and the output FASTA data are in ../data folder.
#
# INSTRUCTIONS: First, download correct version of R in your local computer.
# Then run this file in your R environment. 

## !!WE STRONGLY SUGGEST YOU RUN THIS .R FILE UNDER R CONSOLE. IF YOU WANT TO RUN UNDER ANACODE OR FROM TERMINAL, HERE ARE SOME USEFUL SUGGESTION:

################################################################################################################################################
# PREPARE FOR TERMINAL
##############################################
# If you are on MacOS, try to install udunits first. (i.e. `brew install udunits`)
# Install GDAL from : http://www.kyngchaos.com/software/frameworks/ (before install it, put `echo 'export PATH=/Library/Frameworks/GDAL.framework/Programs:$PATH' >> ~/.bash_profile
# ` in your teriminal)
# download and install xquartz from https://www.xquartz
# Install gfortran : https://gcc.gnu.org/wiki/GFortranBinaries
# You also need to install 'igraph', 'pkg-config', 'svn', 'gdal' ,if you haven't, try with 
# brew install igraph
# brew install pkg-config
# brew install svn
# brew install gdal


# PREPARE FOR R(R version 3.6.3 (2020-02-29))
##############################################
# Step 1. Install udunits2: `install.packages("udunits2")`

# If you see error message like
# ```
# checking for udunits2.h... yes
# checking for ut_read_xml in -ludunits2... no
# -----Error: libudunits2.a not found-----
#      If the udunits2 library is installed in a non-standard location,
#      use --configure-args='--with-udunits2-lib=/usr/local/lib' for example,
#      or --configure-args='--with-udunits2-include=/usr/include/udunits2'
#      replacing paths with appropriate values for your installation.
#      You can alternatively use the UDUNITS2_INCLUDE and UDUNITS2_LIB
#      environment variables.
#      If udunits2 is not installed, please install it.
#      It is required for this package.
# ERROR: configuration failed for package ‘udunits2’
# * removing ‘/Users/guestadmin/opt/anaconda3/lib/R/library/udunits2’
# ```
# That means udunits2 library is installed in a non-standard location. Inside your local udunits2 folder, ./include/udunits2.h and ./lib/libudunits2.a 
# the most important files for correct installation. Above error message has typo too(it should be "configure.args" instead of "configure-args"). If you try
# `install.packages("udunits2",configure.args= "--with-udunits2-include=/usr/local/Cellar/udunits/2.2.27.6_1/include")`
# or `install.packages("udunits2",configure.args= "--with-udunits2-lib=/usr/local/Cellar/udunits/2.2.27.6_1/lib")`
# you will still see the error message, since you will always miss one file. Instead, try
# `install.packages("udunits2",configure.args= c("--with-udunits2-include=/usr/local/Cellar/udunits/2.2.27.6_1/include","--with-udunits2-lib=/usr/local/Cellar/udunits/2.2.27.6_1/lib"))`
# with your correct local udunits2 path. For all package using udunits2, always put correct `configure.args`.

# i.e.
# Anaconda: /Users/guestadmin/opt/anaconda3/lib/R/library
# Termianl: /usr/local/lib/R/library
# udunits2: /usr/local/Cellar/udunits


# Step 2. Install other package: 
# `install.packages(c("KernSmooth","classInt","units","sf","deldir","Matrix","expm","nlme","ade4","ggplot2","ape","seqinr","spdep","vegan","reshape2","igraph","gaston","pegas"))`

# If you see following error:
# ```
# CHOLMOD/Check/cholmod_check.c:170:2: error: no member named 'print_function' in
#       'struct cholmod_common_struct'
#         PRINTVALUE (Xx [p]) ;
#         ^~~~~~~~~~~~~~~~~~~
# CHOLMOD/Check/cholmod_check.c:112:2: note: expanded from macro 'PRINTVALUE'
#         P4 (" %23.15e", value) ; \
#         ^~~~~~~~~~~~~~~~~~~~~~
# CHOLMOD/Check/cholmod_check.c:93:24: note: expanded from macro 'P4'
# #define P4(format,arg) PR(4,format,arg)
#                        ^~~~~~~~~~~~~~~~
# CHOLMOD/Check/cholmod_check.c:84:31: note: expanded from macro 'PR'
#     if (print >= i && Common->print_function != NULL) \
#                       ~~~~~~  ^

# ```
# That means you have issue instailling igraph. To solve this:
# `brew uninstall suite-sparse`
# and `brew uninstall --ignore-dependencies suite-sparse` (https://github.com/igraph/rigraph/issues/173)

################################################################################################################################################


### ---Begin script--- ###

## Step 1. Install package and Read data

install.packages("adegenet", dep=TRUE, repos = "http://cran.us.r-project.org")
install.packages('phangorn', dep=TRUE, repos = "http://cran.us.r-project.org")
# Knitr produces a R session, without a default cran mirror unless you specifically asked for one. We tend to forget we need to set up CRAN for every R session when we use Rstudio because it takes care of it, but only for interactive use, not for knitr.
library("ape")
library("adegenet")
library("phangorn")
dna <- fasta2DNAbin(file="/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/initail_sequence_50_MSA.fasta")
# you may need to replace to the correct location of your data
D <- dist.dna(dna, model="raw")  # other model may have a lot of NAN, which we can not make NJ for the next step
# D <- as.matrix (D)
# dim(D)
# print(unname(D))
tre <- nj(D)

## Step 2. Prunes the longest edges from inital samples tree(i.e. 50 samples as report)


l<-subtrees(tre) # all the subtrees with root of all internal nodes, # of internal node = # of tips -1
initial_clust <- 5 # depend on your flavor or data
edgemax <-order(tre$ edge.length,decreasing=TRUE)[1:(initial_clust-1)] #find the location of max edges
clust <- list() # clust for initial samples

for (i in 1:length(edgemax))
{
    edge_i <- edgemax[i]     # i-th max edge   
    nodemax <- tre$edge[edge_i,] # the nodes of the maxedge
    n <- length(tre$tip)        # # of the all tips

    if (nodemax[1]-n<1)     # which means the max edge connects to a tip
    {
        l1 <- l[[nodemax[2]-n]]$tip
    } else if (nodemax[2]-n<1)  # which means the max edge connects to a tip
    {
        l1 <- l[[nodemax[1]-n]]$tip
    } else
    {
        N1 <- length(l[[nodemax[1]-n]]$tip) 
        N2 <- length(l[[nodemax[2]-n]]$tip)  # # of tips of two ends of the edge
        if (N1>N2)               # find out the node as root of subtree
        {
            nodemax_root <- nodemax[2]
        } else {
            nodemax_root <- nodemax[1]
        }
        l1 <- l[[nodemax_root-n]]$tip # all tips for the subtree
    }

    
    for (j in 1:i)
    {
        if (length(clust)==0 ) # at the first step, clust is an empty list
        {

            L <- tre $tip   # all the tips
            res <- L[is.na(pmatch(L,l1))] # the completion of l1 in L
            clust[[1]] <- l1
            clust[[2]] <- res # add a new member to the clust
            break
        } else 
        {
            if (all(l1 %in% clust[[j]]))
            {
                L <- clust[[j]]    # the union of the partition pd leaves
                res <- L[is.na(pmatch(L,l1))] # the completion of l1 in L
                clust[[j]] <- l1    # replace the older clust member
                clust[[i+1]] <- res # add a new member to the clust
                break
            }
        }
    }
}

clust

# print(tre)
# tre <- ladderize(tre)
# plot(tre, cex=.6)
# title("NJ tree reconstruction for initial 50 samples")


## Step 3. Pick centroid and write into centroid.txt and clust_centroid.txt

prefix = "/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/initial_clust"
clust_txt_suffix = "txt"
clust_centroid_txt_suffix = "centroid.txt"

for (i in 1:initial_clust)
{
    ## find centroid 
    midl <- floor((1+length(clust[[i]]))/2)
    centroid <- clust[[i]][midl]          # we want the centroid to be the one in the middle

    ## write clust.txt
    clust_txt_filename=paste(i, clust_txt_suffix, sep = ".")
    clust_txt_filename=paste(prefix, clust_txt_filename , sep = "_")
    write.table(clust[[i]], clust_txt_filename, append = TRUE, sep = " ", quote = FALSE, dec = ".", row.names = FALSE, col.names = FALSE) # we want to append this txt for new samples


    ## write clust_centroid.txt
    clust_centroid_txt_filename=paste(i, clust_centroid_txt_suffix, sep = "_")
    clust_centroid_txt_filename=paste(prefix, clust_centroid_txt_filename , sep = "_")
    write.table(centroid, clust_centroid_txt_filename, append = FALSE, sep = " ", quote = FALSE, dec = ".", row.names = FALSE, col.names = FALSE) 
}





## Step 4. Using python code write pairwisc fasta file with new samples

install.packages("devtool", repos = "http://cran.us.r-project.org")
install.packages("usethis", repos = "http://cran.us.r-project.org")
install.packages("reticulate", repos = "http://cran.us.r-project.org")
library(devtool)
library(usethis)
library(reticulate)   # to recall .py code
# reticulate::py_config() 
## use above command to find all py environment in your system, if the default one is not the desired one,
## try with command `use_python(/usr/bin/python3)` with your desired path. If it's still not working, try 
## to direct to you home R folder(you can use terminal command "$ which R" to find it out). Then create a 
## `.Renviron` file inside the folder, with one line i.e. `RETICULATE_PYTHON="/Library/Frameworks/Python.framework/Versions/3.9/bin/python3"` 
## with your descired py path. Restart you R and use command `reticulate::py_config()`, you should be able to see
# # i.e.
# # ```
# # python:         /Library/Frameworks/Python.framework/Versions/3.9/bin/python3
# # libpython:      /Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/config-3.9-darwin/libpython3.9.dylib
# # pythonhome:     /Library/Frameworks/Python.framework/Versions/3.9:/Library/Frameworks/Python.framework/Versions/3.9
# # version:        3.9.4 (v3.9.4:1f2e3088f3, Apr  4 2021, 12:32:44)  [Clang 6.0 (clang-600.0.57)]
# # numpy:          /Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/numpy
# # numpy_version:  1.20.2

# # NOTE: Python version was forced by RETICULATE_PYTHON
# # ```
# # ref:  https://stackoverflow.com/questions/50145643/unable-to-change-python-path-in-reticulate-r

py_run_file("/Users/shuqi/Desktop/563/Phylo563_FinalProject/scripts/2_2_write_new_fasta_.py")               # write pairwise fasta file with all centroids



## Step 5. Find the centroid with shortest distance with new samples


for (i in 1:3000)   # sample
{
    D <- list()
    new_member <- list()
    for (j in 1:5)   # clust
    {
        prefix = "/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/pairwise\ with\ centroids"
        middle = "new_sample"
        suffix = "fasta"

        filename = paste(i, suffix, sep = ".")
        filename = paste(middle, filename, sep = "_")
        filename = paste(j, filename, sep = "/")
        filename = paste(prefix, filename, sep = "/")

        dna <- fasta2DNAbin(file=filename)
        d <- dist.dna(dna, model="raw")
        D[j] <- d
        new_member[[j]] <- labels(d)
    }
    clust_num <- which.min(unlist(D))  # find the one with minimum 
    ## write labels into the cluster with smallest distance with centroid
    prefix = "/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/initial_clust"
    clust_txt_suffix = "txt"
    clust_txt_filename=paste(clust_num, clust_txt_suffix, sep = ".")
    clust_txt_filename=paste(prefix, clust_txt_filename , sep = "_")
    write.table(new_member[[clust_num]], clust_txt_filename, append = TRUE, sep = " ", quote = FALSE, dec = ".", row.names = FALSE, col.names = FALSE) 

}

## Step 6. Delete redundant members

py_run_file("/Users/shuqi/Desktop/563/Phylo563_FinalProject/scripts/2_3_final_clust.py")               # generate unrepeated cluster members





### ---End script--- ###







#######################################################  Appendix  ####################################################################################


# Example of plotting subtrees
##############################################

### Random tree with 6 leaves
# phy<-rtree(6)
# par(mfrow=c(1,1))  #plot 1 by 1 picture
# plot(phy, sub="Complete tree")

# ### Extract the subtrees
# l<-subtrees(phy)

# ### plot all the subtrees
# par(mfrow = c(2, 3))  #plot 2 by 3 picture
# # Create the loop.vector (all the columns)
# loop.vector <- 1:5
# for (i in loop.vector) plot(l[[i]], sub=paste("Node", l[[i]]$node.label[1]))
##############################################



# Potential Pasimony method not working on laptop
##############################################

### ---Begin script--- ###

# dna <- fasta2DNAbin(file="/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/sample_6.fasta")
# dna2 <- as.phyDat(dna)
# tre.ini <- nj(dist.dna(dna,model="raw"))
# parsimony(tre.ini, dna2)
# tre.pars <- optim.parsimony(tre.ini, dna2)
# plot(tre.pars, cex=0.6)

### ---End script--- ###


# It will cause error like below:
# ```
# *** caught segfault ***
# address 0x7f86d4563a44, cause 'memory not mapped'

# Traceback:
#  1: new_CppObject_xp(fields$.module, fields$.pointer, ...)
#  2: Rcpp::cpp_object_initializer(.self, .refClassDef, ...)
#  3: .Object$initialize(...)
#  4: initialize(value, ...)
#  5: initialize(value, ...)
#  6: new(Fitch, obj, as.integer(first_1), as.integer(m))
#  7: init_fitch(data, FALSE, FALSE, m = 2L)
#  8: fitch(tree, data, site = site)
#  9: parsimony(tre.ini, dna2)

# Possible actions:
# 1: abort (with core dump, if enabled)
# 2: normal R exit
# 3: exit R without saving workspace
# 4: exit R saving workspace
# ```
# The main reason I guess is due to memory not enough.
##############################################