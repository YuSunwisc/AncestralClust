### 2_3_final_clust.py --- This is the script clean the redundant clust members

# Authors: Yu Sun 
# (Last updated 2021-04-30)
#
# DESCRIPTION: This is the script clean the redundant clust members.
#
# INSTRUCTIONS: Run in Py.


################################################################################################################################################


### ---Begin script--- ###

## Step 1. Reduced duplications by lines


for j in range(1,6):
    list01 = []
    for i in open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/initial_clust_"+str(j)+".txt"):
        if i in list01:
            continue
        list01.append(i)
    with open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/final_clust/final_clust_"+str(j)+".txt", 'w') as handle:
        handle.writelines(list01)






### ---End script--- ###


