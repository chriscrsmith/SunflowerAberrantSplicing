


import sys, os
from operator import itemgetter
import numpy as np
from skbio.stats.composition import ilr
from statsmodels.multivariate.manova import MANOVA
from scipy import stats

trans_path = sys.argv[1]
goodisos_path = sys.argv[2]
consolidation_path = sys.argv[3]
parent_list = sys.argv[4]



##### read in list of transgressive isoforms #####
list_of_trans_isos = {}
list_of_trans_genes = {}
with open(trans_path) as infile:
    for line in infile:
        newline = line.strip().split()
        myiso = newline[0]
        mygene = "_".join(myiso.split("_")[0:4])
        list_of_trans_isos[myiso] = 0
        if mygene not in list_of_trans_genes:
            list_of_trans_genes[mygene] = []
        # unindent
        list_of_trans_genes[mygene].append( myiso )








##### read in list of good isoforms #####   
good_isos = {}
with open(goodisos_path) as infile:
    for line in infile:
        new_gene, isos = line.strip().split()
        isos = isos.split(",")
        for iso in isos:
            if iso in list_of_trans_isos:
                old_gene = new_gene.split(".")[0]
                if old_gene not in good_isos: # some are repeated
                    good_isos[old_gene] = {}
                good_isos[old_gene][new_gene] = [] # adding in the new gene
                for iso in isos:
                    good_isos[old_gene][new_gene].append(iso)




                    

##### read in list of consolidated alleles #####
consolidation = {}
with open(consolidation_path) as infile:
    for line in infile:
        trans_ID, others = line.strip().split()
        consolidation["_".join(trans_ID.split("_")[0:4])] = []
        for other in others.split("|"):
            for indiv in other.split(","):
                consolidation[indiv] = other.split(",")



            


##### read in parent expression data #####
sys.stderr.write("reading in parent data\n")
parents = {} # dictionary of parents, genes within parents, isoform
with open(parent_list) as list_file:
    for line in list_file:
        parent_path = line.strip()
        parent_id = parent_path.split("/")[-2] 
        parents[parent_id] = {}
        with open(line.strip()) as parent_file:
            parent_file.readline() # header
            for line in parent_file:
                newline = line.strip().split() 
                iso = newline[0]
                old_gene = "_".join(iso.split("_")[0:4])
                if old_gene in list_of_trans_genes:
                    tpm = float(newline[3])
                    if old_gene not in parents[parent_id]:
                        parents[parent_id][old_gene] = {}
                    # unindent
                    parents[parent_id][old_gene][iso] = tpm



            


##### now the business #####
for trans_iso in list_of_trans_isos:
    for old_gene in good_isos:
        for new_gene in good_isos[old_gene]:
            if trans_iso in good_isos[old_gene][new_gene]: # finding the correct "new" gene

                
                
                ##### calculate totals #####
                totals = [0,0,0,0,0,0]
                parent_ids = parents.keys() # ['HA89e_out', '1238e_out', 'HA89A_out', 'HA89b_out', '1238A_out', '1238b_out']
                iso_list = good_isos[old_gene][new_gene]
                for parent in range(6):
                    for iso in iso_list:
                        totals[parent] += parents[parent_ids[parent]][old_gene][iso] + 0.000001 # small value added to avoid zeros


                ##### calculate proportions #####
                isoform_clusters = []
                for iso in iso_list:
                    if iso in consolidation:
                        if consolidation[iso] not in isoform_clusters:
                            isoform_clusters.append(consolidation[iso])
                    else:
                        isoform_clusters.append([iso])
                props = [ [0]*len(isoform_clusters) for i in range(6)]
                for parent in range(6):
                    indiv_props = []
                    cluster_totals = [0]*len(isoform_clusters)
                    for cluster in range(len(isoform_clusters)):
                        for allele in isoform_clusters[cluster]: # sum up the alleles of the same splice form
                            cluster_totals[cluster] += parents[parent_ids[parent]][old_gene][allele] + 0.000001
                        # unindent
                        props[parent][cluster] =  cluster_totals[cluster] / totals[parent]  # then divide by the total gene expression for that parent
                count_isos_in_parents = 0 # real quick
                for cluster in range(len(isoform_clusters)):
                    found = False
                    for parent in range(6):
                        if props[parent][cluster] > 1e-05:
                            found = True
                    if found == True:
                        count_isos_in_parents += 1
                print "allele clusters found in parents:", count_isos_in_parents, trans_iso
                            
                        


                    
                ##### ILR transform #####
                ilrs = [None]*6
                for parent in range(6):
                    ilrs[parent] = ilr(props[parent]) # that was easy!  
                ilrs = np.array(ilrs)

                ##### test #####
                print
                print
                print trans_iso
                if len(ilrs.shape) > 1: # if more than one ILR column, use manova
                    pops = np.array([1,0,1,1,0,0]) # ['HA89e_out', '1238e_out', 'HA89A_out', 'HA89b_out', '1238A_out', '1238b_out']
                    manova = MANOVA(endog=ilrs, exog=pops)
                    print manova.mv_test()
                else: # if just one ILR column, usa t.test / anova
                    t = stats.ttest_ind([ilrs[0], ilrs[2], ilrs[3]], [ilrs[1], ilrs[4], ilrs[5]])
                    print t






