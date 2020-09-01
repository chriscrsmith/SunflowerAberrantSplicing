

import sys
import numpy
from copy import deepcopy
import time
import rpy2.robjects as ro
import random





# parameters
pheno_number = int(sys.argv[3])-1 # -1 for 0-indexing
if pheno_number < 0:
    sys.stderr.write("pheno number needs to be >0\n")
    sys.exit()
locus1_number = sys.argv[4] # needs to be an index >0, or "all"
if locus1_number != "all":
    if int(locus1_number) < 1:
        sys.stderr.write("locus number needs to be >0\n")
        sys.exit()
permute = sys.argv[5]  
if permute != "yes" and permute != "no" and permute != "stratified" and permute != "count":
    sys.stderr.write("arg5 needs to be yes or no or stratified or count\n")
    sys.exit()
if locus1_number == "all" and permute == "stratified":
    sys.stderr.write("you don't want to do stratified permutations while doing a whole genome two-locus scan\n")
    sys.exit()
min_samples = 3





# read in genotypes
with open(sys.argv[1]) as infile:
    header_geno = infile.readline().strip().split()
    chroms = infile.readline().strip().split()
    positions = infile.readline().strip().split()
    num_genos = len(chroms)
    genotype_data = [ [ ] for i in range(num_genos)] # initialize giant data matrix
    RIL_order_genos = [] # keeping track of the correct RIL order
    for line in infile:
        newline = line.strip().split()
        rilID, genos = newline[0], newline[1:]
        RIL_order_genos.append(rilID)
        for g in range(num_genos):
            datum = str(genos[g])
            if datum == "HA89":
                datum = 0
            elif datum == "1238":
                datum = 2
            elif datum == "hetero":
                datum = 1
            genotype_data[g].append(datum)





            
# read in phenotypes
with open(sys.argv[2]) as infile:
    header_pheno = infile.readline().strip().split()
    num_phenos = len(header_pheno)-1
    phenotype_data = [ [ ] for i in range(num_phenos) ]
    RIL_order_phenos = []
    for line in infile:
        newline = line.strip().split()
        rilID, phenos = newline[0], newline[1:]
        RIL_order_phenos.append(rilID)
        for p in range(num_phenos):
            phenotype_data[p].append(phenos[p])




            
# reorder the phenotypes to get the RILS in the same order
if permute == "no" or permute == "count":
    temp_data = list(phenotype_data) # copying the phenotype data
    phenotype_data = [ [ ] for i in range(num_phenos) ] # new empty matrix
    reorder_indices = {}
    for ril in range(100): # looking through each ril number- in the order of the RILs in the genotype file
        ril_index = RIL_order_phenos.index(RIL_order_genos[ril]) # this is the index of the RIL in the phenotype (presumably different)
        reorder_indices[ril] = ril_index    
    for pheno in range(num_phenos):
        for ril in range(100):
            phenotype_data[pheno].append( temp_data[pheno][reorder_indices[ril]] )




            
# permute phenotypes
if permute == "yes":
    old_mat = list(phenotype_data)
    phenotype_data = [ None for i in range(num_phenos) ] # new empty matrix
    for row in range(num_phenos):
        temp_row = list(old_mat[row])
        random.shuffle(temp_row)
        phenotype_data[row] = list(temp_row)




        
# stratified permutation
if permute == "stratified":
    geno1 = genotype_data[int(locus1_number)-1] # -1 for 0-indexing
    transposed_phenotype_data = numpy.array(phenotype_data).T.tolist() # transpose the phenotype matrix to get 100 rows.
    for group in range(3): # 3 strata. Within each genotype group (stratum), permute the phenotypes.
        temp_matrix = [] 
        for s in range(len(geno1)):
            if geno1[s] == group: # if this individual is in the current group, add to the temperary matrix
                temp_matrix.append( list(transposed_phenotype_data[s]) )
        random.shuffle(temp_matrix) # shuffle the matrix
        sample_counter = 0
        for s in range(len(geno1)): # loop back through the rows of the matrix, and add them back to the stratum positions in the original matrix
            if geno1[s] == group:
                transposed_phenotype_data[s] = list(temp_matrix[sample_counter])
                sample_counter += 1
    # unindent
    phenotype_data = numpy.array(transposed_phenotype_data).T.tolist()
    # (this approach didn't touch the samples with missing data. The lm() function should ignore/remove these samples completely.)




    
##### begin scan #####
pures = [0,3]
mixed_indices = [1,2,4,5,6,7,8] 
mixed_set_1 = [1,2,8]
mixed_set_2 = [4,5,7]
pures_and_mixed_set_1 = pures + mixed_set_1
pures_and_mixed_set_2 =	pures +	mixed_set_2
count_total_tests = 0
count_before_dimsum = 0
t0 = time.time()
min_p = 1
pheno = phenotype_data[pheno_number]
for p in range(len(pheno)):
    datum = pheno[p]
    if datum == "NA":
        pass
        #datum = None       
    else:
        datum = float(datum)
    # unindent
    pheno[p] = datum
geneID = header_pheno[pheno_number+1]
if locus1_number == "all":
    locus1_range = [0, num_genos-1]
else:
    locus1_range = [int(locus1_number)-1, int(locus1_number)]
for locus1 in range(locus1_range[0], locus1_range[1]):
    sys.stderr.write("locus1: " + str(locus1+1) + "\n")  
    chrom1 = chroms[locus1]
    pos1 = positions[locus1]
    geno1 = genotype_data[locus1]
    if locus1_number == "all":
        locus2_range = [locus1+1, num_genos] # if testing all pairs of snps, want to avoid double-checking 
    else:
        locus2_range = [0, num_genos] # if locus1 is a particular snp, testing it against every other locus in the genome. 
    for locus2 in range(locus2_range[0], locus2_range[1]):
        geno2 = genotype_data[locus2]
        chrom2 = chroms[locus2]
        if chrom1 != chrom2:


            


            # count each genotype group with non missing data
            genotype_counts = [0,0,0,0,0,0,0,0,0]
            phenotype_totals = [0,0,0,0,0,0,0,0,0]
            for ril_index in range(100):
                if pheno[ril_index] != "NA": # if non-missing data
                    if geno1[ril_index] == 0 and geno2[ril_index] == 0: # (0 -P)
                        genotype_counts[0] += 1
                        phenotype_totals[0] += pheno[ril_index]
                    elif geno1[ril_index] == 0 and geno2[ril_index] == 2: # (1 -CM)
                        genotype_counts[1] += 1
                        phenotype_totals[1] += pheno[ril_index]
                    elif geno1[ril_index] == 0 and geno2[ril_index] == 1: # (2)
                        genotype_counts[2] += 1
                        phenotype_totals[2] += pheno[ril_index]
                    elif geno1[ril_index] == 2 and geno2[ril_index] == 2: # (3 -P)
                        genotype_counts[3] += 1
                        phenotype_totals[3] += pheno[ril_index]
                    elif geno1[ril_index] == 2 and geno2[ril_index] == 0: # (4 -CM)
                        genotype_counts[4] += 1
                        phenotype_totals[4] += pheno[ril_index]
                    elif geno1[ril_index] == 2 and geno2[ril_index] == 1: # (5)
                        genotype_counts[5] += 1
                        phenotype_totals[5] += pheno[ril_index]
                    elif geno1[ril_index] == 1 and geno2[ril_index] == 1: # (6)
                        genotype_counts[6] += 1
                        phenotype_totals[6] += pheno[ril_index]
                    elif geno1[ril_index] == 1 and geno2[ril_index] == 0: # (7)
                        genotype_counts[7] += 1
                        phenotype_totals[7] += pheno[ril_index]
                    elif geno1[ril_index] == 1 and geno2[ril_index] == 2: # (8)
                        genotype_counts[8] += 1
                        phenotype_totals[8] += pheno[ril_index]





            # filter
            analyze, count_it = False, False # default
            most_transgressive = [0,0] # default
            if genotype_counts[0] >= min_samples and genotype_counts[3] >= min_samples:
                for group in mixed_indices: # for each mixed genotype
                    if genotype_counts[group] > 0: # >0 is enough to calculate mean
                        me = phenotype_totals[group] / genotype_counts[group] # mean
                        if me > most_transgressive[1]: # if mean for this group is the most transgressive so far, save.
                            most_transgressive = [genotype_counts[group], me]
                        if genotype_counts[group] >= min_samples:
                            count_it = True
                if most_transgressive[0] >= min_samples: # if the most transgressive has >3 samples, then do a test.
                    analyze = True # default
            if count_it == True:
                count_before_dimsum += 1




            # prep data for lm() 
            if analyze == True and permute != "count": # if simply counting the tests, skip running lm() in R.
                ro.globalenv['df'] = pheno + geno1 + geno2 # the order, here, is important (I've messed it up before)
                ro.r( 'df = matrix(as.numeric(unlist(df)), ncol = 3, nrow = 100)' )
                ro.r( "lm1 = summary(lm( df[,1] ~ df[,2]*df[,3] ))$coefficients" )
                dim_sum = int(numpy.array( ro.r( "dim(lm1)[1]" ) )) # haha





                # output
                if dim_sum == 4: # sometimes lm() doesn't output an interaction term for probably many different reasons
                    ro.r( "p.lm = lm1[4, 4]" )
                    p = float(numpy.array(ro.r("p.lm"))[0])
                    count_total_tests += 1
                    if permute == "yes" or permute == "stratified":
                        if p < min_p:
                            min_p = float(p)
                    else:
                        if (p < 0.05):
                            pos2 = positions[locus2]
                            outline = [geneID, 
                                       pheno_number+1, 
                                       locus1+1, 
                                       chrom1, 
                                       pos1, 
                                       locus2+1, 
                                       chrom2, 
                                       pos2, 
                                       p]
                            outline = "\t".join(map(str,outline))
                            print outline



                    
# unindent all the way
t1 = time.time()
if permute == "yes" or permute == "stratified":
    print pheno_number+1, locus1+1, min_p
elif permute == "count":
    print "total tests:", count_before_dimsum
elif permute == "no":
    print "total tests:", count_total_tests
sys.stderr.write(str(t1-t0) + " seconds\n")


