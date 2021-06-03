
# quick filter to remove genes with no splicing after consolidating alleles from same isoform 
# example: python post_consolidate_filter.py good_isos.txt consolidatedIsos.txt > good_isos_filtered.txt


import sys


### read in the "good_isos" list
good_isos = {}
with open(sys.argv[1]) as infile:
    for line in infile:
        new_gene, isos = line.strip().split()
        isos = isos.split(",")
        for iso in isos:
            old_gene = new_gene.split(".")[0]
            if old_gene not in good_isos: # some are repeated                                                                                                                          
                good_isos[old_gene] = {}
            good_isos[old_gene][new_gene] = [] # adding in the new gene                                                                                                                
            for iso in isos:
                good_isos[old_gene][new_gene].append(iso)

                

### read in consolidated isoforms
with open(sys.argv[2]) as infile:
    for line in infile:
        newline = line.strip().split()
        iso = newline[0]
        old_gene = "_".join(iso.split("_")[0:4])
        clusters = newline[1].split("|")
        keep = True # default
        for new_gene in good_isos[old_gene]:
            if iso in good_isos[old_gene][new_gene]:
                for cluster in clusters:
                    if len(cluster.split(",")) == len(good_isos[old_gene][new_gene]):
                        keep = False
                if keep == True:
                        print (line.strip())
                else:
                    sys.stderr.write("\tfiltered " + "_".join(new_gene.split("_")[0:4]) + "\n")

           
