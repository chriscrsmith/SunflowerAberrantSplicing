
# quick filter to remove genes with no splicing after consolidating alleles from same isoform 
# example: python post_consolidate_filter.py good_isos.txt consolidatedIsos.txt > good_isos_filtered.txt


import sys



### read in consolidated isoforms
consolidated = {}
with open(sys.argv[2]) as infile:
    for line in infile:
        newline = line.strip().split()
        iso = newline[0]
        old_gene = "_".join(iso.split("_")[0:4])
        clusters = newline[1].split("|")
        if old_gene not in consolidated:
            consolidated[old_gene] = list(clusters)



### read in the "good_isos" list
with open(sys.argv[1]) as infile:
    for line in infile:
        new_gene, isos = line.strip().split()
        isos = isos.split(",")
        old_gene = new_gene.split(".")[0]
        if old_gene in consolidated:
            keep = True # (default)
            for cluster in consolidated[old_gene]:
                if len(isos) <= len(cluster.split(",")):
                    keep = False
                    sys.stderr.write("\tfiltered " + old_gene + "\n")
            if keep == True:
                print(line.strip())
        else:
            print(line.strip())    
