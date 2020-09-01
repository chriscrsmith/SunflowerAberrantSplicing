



import sys, decimal
from copy import deepcopy

maxDistance = float(sys.argv[2]) # in cM



with open(sys.argv[1]) as infile:

    # printing, and saving, header line data
    header = infile.readline().strip()
    print header
    chroms = infile.readline().strip("\n")
    print chroms
    chroms = chroms.split()
    cM_positions = infile.readline().strip("\n")
    print cM_positions
    cM_positions = map( decimal.Decimal, cM_positions.split() )
    
    # initialize a dictionary of 17 chromsomes
    genotypes = {}
    for field in range(len(chroms)): # for the number of fields in the header line
        if chroms[field] not in genotypes:
            genotypes[chroms[field]] = {}
        genotypes[chroms[field]][cM_positions[field]] = None
        if "cp" in chroms[field]:
            sys.stderr.write("\n\nwith cM=0.5 for each SNP, the code can't differentiate between cpDNA SNPs\n\n")
            1/0
            
    # save sorted list of snp locations on each chrom
    locations = {}
    for chrom in genotypes:
        locs = genotypes[chrom].keys()
        locs.sort()
        locations[chrom] = list(locs)
        
    # loop through each sample, impute missing genotypes
    for line in infile:
        newline = line.strip().split()
        rilID, genos = newline[0], newline[1:]
        sys.stderr.write(rilID + "\n")
        currentGenos = deepcopy(genotypes) # make a copy of the empty dictionary w/ 17 chromosomes
        for col in range(len(genos)): # assigning the genos to dictionary with positions as keys
                                      # this is necessary because the snps are totally out of order, but we want them in order
            currentGenos[chroms[col]][cM_positions[col]] = genos[col]
            
        ##### the meat: slide along each chromosome and impute #####
        for lg in currentGenos:            
            for snp in range(1, len(locations[lg])-1 ): # skipping the first, and last, snps; don't want to impute them
                currentPos = locations[lg][snp]
                genotype = currentGenos[lg][currentPos]
                beforePos = locations[lg][snp-1]
                afterPos = locations[lg][snp+1]
                if (currentPos-beforePos <= maxDistance) and (afterPos-currentPos <= maxDistance): # two neighboring snps on either side
                    beforeGeno = currentGenos[lg][beforePos]
                    afterGeno = currentGenos[lg][afterPos]
                    if (beforeGeno == afterGeno) and (beforeGeno != "NA"): # the two neighbor snps are the same genotype (and not NA)
                        print "opp"
                        if genotype == "NA": # if missing data, fill it in
                            print "impute"
                            currentGenos[lg][currentPos] = str(beforeGeno)
                        elif beforeGeno == "hetero": # if the neighbor snps are heterozygous, make the middle heterozygous
                            if genotype != "hetero":
                                print "dropin"
                                currentGenos[lg][currentPos] = "hetero"

        # print
        outline = [rilID]
        for i in range(len(chroms)):
            outline.append(currentGenos[chroms[i]][cM_positions[i]])
        print "\t".join(outline)













        


        
