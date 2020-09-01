

minRILs = 35


import sys, decimal



with open(sys.argv[1]) as infile:
    header = infile.readline().strip().split()[1:] # cutting off the first field, which is "RIL"
    chroms = infile.readline().strip().split()
    cM_positions = infile.readline().strip().split()
    num_snps = len(chroms)
    RIL_IDS = []
    
    # initialize data matrix
    genotypes = []
    for snp in range(num_snps):
        genotypes.append([])
        
    # loop through each sample and fill the matrix
    for line in infile:
        newline = line.strip().split()
        rilID, genos = newline[0], newline[1:]
        RIL_IDS.append(rilID)
        for geno in range(num_snps):
            genotypes[geno].append( genos[geno] )


            
# go back through and filter some snps
new_header = ["RIL"]
new_chroms = [""]
new_positions = [""]
new_genotypes = [ [ ] for i in range(100) ]
for field in range(num_snps):
    numMissing = genotypes[field].count("NA")
    num1238 = genotypes[field].count("1238")
    numHA89 = genotypes[field].count("HA89")
    if numMissing <= (100-minRILs):
        if num1238 >= 5:
            if numHA89 >= 5:
                new_header.append(header[field])
                new_chroms.append(chroms[field])
                new_positions.append(cM_positions[field])
                for ril in range(100):
                    new_genotypes[ril].append(genotypes[field][ril])


                
# print
print "\t".join(new_header)
print "\t".join(new_chroms)
print "\t".join(new_positions)
for ril in range(100):
    print "\t".join( [RIL_IDS[ril]] + new_genotypes[ril] )










        
