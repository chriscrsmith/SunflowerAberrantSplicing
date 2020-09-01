

import sys
import numpy as np
import numpy.ma as ma
from scipy import stats
fixedSNPs_path = sys.argv[1]
RILlist_path = sys.argv[2]
RIL_vcf_path = sys.argv[3]
minDepth = float(sys.argv[4])
nuc_list = ['A','T','C','G']




# read in the list of fixed snps with genomic (bp & cM) positions
fixedSNPs = {}
cmPositions = {}
with open(fixedSNPs_path) as infile:
    infile.readline() # header
    for line in infile:
        newline = line.strip().split()
        contig, pos, chrom, bp, cM, a1, a2 = newline[:]
        if contig not in fixedSNPs:
            fixedSNPs[contig] = {}
        fixedSNPs[contig][pos] = [chrom, bp] # this dictionary represents the translation from transcriptome position to genome position
        if chrom not in cmPositions:
            cmPositions[chrom] = {}
        cmPositions[chrom][bp] = cM


        
            
# read RIL list (100 RILs)
RILlist = []
with open(RILlist_path) as infile:
    for line in infile:
        newline = line.strip()
        if newline != "RIL": # header
            RILlist.append(newline)


            
            

# go through the RIL vcfs; filter genotypes with insufficient read depth
alleles = {}
numRILs = len(RILlist)
for contig in fixedSNPs: # pre-allocate the dictionary; also, converting to bp position
    for pos in fixedSNPs[contig]:
        chrom, bp = fixedSNPs[contig][pos][:]
        if chrom not in alleles:
            alleles[chrom] = {}
        alleles[chrom][bp] = [np.NaN]*numRILs
for ril in range(len(RILlist)):
    sys.stderr.write(RILlist[ril]+"\n")
    with open(RIL_vcf_path+"/"+RILlist[ril]+".vcf") as infile:
        for line in infile:
            if line[0] != "#":
                newline = line.strip().split()
                contig, pos = newline[0:2]
                a1,a2 = newline[3], newline[4]
                if contig in fixedSNPs:
                    if pos in fixedSNPs[contig]:
                        if a1 in nuc_list and a2 in nuc_list: # necessary check. Although parents fixed, the RILs can bring third allele.
                            RIL_data = newline[9].split(":")
                            RIL_DP = int(RIL_data[2]) 
                            if RIL_DP >= minDepth:
                                chrom, bp = fixedSNPs[contig][pos][:]
                                RIL_geno = RIL_data[0]
                                if RIL_geno == "0/1":
                                    alleles[chrom][bp][ril] = "hetero"
                                else:
                                    geno_1238 = newline[10].split(":")[0] 
                                    if RIL_geno == geno_1238:
                                        alleles[chrom][bp][ril] = "1238"
                                    else:
                                        alleles[chrom][bp][ril] = "HA89"




            
# print header
header, row2, row3 = ["RIL"], [""], [""]
for chrom in alleles:
    for pos in alleles[chrom]:
        cm_position = cmPositions[chrom][pos]
        binID = chrom + "_" + pos + "_" + cm_position
        header.append(binID)
        row2.append(chrom)
        row3.append(cm_position)
print "\t".join(header)
print "\t".join(row2)
print "\t".join(row3)
        




# print genotypes
numRILs = int(numRILs)
for ril in range(numRILs):
    newline = [RILlist[ril]]
    for chrom in alleles:
        for pos in alleles[chrom]:
            geno = str(alleles[chrom][pos][ril])
            if geno == "nan":
                geno = "NA"
            newline.append(geno)
    print "\t".join(newline)






            
