


import sys
import rpy2.robjects as ro

splicingPath = sys.argv[1]
transgressive_list_path = sys.argv[2]
ril_list_path = sys.argv[3]





# read in the splicing data; all we need for this is a dictionary of old* isoform ids, and the associated new* geneIDs
splicing = {}
with open(splicingPath) as infile:
    for line in infile:
        newline = line.strip().split()
        new_gene_id = newline[0]
        isos = newline[1].split(",")
        for iso in isos:
            if iso in splicing:
                1/0 # (surprising)
            else:
                splicing[iso] = new_gene_id # dictionary of verified isoforms (regardless whether they are transgressive)

                

            
# read in list of transgressive transcripts
transgressive_list = {}
alternate_alleles = {}
with open(transgressive_list_path) as infile:
    for line in infile: 
        newline = line.strip().split()
        new_iso_id = newline[0]
        gene = splicing[new_iso_id]
        if gene not in transgressive_list:
            transgressive_list[gene] = []
        transgressive_list[gene].append(new_iso_id) # dictionary of GENES, with values being their corresponding transgressive isoform(s) 
        if newline[5] != "NA": # consolidating similar alleles
            alternate_alleles[new_iso_id] = newline[5].split(",")
            

        


# read in RIL data; need to organize this using the updated gene and isoform assignments
rildata = {}
ril_order = [] # keeping the order straight
with open(ril_list_path) as listfile:
    for filepath in listfile:
        ril = "_".join(filepath.split("/")[-2].split("_")[:-1]) 
        sys.stderr.write(ril+"\n")
        rildata[ril] = {}
        ril_order.append(ril)
        with open(filepath.strip()) as rilfile:
            rilfile.readline() # header
            for line in rilfile:
                newline = line.strip().split()
                isoform, tpm = newline[0], newline[3]
                if isoform in splicing: 
                    gene = splicing[isoform]
                    if gene in transgressive_list: # if the gene the isoform came from is transgressive, keep track of this isoform
                        if gene not in rildata[ril]:
                            rildata[ril][gene] = {}
                        rildata[ril][gene][isoform] = float(tpm) # + 0.01



                        
# print some headers real quick
outline = ["RIL"]
for gene in rildata[ril]:
    for trans in transgressive_list[gene]:
        outline.append(trans)
print "\t".join(outline)



                    
# go back through and find the proportion of each transgressive isoform
ro.r( "library(compositions)" ) # R library; doing ILR transforms in R
for ril in ril_order:
    outline = [ril]
    for gene in rildata[ril]: # these are the newly-assigned gene designations
        totalTPM = 0
        for iso in rildata[ril][gene]:
            totalTPM += rildata[ril][gene][iso] 
        for trans in transgressive_list[gene]: # (here is a loop; but, currently, there is only a single isoform per gene so it's not necessary)
            trans_expression = rildata[ril][gene][trans]
            #
            if trans in alternate_alleles:
                for consol_iso in alternate_alleles[trans]:
                    trans_expression += rildata[ril][gene][consol_iso] # now, we've counted all "alleles" with same splice type towards the transgressive expression value
            #
            if totalTPM < 1: # first of all, I don't want to consider the proportion transgressive, if it doesn't have substantial expression (>1tpm) at the gene level
                outline.append("NA")
            elif trans_expression == 0: # if zero transgressive, assign 0 (can't divide by 0)
                outline.append(str(0))
            elif trans_expression > 0:
                prop = trans_expression / totalTPM
                outline.append(str(prop))
            else:    
                1/0 # would be surprising
    print "\t".join( outline )






        
