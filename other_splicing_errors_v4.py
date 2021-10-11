
# checking for isoforms in parents not expressed in hybrids


import sys


# read in parent data
parent_data = {}
with open(sys.argv[1]) as infile:
    for f in infile:
        with open(f.strip()) as tpm_file:
            tpm_file.readline() # header
            for line in tpm_file:
                datum = line.strip().split()
                isoid = datum[0]
                geneid = "_".join(isoid.split("_")[0:4])
                if geneid not in parent_data: # saving all genes this time
                    parent_data[geneid] = {}
                if isoid not in parent_data[geneid]:
                    parent_data[geneid][isoid] = []
                parent_data[geneid][isoid].append(float(datum[3]))

                    
# go back through and see which genes have substantial expression in all 6 parents for ALL isoforms.
parent_alternative_splicing = {}
for geneid in parent_data:
    if len(parent_data[geneid]) > 1:
        count_big_isoforms = 0
        for isoid in parent_data[geneid]: 
            if min(parent_data[geneid][isoid]) >= 1: # if the isoform is >1 TPM in all six parents
                count_big_isoforms += 1 # count it
        if count_big_isoforms == len(parent_data[geneid]): # if ALL isoforms have substantial expression in all parents
            if geneid not in parent_alternative_splicing:
                parent_alternative_splicing[geneid] = []
            for isoid in parent_data[geneid]:
                parent_alternative_splicing[geneid].append(isoid)
sys.stderr.write("In the parents, there are " + str(len(parent_alternative_splicing)) + " genes where all parents have substantial expression of all isoforms (>1 TPM)\n")

                




# read in RILs
ril_data = {}
with open(sys.argv[2]) as infile:
    for f in infile:
        sys.stderr.write(f)
        ril = f.split("/")[9]
        ril_data[ril] = {}
        with open(f.strip()) as tpm_file:
            tpm_file.readline() # header    
            for line in tpm_file:
                datum = line.strip().split()
                isoid = datum[0]
                geneid = "_".join(isoid.split("_")[0:4])
                if geneid in parent_alternative_splicing:
                    if geneid not in ril_data[ril]:
                        ril_data[ril][geneid] = {}
                    ril_data[ril][geneid][isoid] = float(datum[3])


                    

# go back through and check RILs for errors
splicing_errors = {}
for geneid in parent_alternative_splicing:
    for ril in ril_data:
        total_gene_exp = 0 
        for isoid in ril_data[ril][geneid]: # first, loop through all isos to get total gene exp.
            total_gene_exp += ril_data[ril][geneid][isoid]
        if total_gene_exp >= 1:
            error = False # default setting, before checking.
            for isoid in parent_alternative_splicing[geneid]: # second, loop through isoforms to see if error
                if ril_data[ril][geneid][isoid] == 0:
                    error = True
            if error == True:
                if geneid not in splicing_errors:
                    splicing_errors[geneid] = 0
                splicing_errors[geneid] += 1
    
                    
            
    
                    
sys.stderr.write("there were " + str(len(splicing_errors)) + " splicing errors\n" )
for thing in splicing_errors:
    print thing, len(parent_alternative_splicing[thing]), splicing_errors[thing]










        
