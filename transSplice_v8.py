

import sys, os
from operator import itemgetter

min_expression_next_closest_isoform = 1 # requiring each parent to have this much of the next-most-similar isoform to the transgressive isoform
min_exon_size = 26 # requiring an exon of this length to differentiate the transgressive and the next most similar isoform. 

goodisos_path = sys.argv[1]
parent_list = sys.argv[2] # list of paths to parent quant.sf files
RIL_list = sys.argv[3] # list of paths to RIL quant.sf files



##### read in list of good isoforms #####   
good_isos = {}
with open(goodisos_path) as infile:
    for line in infile:
        new_gene, isos = line.strip().split()
        old_gene = new_gene.split(".")[0]
        if old_gene not in good_isos: # some are repeated
            good_isos[old_gene] = {}
        good_isos[old_gene][new_gene] = [] # adding in the new gene
        for iso in isos.split(","):
            good_isos[old_gene][new_gene].append(iso)




##### read in the parent expression data #####   
sys.stderr.write("reading in parent data\n")
parents = {} # dictionary of parents, genes within parents, isoforms within genes
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
                if old_gene in good_isos:
                    tpm = float(newline[3])
                    for new_gene in good_isos[old_gene]: # looping through the new genes to find the one with this particular isoform
                        if iso in good_isos[old_gene][new_gene]: # the iso was associated with this new gene ID
                            if new_gene not in parents[parent_id]:
                                parents[parent_id][new_gene] = {}
                            parents[parent_id][new_gene][iso] = tpm




                            
                            
##### quick check: how many transcripts, and how many isoforms, are not found in the parents at all? #####   
sys.stderr.write("checking parent data\n")
genes_significant_expression_all_parents = {} # recording which GENES have >1 tpm in all six parents
isoform_significant_expression_all_parents = {}
zero_expression_in_parents = {} # recording which isoforms have absolutely zero expression in all six parents
for old_gene in good_isos: # for each (old) gene with verified splicing
    for new_gene in good_isos[old_gene]: # for each new gene associated with the old gene id

        # check total expression in each parent: >1 
        keep = True # default setting only: keep the gene. If we see that any parent has <1 tpm, then exclude the gene.
        for parent in parents:
            total_exp = 0 # counting the total expression for each parent
            for iso in good_isos[old_gene][new_gene]: # loop through each iso and add to total expression for this parent
                total_exp += parents[parent][new_gene][iso]
            if total_exp < 1: # if ANY parent has less than 1 total tpm, then exclude for all parents
                keep = False
        if keep:
            genes_significant_expression_all_parents[new_gene] = 0

        # check if individual isoforms have 0 expression in the parents
        for iso in good_isos[old_gene][new_gene]: # for EACH individual isoform, do any parents have any?
            transgressive = True # default setting only
            expressed = True 
            for parent in parents:
                if parents[parent][new_gene][iso] > 0: 
                    transgressive = False # if any parent has any reads from this isoform, not transgressive 
                if parents[parent][new_gene][iso] < min_expression_next_closest_isoform: # if any parent has <1 tpm, then we're counting it as not substantially expressed
                    expressed =	False
            if transgressive:
                zero_expression_in_parents[iso] = 0
            elif expressed:
                isoform_significant_expression_all_parents[iso] = 0
                







##### read in the rils #####      
sys.stderr.write("reading in rils\n")
trans_dict = {}
with open(RIL_list) as list_file:
    for line in list_file:
        sys.stderr.write(line.strip()+"\n")
        with open(line.strip()) as ril_file:
            ril_file.readline() # header
            for line in ril_file:
                newline = line.strip().split()
                iso = newline[0]
                old_gene = "_".join(iso.split("_")[0:4])
                if old_gene in good_isos: # remember not all genes in the Trinity output were verified by my preceding script
                    for new_gene in good_isos[old_gene]: # looping through the new genes to find the one with this particular isoform 
                        if iso in good_isos[old_gene][new_gene]: # this is the new gene associated with the current iso
                            if (new_gene in genes_significant_expression_all_parents) and (iso in zero_expression_in_parents): 
                                tpm = float(newline[3])
                                if tpm >= 1: 
                                    if iso not in trans_dict:
                                        trans_dict[iso] = 0
                                    trans_dict[iso] += 1





##### Next, I want to create a fasta file with all the (verified) isoforms and align them #####     
for trans_iso in trans_dict: 
    old_gene = "_".join(trans_iso.split("_")[0:4])
    for new_gene in good_isos[old_gene]: # again, looping through to find the correct new_gene
        if trans_iso in good_isos[old_gene][new_gene]: # this is the correct one
            os.system( "rm temp1.fa temp2.fa 2> stderr.txt" ) # removing these up front, just in case they were left over from before

            
            ###### EMBOSS needle (pairwise sequence alignment) #####
            if len(good_isos[old_gene][new_gene]) == 2:
                ids = list(good_isos[old_gene][new_gene])
                os.system( "grep -w -A 1 " + ids[0] + " /media/baswhd/ChrisData/Sunflowers/2020/Transcriptome/NewAssembly/Trinity_2019.fasta >> temp1.fa" )
                os.system( "grep -w -A 1 " + ids[1] + " /media/baswhd/ChrisData/Sunflowers/2020/Transcriptome/NewAssembly/Trinity_2019.fasta >> temp2.fa" )
                os.system( "needle temp1.fa temp2.fa -gapopen 10 -gapextend 0.5 -outfile temp1.needle 2> stderr.txt" )
                seqs = {} # initialize empty sequence dictionary
                with open("temp1.needle") as infile:
                    line_number = 0
                    seq1 = ""
                    seq2 = ""
                    for line in infile:
                        if line[0] != "#" and line != "\n" and "|" not in line:
                            newline = line.strip().split()
                            if newline != []: # this is necessary (filters out stretches of nonmatching positions)
                                line_number += 1
                                if line_number % 2 == 1: # odd numbered lines                                                                                               
                                    seq1 += line.split()[2]
                                else:                    # even numbered lines     
                                    seq2 += line.split()[2]
                seqs[ids[0]] = str(seq1)
                seqs[ids[1]] = str(seq2)
                os.system( "rm temp1.needle" ) # clean up   
            
                
            

            ##### MUSCLE (multiple sequence alignment) #####
            else:
                for iso in good_isos[old_gene][new_gene]:
                    os.system( "grep -w -A 1 " + iso + " /media/baswhd/ChrisData/Sunflowers/2020/Transcriptome/NewAssembly/Trinity_2019.fasta >> temp1.fa" )
                # unindent; only want to run muscle once per gene
                os.system( "muscle3.8.31_i86linux64 -in temp1.fa -out temp1.muscle 2> stderr.txt" )
                with open( "temp1.muscle" ) as infile:
                    seqs = {} # initialize empty sequence dictionary
                    seq = "" # initialize first sequence
                    for line in infile:
                        if line[0] == ">":
                            if seq != "": # if a sequence just ended and beginning a new sequence
                                seqs[iso_id] = seq
                                seq = ""
                            iso_id = line.strip().split()[0].split(">")[1]
                        else:
                            seq += line.strip()
                if seq != "":
                    seqs[iso_id] = seq # adding the last transcript in
                os.system( "rm temp1.muscle" ) # clean up



            ##### comparing isoform alignments #####
            trans_seq = seqs[trans_iso]
            iso_list = seqs.keys()
            iso_list.remove(trans_iso)
            alternative_splice_forms = [] # two lists. One for transcripts that actually look like splicing isoforms; one for mutations
            alternative_alleles = [trans_iso] # including the current trans isoform as the first allele
            for transcript in iso_list: # loop through each isoform and compare to the transgressive isoform
                mutation = False # default setting only
                exon_spliced = False 
                run = 0 # measuring the length of each "run" of differences (indel or substitution)
                longest_run = 0
                for p in range(len(trans_seq)): # looping through each bp position in the alignment  
                    if trans_seq[p] != seqs[transcript][p]:
                        run += 1
                    else:
                        if run > 0:
                            if run >= min_exon_size: # checking if this isoform appears spliced from the transgressive isoform (contains >=10bp indel(s)) 
                                exon_spliced = True
                                if run > longest_run:
                                    longest_run = int(run)
                            else: # these are cases where there is a presumptive mutation (<10bp difference)
                                mutation = True
                            # unindent    
                            run = 0
                # (unindent)
                if run > 0: # this is getting the final run after the scan finishes
                    if run >= min_exon_size:
                        exon_spliced = True
                        if run > longest_run:
                            longest_run = int(run)
                    else:
                        mutation = True
                if (exon_spliced == True) and (mutation == False): # looks like a splicing isoform
                    if transcript in isoform_significant_expression_all_parents: # if the other one is expressed significantly in the parents
                        alternative_splice_forms.append( [transcript, longest_run] )
                elif (exon_spliced == False) and (mutation == True): # looks like a mutation on the same isoform
                    alternative_alleles.append(transcript)



                 
            ##### final filters and output #####
            if len(alternative_splice_forms) > 0: # this means there was at least one potential "parent" isoform without mutations
                other_alleles_also_transgressive = True # default setting only
                for allele in alternative_alleles:
                    if allele not in zero_expression_in_parents:
                        other_alleles_also_transgressive = False # an alternate allele for the transgressive looking isoform was found in the parents.
                if other_alleles_also_transgressive == True: # all the alternative forms (if any) check out
                    outline = trans_iso + "\t"
                    outline += str(trans_dict[trans_iso]) + " RILs\t"
                    outline += "alternate alleles: "
                    if len(alternative_alleles) > 1:
                        outline += ",".join(alternative_alleles[1:]) + "\t"
                    else:
                        outline += "NA\t"
                    outline += "parent isoform(s): "
                    parent_isos = []
                    for pi in alternative_splice_forms:
                        parent_isos.append( "|".join(map(str,pi)) )
                    outline += ",".join(parent_isos)
                    print outline



                    
# unindent all the way
os.system( "rm temp1.fa temp2.fa stderr.txt" )







