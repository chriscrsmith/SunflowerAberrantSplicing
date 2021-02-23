

import sys, os
from operator import itemgetter

min_expression_next_closest_isoform = 1 # important parameter; I'm requiring each* parent to have this much of the next-most-similar isoform to the transgressive isoform
min_exon_size = 25 # requiring a >=10bp exon to differentiate the transgressive and the next most similar isoform. 
max_polymorphisms_per_100bp = 5

trans_path = sys.argv[1]
goodisos_path = sys.argv[2]
trinity_path = sys.argv[3]




##### function to consolidate alleles #####
def consolidate(s1, s2):
    mutation_1 = False # default setting only                                               
    exon_spliced_1 = False
    run_diff = 0 # measuring the length of each "run" of differences (indel or substitution)
    longest_run_diff_1 = 0
    diff_locs = []
    for p in range(len(s1)): # looping through each bp position in the alignment            
        if s1[p] != s2[p]:
            run_diff += 1
        else:
            if run_diff > 0:
                diff_locs.append(p)
                if run_diff > min_exon_size: # checking if this isoform appears spliced from the transgressive isoform (contains >25bp indel(s))               
                    exon_spliced_1 = True
                    if run_diff > longest_run_diff_1:
                        longest_run_diff_1 = int(run_diff)
                else:
                    mutation_1 = True
                # unindent
                run_diff = 0
    # (unindent)
    if run_diff > 0: # this is getting the final run after the scan finishes
        diff_locs.append(p)
        if run_diff > min_exon_size:
            exon_spliced_1 = True
            if run_diff > longest_run_diff_1:
                longest_run_diff_1 = int(run_diff)
        else:
            mutation_1 = True
    if len(diff_locs) > max_polymorphisms_per_100bp: # using a sliding window model, if too many snps per 100bp, not an allele, but gene duplication or something.
        for d in range(len(diff_locs)-max_polymorphisms_per_100bp):
            if diff_locs[d+max_polymorphisms_per_100bp] - diff_locs[d] < 100:
                mutation_1 = False
    return mutation_1, exon_spliced_1, longest_run_diff_1













##### read in list of transgressive isoforms #####
list_of_trans = {}
with open(trans_path) as infile:
    for line in infile:
        newline = line.strip().split()
        myiso = newline[0]
        list_of_trans[myiso] = 0








##### read in list of good isoforms #####   *double checked this section 041620
sys.stderr.write("reading in good isos\n")
good_isos = {}
with open(goodisos_path) as infile:
    for line in infile:
        new_gene, isos = line.strip().split()
        isos = isos.split(",")
        for iso in isos:
            if iso in list_of_trans:
                old_gene = new_gene.split(".")[0]
                if old_gene not in good_isos: # some are repeated
                    good_isos[old_gene] = {}
                good_isos[old_gene][new_gene] = [] # adding in the new gene
                for iso in isos:
                    good_isos[old_gene][new_gene].append(iso)




                    

##### Next, I want to create a fasta file with all the (verified) isoforms and align them #####     *double checked this section 041620 
##### (I've thought about just looking at the Trinity node information instead, but there could be both 10bp indels and 1bp SNPs within nodes) #####
sys.stderr.write("doing alignments\n")
for trans_iso in list_of_trans: 
    old_gene = "_".join(trans_iso.split("_")[0:4])
    for new_gene in good_isos[old_gene]: # again, looping through to find the correct new_gene
        if trans_iso in good_isos[old_gene][new_gene]: # this is the correct one
            os.system( "rm temp1.fa temp2.fa 2> stderr.txt" ) # removing these up front, just in case they were left over from before

            
            ###### EMBOSS needle (pairwise sequence alignment) #####
            if len(good_isos[old_gene][new_gene]) == 2:
                ids = list(good_isos[old_gene][new_gene])
                os.system( "grep -w -A 1 " + ids[0] + " " + trinity_path + " >> temp1.fa" )
                os.system( "grep -w -A 1 " + ids[1] + " " + trinity_path + " >> temp2.fa" )
                os.system( "needle temp1.fa temp2.fa -gapopen 10 -gapextend 0.5 -outfile temp1.needle 2> stderr.txt" )
                seqs = {} # initialize empty sequence dictionary
                with open("temp1.needle") as infile:
                    line_number = 0
                    seq1 = ""
                    seq2 = ""
                    for line in infile:
                        if line[0] != "#" and line != "\n" and "|" not in line and "." not in line:
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
                    os.system( "grep -w -A 1 " + iso + " " + trinity_path + " >> temp1.fa" )
                # unindent; only want to run muscle once per gene
                os.system( "muscle -in temp1.fa -out temp1.muscle 2> stderr.txt" )
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
            iso_list = seqs.keys()
            alternative_alleles = [] # including the current trans isoform as the first allele
            for t1 in range(len(iso_list)-1): # loop through each isoform and compare to the transgressive isoform
                for t2 in range(t1+1, len(iso_list)):
                    mutation, exon_spliced, longest_run_diff = consolidate(seqs[iso_list[t1]], seqs[iso_list[t2]])
                    #print mutation, exon_spliced, longest_run_diff
                    if (exon_spliced == False) and (mutation == True): # looks like a mutation on the same isoform
                        alternative_alleles.append([iso_list[t1], iso_list[t2]])
                    else: # looks like a splicing isoform, gene duplication, or something else besides alternate alleles
                        pass


            ##### find pairs that are overlapping #####
            new_set = [] 
            while len(alternative_alleles) > 0:
                delete_list = [ alternative_alleles[0] ] 
                new_cluster = list(alternative_alleles[0]) 
                for pair in range(1, len(alternative_alleles)):
                    for iso in new_cluster:
                        if iso in alternative_alleles[pair]:
                            new_cluster = new_cluster + alternative_alleles[pair]
                            delete_list.append(alternative_alleles[pair])
                new_cluster = list(set(new_cluster))
                new_set.append(new_cluster)
                for pair in delete_list:
                    if pair in alternative_alleles:
                        alternative_alleles.remove(pair)
                    




            ##### output #####
            if len(new_set) > 0:
                outline = []
                for cluster in new_set:
                    outline.append(",".join(cluster))
                print trans_iso + "\t" + "|".join(outline)



                    
# unindent all the way
os.system( "rm temp1.fa temp2.fa stderr.txt" )







