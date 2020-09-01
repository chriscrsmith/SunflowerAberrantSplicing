

import sys, os
from operator import itemgetter

trans_path = sys.argv[1]
fastas = sys.argv[2:]

# first read in the list of transgressive isoforms
trans = {}
with open(trans_path) as infile:
    for line in infile:
        old_iso = line.strip().split()[0]
        gene = "_".join(old_iso.split("_")[0:4])
        if gene not in trans: 
            trans[gene] = []
        trans[gene].append(old_iso)

        
# read in the orfs
orfs = {}
for gene in trans:
    for fasta in fastas:
        iso_id = fasta.split(".")[0].split("/")[-1]
        orfs[iso_id] = []
        with open(fasta) as infile:
            seq = ""
            for line in infile:
                if line[0] == ">":
                    if seq != "":
                        orfs[iso_id].append( [header, seq] )
                        seq = ""
                    # unindent
                    header = line.strip()
                else:
                    seq += line.strip()
            # unindent
            if seq != "": # including this because there are empty ORF files causing issues
                orfs[iso_id].append( [header, seq] ) # this is adding in the last line/orf
trans_iso = trans[gene] # (it looks like all genes are accounted for at this point)


if len(trans_iso) == 1: # some genes have two trans isos... ignoring those for now
    iso_list = orfs.keys()
    trans_index = iso_list.index(trans_iso[0]) # what's the index of the transrgessive iso
    non_trans_index = abs(trans_index - 1) # one liner for getting the other index (i.e. if trans is 0, non trans is 1; if trans is 1, non trans is 0)

    # starting with the non-transgressive isoform, what's the longest ORF?
    if len(orfs[iso_list[non_trans_index]]) == 0:
        longest1 = None
    else:
        longest1 = orfs[iso_list[non_trans_index]][0] # default
        for orf in orfs[iso_list[non_trans_index]][1:]:
            if len(orf[1]) > len(longest1[1]):
                longest1 = list(orf)

    # now the transgressive isoform
    longest2 = [None] # default setting 
    if longest1 and len(orfs[iso_list[trans_index]]) > 0:
        for orf in orfs[iso_list[trans_index]]:
            if orf[1] == longest1[1]: # exact same ORF
                longest2 = list(orf) + ["same"]
            elif len(longest2) == 3: # this occurs when you found an identical ORF already; skip doing the alignment
                pass
            else:
                
                # now, going to want to compare the ORFs with needle alignment
                os.system( "echo \>parent > parent_orf.fa" )
                os.system( "echo " + longest1[1] + " >> parent_orf.fa" )
                os.system( "echo \>trans > trans_orf.fa" )
                os.system( "echo " + orf[1] + " >> trans_orf.fa" )
                os.system( "needle parent_orf.fa trans_orf.fa -gapopen 10 -gapextend 0.5 -outfile orf.needle 2> stderr.txt" )
                seqs = {} # initialize empty sequence dictionary
                with open("orf.needle") as infile:
                    line_number = 0
                    seq1 = ""
                    seq2 = ""
                    for line in infile:
                        if line[0] != "#" and line != "\n" and "|" not in line and "." not in line and ":" not in line:
                            newline = line.strip().split()
                            if newline != []: # this is necessary (filters out stretches of nonmatching positions)
                                line_number += 1
                                if line_number % 2 == 1: # odd numbered lines                                     
                                    seq1 += line.split()[2]
                                else:                    # even numbered lines                                    
                                    seq2 += line.split()[2]
                seqs["parent"] = str(seq1)
                seqs["trans"] = str(seq2)
                os.system( "rm parent_orf.fa trans_orf.fa orf.needle" ) # clean up
                longest_run = 0 # I'm interested in both the longest run of differences, and the total number of differences, at least for now
                current_run = 0
                total_differences = 0 
                for p in range(len(seqs["parent"])):
                    if seqs["parent"][p] != seqs["trans"][p]:
                        current_run += 1
                        total_differences += 1
                    else:
                        if current_run > longest_run:
                            longest_run = int(current_run)
                        # unindent
                        current_run = 0
                if current_run > longest_run: # check the final run of differences at the end of the sequence
                    longest_run = int(current_run)
                if longest2[0] == None:
                    longest2 = list(orf) + [longest_run, total_differences]
                else:
                    if total_differences < longest2[3]:
                        longest2 = list(orf) + [longest_run, total_differences]



    # output
    if longest1 and longest2[0] != None:
        if longest2[2] == "same":
            print trans_iso[0], "relative to", iso_list[non_trans_index], "same longest open reading frame"
        elif len(longest2[1]) < len(longest1[1]): 
            print trans_iso[0], "relative to", iso_list[non_trans_index], "parent ORF shortened by", str(len(longest1[1]) - len(longest2[1])), "aa's, ", str(longest2[3]), "aa total differences"
        elif len(longest2[1]) > len(longest1[1]):
            print trans_iso[0], "relative to", iso_list[non_trans_index], "parent ORF elongated by", str(len(longest2[1]) - len(longest1[1])), "aa's, ", str(longest2[3]), "aa total differences"
        else:
            print trans_iso[0], "relative to", iso_list[non_trans_index], "same length as parent ORF but",  str(longest2[2]), "aa run of mismatches", str(longest2[3]), "aa total differences"
    elif longest1:
        print trans_iso[0], "relative to", iso_list[non_trans_index], "transgressive isoform doesn't have ORF"
    else:
        print trans_iso[0], "relative to", iso_list[non_trans_index], "no ORF in parent isoform"
else:
    sys.stderr.write("\n\ncode broke\n\n")
    1/0









