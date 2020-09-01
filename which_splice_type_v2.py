


import sys, os
from operator import itemgetter

min_exon_length = 25

trans_isos_path = sys.argv[1]
consolidation_path = sys.argv[2]
good_iso_path = sys.argv[3]
parent_list = sys.argv[4]
blast_path = sys.argv[5]



### read in the short list of trans isos
trans_isos = {}
trans_genes = {}
with open(trans_isos_path) as infile:
    for line in infile:
        newline = line.strip().split()
        the_iso, parent_isos = newline[0], newline[-1]
        trans_isos[the_iso] = []
        for iso in parent_isos.split(","):
            trans_isos[the_iso].append(iso.split("|")[0])
        trans_gene = "_".join(the_iso.split("_")[0:4])
        trans_genes[trans_gene] = the_iso


        

### read in consolidated isoforms
consolidation = {}
with open(consolidation_path) as infile:
    for line in infile:
        newline = line.strip().split()
        trans, others = newline[:]
        old_gene = "_".join(trans.split("_")[0:4])
        consolidation[old_gene] = []
        for other in others.split("|"):
            consolidation[old_gene].append(other.split(","))




### read in list of good_isos
good_isos = {}
with open(good_iso_path) as infile:
    for line in infile:
        new_gene, isos = line.strip().split()
        old_gene = new_gene.split(".")[0]
        if old_gene in trans_genes:
            if old_gene not in good_isos: 
                good_isos[old_gene] = {}
            good_isos[old_gene][new_gene] = [] 
            for iso in isos.split(","):
                good_isos[old_gene][new_gene].append(iso)

            


### read in parent data
parents = {}
parent_isos = {} # list of all parent isoforms with >0 expression
with open(parent_list) as infile:
    for line in infile:
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
                    for new_gene in good_isos[old_gene]:
                        if iso in good_isos[old_gene][new_gene]:
                            if new_gene not in parents[parent_id]:
                                parents[parent_id][new_gene] = {}
                            if old_gene in consolidation:
                                added = False
                                for cluster in consolidation[old_gene]:
                                    if iso in cluster:
                                        if cluster[0] not in parents[parent_id][new_gene]:
                                            parents[parent_id][new_gene][cluster[0]] = 0 # going forward, just using the first isoform ID, cluster[0]
                                            parent_isos[cluster[0]] = 0
                                        # unindent
                                        parents[parent_id][new_gene][cluster[0]] += tpm 
                                        added = True
                                if added == False: # this iso wasn't part of a cluster
                                    parents[parent_id][new_gene][iso] = tpm
                                    parent_isos[iso] = 0
                            # unindent
                            else:
                                parents[parent_id][new_gene][iso] = tpm
                                parent_isos[iso] = 0




            
### read in the blast output
blast_data = {}
with open(blast_path) as infile:
    for line in infile:
        newline = line.strip().split() 
        iso = newline[1]
        if iso in trans_isos or iso in parent_isos:
            if iso not in blast_data:
                blast_data[iso] = []
            blast_data[iso].append(newline)



            
            
### compare alignments
sys.stderr.write("")
for tran_iso in trans_isos:
    old_g = "_".join(tran_iso.split("_")[0:4])
    for new_g in good_isos[old_g]: # looping through to find the right new gene       
        if tran_iso in good_isos[old_g][new_g]:
            for par_iso in good_isos[old_g][new_g]:
                if par_iso in parent_isos:
                    found_in_parents = False # quickly check that this iso was expressed at some level in the parents
                    for pare in parents:
                        if parents[pare][new_g][par_iso] > 0:
                            found_in_parents = True
                    if found_in_parents == True:
                        splice_types = []

                        # go through each transgressive exon
                        for hit_num_t in range(len(blast_data[tran_iso])):
                            exclusive = True # default
                            for hit_num_p in range(len(blast_data[par_iso])):
                                ref_start_trans = int(blast_data[tran_iso][hit_num_t][8])
                                ref_start_par = int(blast_data[par_iso][hit_num_p][8])
                                ref_end_trans = int(blast_data[tran_iso][hit_num_t][9])
                                ref_end_par = int(blast_data[par_iso][hit_num_p][9])

                                if abs(ref_start_trans - ref_start_par) < min_exon_length and abs(ref_end_trans - ref_end_par) < min_exon_length: # same reference positions; same exon
                                    exclusive = False

                                elif abs(ref_start_trans - ref_start_par) >= min_exon_length and abs(ref_end_trans - ref_end_par) < min_exon_length: # starts different, ends the same
                                    exclusive =	False
                                    if abs(ref_start_trans - ref_end_trans) > abs(ref_start_par - ref_end_par): # if trans exon longer, see if there is another parent exon with the same start position
                                        alt_splice_site = True # default 
                                        for hit_num_p2 in range(len(blast_data[par_iso])):
                                            if abs(ref_end_trans - int(blast_data[par_iso][hit_num_p2][9])) < min_exon_length:
                                                splice_types.append("intron_retention_trans_retained")
                                                alt_splice_site = False
                                        if alt_splice_site == True:
                                            splice_types.append("alt_splice_site_trans_longer")

                                elif abs(ref_start_trans - ref_start_par) < min_exon_length and abs(ref_end_trans - ref_end_par) >= min_exon_length: # ends different, starts the same
                                    exclusive =	False
                                    if abs(ref_start_trans - ref_end_trans) > abs(ref_start_par - ref_end_par): # if trans exon longer, see if there is another parent exon with the same end position
                                        alt_splice_site = True # default
                                        for hit_num_p2 in range(len(blast_data[par_iso])):
                                            if abs(ref_end_trans - int(blast_data[par_iso][hit_num_p2][9])) < min_exon_length:
                                                splice_types.append("intron_retention_trans_retained")
                                                alt_splice_site = False
                                        if alt_splice_site == True:
                                            splice_types.append("alt_splice_site_trans_longer")

                            if exclusive == True:
                                splice_types.append("exclusive_trans_exon")

                        # go through each parent exon
                        for hit_num_p in range(len(blast_data[par_iso])):
                            exclusive = True # default                                                                                                                                            
                            for hit_num_t in range(len(blast_data[tran_iso])):
                                ref_start_trans = int(blast_data[tran_iso][hit_num_t][8])
                                ref_start_par = int(blast_data[par_iso][hit_num_p][8])
                                ref_end_trans = int(blast_data[tran_iso][hit_num_t][9])
                                ref_end_par = int(blast_data[par_iso][hit_num_p][9])
                                if abs(ref_start_trans - ref_start_par) < min_exon_length and abs(ref_end_trans - ref_end_par) < min_exon_length: # same reference positions; same exon           
                                    exclusive = False
                                elif abs(ref_start_trans - ref_start_par) >= min_exon_length and abs(ref_end_trans - ref_end_par) < min_exon_length: # starts different, ends the same            
                                    exclusive = False
                                    if abs(ref_start_trans - ref_end_trans) < abs(ref_start_par - ref_end_par): # if trans exon shorter, see if there is another trans exon with the same start position
                                        alt_splice_site = True # default                                                                                                                          
                                        for hit_num_t2 in range(len(blast_data[tran_iso])):
                                            if abs(ref_end_par - int(blast_data[tran_iso][hit_num_t2][9])) < min_exon_length:
                                                splice_types.append("intron_retention_parent_retained")
                                                alt_splice_site = False
                                        if alt_splice_site == True:
                                            splice_types.append("alt_splice_site_parent_longer")
                                elif abs(ref_start_trans - ref_start_par) < min_exon_length and abs(ref_end_trans - ref_end_par) >= min_exon_length: # ends different, starts the same            
                                    exclusive = False
                                    if abs(ref_start_trans - ref_end_trans) < abs(ref_start_par - ref_end_par): # if trans exon shorter, see if there is another trans exon iso with the same end position
                                        alt_splice_site = True # default                                                                                                                          
                                        for hit_num_t2 in range(len(blast_data[tran_iso])):
                                            if abs(ref_end_par - int(blast_data[tran_iso][hit_num_t2][9])) < min_exon_length:
                                                splice_types.append("intron_retention_parent_retained")
                                                alt_splice_site = False
                                        if alt_splice_site == True:
                                            splice_types.append("alt_splice_site_parent_longer")
                            if exclusive == True:
                                splice_types.append("exclusive_parent_exon")




                        if len(splice_types) == 0:
                            splice_types.append("other")
                        print tran_iso + "\t" + par_iso + "\t" + ",".join(set(splice_types))


