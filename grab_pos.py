


import sys

bp_map = {}
with open(sys.argv[1]) as infile:
    infile.readline()
    for line in infile:
        newline = line.strip().split()
        if newline[0] == "not_aligned":
            pass
        elif newline[0] == "good_alignment":
            pass
        elif newline[0] == "hit":
            transcriptID = newline[1]
            newline[6] = int(newline[6])
            newline[7] = int(newline[7])
            newline[8] = int(newline[8])
            newline[9] = int(newline[9])
            if transcriptID not in bp_map:
                bp_map[transcriptID] = []
            bp_map[transcriptID].append(newline)
        else:
            1/0

with open(sys.argv[2]) as infile:
    infile.readline()
    for line in infile:
        newline = line.strip().split()
        transcriptID, pos, a1, a2 = newline[:]
        pos = int(pos)
        output = False
        exon_aligned = False
        if transcriptID in bp_map:
            for hit in bp_map[transcriptID]: # find which blast hit (exon) the snp is on
                tranStart, tranEnd = hit[6:8]
                if pos >= tranStart and pos <= tranEnd: # this is the right blast hit
                    refStart, refEnd = hit[8:10]
                    exon_aligned = True
                    if refStart < refEnd: # forward-orientation
                        output = [ hit[2], refStart+(pos-tranStart) ]
                    else: # reverse compliment
                        output = [ hit[2], refStart-(pos-tranStart) ]
                        
            #            
            if exon_aligned == False: # here, the transcript with the snp aligned "well", but the specific (small) exon the snp was on did not align
                for exon in bp_map[transcriptID]:
                    if pos < exon[6]: # going through each exon in order, checking if the SNP position comes before* the start of the exon ailgnment
                        tranStart, tranEnd = exon[6:8]   
                        refStart, refEnd = exon[8:10]
                        if refStart < refEnd:
                            output = [ exon[2], refStart-(tranStart-pos) ]
                        else:
                            output = [ exon[2], refStart+(tranStart-pos) ]
                        exon_aligned = True
                        break # break from the for-loop

                #    
                if exon_aligned == False: # here, the SNP position comes after* the last exon
                    last_exon = bp_map[transcriptID][-1]
                    tranStart, tranEnd = last_exon[6:8]
                    refStart, refEnd = last_exon[8:10]
                    if refStart < refEnd:
                        output = [ exon[2], refEnd+(pos-tranEnd) ]
                    else:
                        output = [ exon[2], refEnd-(pos-tranEnd) ]

        #
        if output:
            print "\t".join( [transcriptID, str(pos), output[0], str(output[1]), a1, a2] )
        else:
            print "\t".join( [transcriptID, str(pos), "NA", "NA", "NA", "NA"] )


    
