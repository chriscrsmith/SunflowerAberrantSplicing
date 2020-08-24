

# updated 030620

import sys
from operator import itemgetter

transcriptPath = sys.argv[1] 
maxIntronSize = 100000


# first, simply load in the trancript alignment data
sys.stderr.write("reading in transcript alignment data\n")
transcripts = {} # dictionary full of all the transcripts we built
with open(transcriptPath, 'r') as infile:
    for line in infile:
        newline = line.strip().split()
        if newline[0] == "hit":
            chrom = newline[2] 
            transcriptID = newline[1]
            geneID = "_".join(transcriptID.split("_")[:-1])
            newline[8] = int(newline[8]) # changing the reference positions to int() type in the list
            newline[9] = int(newline[9])
            if newline[8] > newline[9]: # if the end position comes first, go ahead and switch (the query start/end are already good)
                newline[8],newline[9] = newline[9],newline[8]
            if geneID not in transcripts:
                transcripts[geneID] = {}
            if transcriptID not in transcripts[geneID]:
                transcripts[geneID][transcriptID] = []
            transcripts[geneID][transcriptID].append(newline[1:-2]) # cutting out some extraneous fields ("hit", and the e value and bit score)
deleteList = [] # this is for removing genes with only a single transcript (no splicing)
for geneID in transcripts:
    if len(transcripts[geneID]) == 1: 
        deleteList.append(geneID)
for geneID in deleteList:
    del transcripts[geneID]
for geneID in transcripts: # for each transcript, sorting the blast hits by reference start position       
    for transcriptID in transcripts[geneID]:
        orderedList = sorted(transcripts[geneID][transcriptID], key=itemgetter(7))
        transcripts[geneID][transcriptID] = list(orderedList)


    



    
##### check: do transcripts align to the same genomic locus/region, or not? #####
sys.stderr.write("checking that transcripts align to the same genomic locus/region\n")
expandedTranscripts = {} # "expanded", because in some cases we want to break up what trinity calls a "gene" into multiple genes because they align different places
for gene in transcripts:
    numIsos = len(transcripts[gene])
    chroms = {} # keeping track of which chromosome(s) transcripts are aligning to
    for isoform in transcripts[gene]:
        chrom = transcripts[gene][isoform][0][1] # (note: there is one mtDNA gene, at this point)
        if chrom not in chroms:
            chroms[chrom] = []
        chroms[chrom].append(isoform)
    numDiffChroms = len(chroms)
    geneNumber = -1 # intitializing to "-1", so that the first gene ID is "0"
    for chrom in chroms:
        if len(chroms[chrom]) == 1:
            pass # only one transcript aligns to this chromosome; no splicing happening
        else:
            geneNumber += 1 # begin a new geneID for each new chromosome
            positions = [] # go through the isos aligning to the current chrom, and collect their reference positions
            for iso in chroms[chrom]:
                positions.append([ iso, transcripts[gene][iso][0][7] ])
            positions = sorted(positions, key=itemgetter(1))
            previousIsoform = positions[0] # start with the first one, before looping through the rest
            expandedTranscripts[ gene+"."+str(geneNumber) ] = {} # this is the first gene; we may expand it to be more than one gene
            expandedTranscripts[ gene+"."+str(geneNumber) ][ positions[0][0] ] = list(transcripts[gene][positions[0][0]]) # default first transcript for this gene
            for iso in range(1, len(positions)): # loop through the remaining isoforms
                currentIsoform = list(positions[iso])
                if abs(currentIsoform[1] - previousIsoform[1]) < maxIntronSize: # aligned to same region if less than 100000bp away
                    expandedTranscripts[ gene+"."+str(geneNumber) ][ positions[iso][0] ] = list(transcripts[gene][positions[iso][0]])
                else:
                    geneNumber += 1 # isoforms align to different locations on the same chromosome; assign new geneID
                    expandedTranscripts[ gene+"."+str(geneNumber) ] = {}
                    expandedTranscripts[ gene+"."+str(geneNumber) ][ positions[iso][0] ] = list(transcripts[gene][positions[iso][0]])
                previousIsoform = list(currentIsoform)
deleteList = [] # remove new gene-clusters with only a single isoform
for gene in expandedTranscripts:
    if len(expandedTranscripts[gene]) == 1:
        deleteList.append(gene)
for gene in deleteList:
    del expandedTranscripts[gene]

    
    

##### output #####
for gene in expandedTranscripts:
    outline = [gene]
    isos = []
    for iso in expandedTranscripts[gene]:
        isos.append( iso ) 
    outline.append(",".join(isos))
    print "\t".join(outline) 

            









