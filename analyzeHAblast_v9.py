

import sys, numpy
from operator import itemgetter

minProportionTranscriptAligned = 0.66 # important parameter: how much of the transcript must align to output?
minPercentID = float(sys.argv[3]) # min %id to consider individual blast hits
maxIntronSize = 100000 # maximum space between blast hits 
maxQueryGap = 10 # used when putting together blast hits: maximum gap OR OVERLAP between two query* start/ends of two blast hits
minExonLength = 10 # used to filter out isoforms if they look the same as other isoforms. If they are almost the same length, then not outputting it.



#############################################
# read in transcriptome and find longest iso
#############################################
lengthDictionary = {}
with open(sys.argv[1], "r") as infile:
    for line in infile:
        if line[0] == ">":    
            newline = line.strip().split()
            gene = "_".join(newline[0].split(">")[1].split("_")[0:4])
            iso = newline[0].split(">")[1].split("_")[4]
        else:
            length = len(line.strip())
            if gene not in lengthDictionary:
                lengthDictionary[gene] = [iso, length]
            else:
                if length > lengthDictionary[gene][1]:
                    lengthDictionary[gene] = [iso, length]




                    
                           
#########################################################
# read in and organize all the blast data by chromosome
#########################################################
genes = {} # dictionary of isos with corresponding blast hits (confusing variable names, sorry)
with open(sys.argv[2], 'r') as infile:
    for line in infile:
        newline = line.strip().split() 
        iso = newline[0]
        gene = '_'.join(iso.split('_')[:-1])
        chrom = newline[1]
        percentID = float(newline[2])
        if "00c" in chrom: # these appear to be contigs unanchored to any chromosome; skippem!
            pass
        elif percentID < minPercentID:
            pass
        else:
            if iso.split('_')[-1] == lengthDictionary[gene][0]: # if there are blast hits for the longest isoform, save the data. Otherwise, we're skipping the blast hit
                if iso not in genes:
                    genes[iso] = {}
                if chrom not in genes[iso]:
                    genes[iso][chrom] = []
                genes[iso][chrom].append(newline)
                    







############################################################################
# now go through the genes individually and find good alignments
############################################################################
for gene in genes:

    ##### go through the hits, and put together exons #####    
    transcripts = {} # building multiple theoretical transcripts aligning to different places
    hitDict = genes[gene] # creating new variable with blast hits for the current gene/iso
    for chromosome in hitDict:
        transcripts[chromosome] = []
        for currentHit in hitDict[chromosome]:
            hits = list( hitDict[chromosome] ) # use list() or else python links the variables
            transcript = list( [currentHit] )  # this is the transcript we will start building
                                               # The basic idea here, is we want all exons to align, and close enough together
                                               # This loop tries "building" different transcripts, starting the transcript from each different isoform
                                               # I think I chose to do this because sometimes transcripts align to various places on a chromosome
            if ( int(transcript[-1][8]) - int(transcript[0][7])) > 0:             
                strand = '+'
            else:
                strand = '-'
            # recursively add on exons to the growing transcript 
            removeList = ['space_holder'] # need to define this before the following while loop
            while (len(removeList) > 0): # things get added to remove list as they are added to the growing transcript
                removeList = []
                for hit in hits:
                    exonRefStart = int(hit[7])
                    exonRefEnd = int(hit[8])
                    if ( int(exonRefEnd) - int(exonRefStart) ) > 0: # hit is positive (+) strand
                        hitStrand = '+'
                    else:
                        hitStrand = '-'
                    if hitStrand == strand: # (not sure why I never reverse complimented, I guess the hits only occur one way or the other)
                        exonQueryStart = int(hit[5])
                        exonQueryEnd = int(hit[6])
                        transcriptRefStart = int(transcript[0][7])
                        transcriptRefEnd = int(transcript[-1][8])
                        transcriptQueryStart = int(transcript[0][5])
                        transcriptQueryEnd = int(transcript[-1][6])
                        if abs(exonRefStart - transcriptRefEnd) < maxIntronSize: # making sure hits do not have too-large of gaps between them 
                            if abs(exonQueryStart - transcriptQueryEnd) <= maxQueryGap: # if exon start is near transcript end
                                if strand == '+':
                                    if exonRefStart > transcriptRefEnd: # ref positions do NOT overlap, because separated by intron
                                        transcript.append(hit)
                                        removeList.append(hit)
                                else:
                                    if exonRefStart < transcriptRefEnd:
                                        transcript.append(hit)
                                        removeList.append(hit)
                            elif abs(exonQueryEnd - transcriptQueryStart) <= maxQueryGap: # if exon end is near transcript start
                                if abs(exonQueryEnd - transcriptQueryStart) < maxIntronSize:
                                    if strand == '+':
                                        if exonRefEnd < transcriptRefStart:
                                            transcript.insert(0,hit)
                                            removeList.append(hit)
                                    else:
                                        if exonRefEnd > transcriptRefStart:
                                            transcript.insert(0,hit)
                                            removeList.append(hit)
                for item in removeList:
                    hits.remove(item)
            transcripts[chromosome].append(transcript) # adding transcript, even if single hit



            
    ################################################################################################################################       
    # now go back through all "assembled" transcripts, and choose the best one, or if ambiguous do not output anything for this gene
    ################################################################################################################################
    output = True # default
    best = None
    topScore = None
    for chrom in transcripts:
        for trans in transcripts[chrom]:
            totalAlignLength = 0
            for hit in trans:
                totalAlignLength += int(hit[4]) # total alignment length
            if best == None:
                best = list(trans)
                topScore = int(totalAlignLength)
            elif totalAlignLength > topScore: # new best transcript has longest alignment length totalled across exons
                if (totalAlignLength - topScore) < minExonLength: # if the alignment lengths differ by fewer than 10bp, then one isn't actually better, the same number of exons is aligning. 
                    output = False
                elif output == False:
                    output = True  # up to this point the best length aligned to two different chromosomes, but now we found a better one 
                best = list(trans)
                topScore = totalAlignLength
            elif totalAlignLength < topScore: # not as good as previous best hits; don't update best score                
                if (topScore - totalAlignLength) < minExonLength: # again checking if the alignment lengths are too close
                    output = False
            else: # two transcripts have the exact same aligment length
                if trans == best: # this only happens for the first transcript
                    pass 
                else: 
                    output = False
                
    if topScore < (lengthDictionary["_".join(gene.split("_")[0:4])][1] * minProportionTranscriptAligned): # make sure enough of the exons aligned
        output = False
    if best != None:
        if best[0][1] == "Ha0_73Ns": # counts as unaligned/ambiguous
            output = False
        






            

#################################################################
# output
################################################################# 
    if output == True:
        print "good_alignment"
        for hit in best:
            print "\t".join(["hit"]+hit)
    else:
        print "not_aligned " + gene







        
