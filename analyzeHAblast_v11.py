
# last updated: 040920

import sys, numpy
sys.stderr.write("\nWARNING: This script uses tons of memory (potentially >100Gb) and doesn't have a built in fail safe. Need to watch it closely \n")
from operator import itemgetter
from copy import deepcopy

minProportionTranscriptAligned = 0.90 # important parameter: how much of the transcript must align to output it?
maxIntronSize = 100000 # i.e. maximum distance between two blast hits
maxQueryGap = 10 # used when putting together OBblast hits: maximum gap OR OVERLAP between two query start/ends of two blast hits
minExonLength = 10 # used to filter out isoforms if they look the same as other isoforms. If they are almost the same length, then not outputting it.
minPercentID = float(sys.argv[2]) # min %ID per blast hit




                           
#########################################################
# read in and organize all the blast data by chromosome
#########################################################
isos = {} # dictionary of isoforms with corresponding blast hits
sys.stderr.write("reading in blast data\n")
with open(sys.argv[1], 'r') as infile:
    for line in infile:
        newline = line.strip().split()
        iso = newline[0]
        chrom = newline[1]
        percentID = float(newline[2])
        if "00c" in chrom: # these appear to be contigs not anchored to any chromosome
            pass
        elif percentID < minPercentID:
            pass
        else:
            if iso not in isos:
                isos[iso] = {}
            if chrom not in isos[iso]:
                isos[iso][chrom] = []
            isos[iso][chrom].append(newline)







############################################################################
# now go through the genes individually and find good alignments
############################################################################
sys.stderr.write("going through blast hits, finding alignments\n")
for iso in isos:

    ##### go through the hits, and put together exons #####    
    transcripts = {} # building multiple theoretical transcripts aligning to different places; those are stored here.
    hitDict = deepcopy(isos[iso]) # creating new variable with blast hits for the current gene/iso
    for chromosome in hitDict:
        transcripts[chromosome] = []
        for currentHit in hitDict[chromosome]:
            hits = list( hitDict[chromosome] ) # use list() or else python links the variables.
            transcript = list( [currentHit] )  # this is the transcript we will start building.
                                               # The basic idea here, is we want all exons to align, and align close enough together.
                                               # The below code tries "building" different transcripts recursively.
            if ( int(transcript[-1][8]) - int(transcript[0][7])) > 0:
                strand = '+'
            else:
                strand = '-'
            removeList = ['space holder'] # need to define this before the following while loop
            while (len(removeList) > 0): # things get added to remove list as they are added to the growing transcript
                removeList = []
                for hit in hits: # for each hit on the current chromosome (for this isoform)
                    exonRefStart = int(hit[7]) # these are the start/end reference genome positions from the current blast hit we're taking a look at
                    exonRefEnd = int(hit[8])
                    if ( int(exonRefEnd) - int(exonRefStart) ) > 0: # hit is positive (+) strand
                        hitStrand = '+'
                    else:
                        hitStrand = '-'
                    if hitStrand != strand: # if the hits occur in different orientations, then I'm assuming they don't go together
                        pass
                    else: # going to loop through all the blast hits in order, comparing to the ends (only) of the growing transcript
                        exonQueryStart = int(hit[5]) # these are the start/end trinity transcript positions on the current blast hit we're taking a look at
                        exonQueryEnd = int(hit[6])   # these are always in order: start smaller position, end larger position
                        transcriptRefStart = int(transcript[0][7]) # reference genome start position of the FIRST blast hit in the growing transcript
                        transcriptRefEnd = int(transcript[-1][8]) # reference genome end position of the LAST blast hit in the growing transript
                        transcriptQueryStart = int(transcript[0][5]) # trinity transcript start position on the first blast hit in the growing transcript
                        transcriptQueryEnd = int(transcript[-1][6]) # trinity transcript end position on the last blast hit in the growing transcript
                        if abs(exonRefStart - transcriptRefStart) < maxIntronSize: # first, make sure the new hit is in the vicinity (100,000bp?) of the growing transcript
                            if abs(exonQueryStart - transcriptQueryEnd) <= maxQueryGap: # on the trinity transcript, checking the bp gap between new hit and growing transcrpt.
                                                                                        # we want this to be quite small, or else the whole transcript did not align.
                                                                                        # if the gap is small, then this looks like a good blast hit to add onto the transcript.
                                if strand == '+': # check strand: the ref start/end CAN be different orientation than transcript start/end
                                    if exonRefStart <= transcriptRefEnd:
                                        pass # here, the new hit overlaps with the end of the growing transcript.
                                             # if they overlapped substantially, that obviously wouldn't be good.
                                             # But what if they overlapped by 1 bp, or a few?   I'm saying skip it.
                                             # Because, these are reference positions, and blast wouldn't have broken them up into different "hits" if they were connected.
                                    else: # ref positions do NOT overlap, because separated by intron
                                        transcript.append(hit) # add to growing transcript
                                        removeList.append(hit) # remove from the bag-full-of-hits to add on
                                else: # "-" orientation (for both the growing transcript, and the hits being added on)
                                    if exonRefStart < transcriptRefEnd: 
                                        transcript.append(hit) 
                                        removeList.append(hit)
                            elif abs(exonQueryEnd - transcriptQueryStart) <= maxQueryGap: # if exon end is near transcript start
                                if strand == '+':
                                    if exonRefEnd < transcriptRefStart:
                                        transcript.insert(0,hit) # insert to first index (instead of appending to end)
                                        removeList.append(hit)
                                else:
                                    if exonRefEnd > transcriptRefStart:
                                        transcript.insert(0,hit)
                                        removeList.append(hit)           # P.S. I drew each out each of these scenarios, and I believe this is correct.
                            
                for item in removeList:
                    hits.remove(item)
            transcripts[chromosome].append(transcript) # adding transcript, even if single hit


            
    ################################################################################################################################       
    # now go back through all "assembled" transcripts, and choose the best one, or if ambiguous do not output anything for this gene
    ################################################################################################################################
    output = True # default setting: output (one) transcript for the isoform
    best = None 
    topScore = None
    for chrom in transcripts: # remember, there could be good looking alignments on multiple chromosomes (or on the same chromosome)
        for trans in transcripts[chrom]:
            totalAlignLength = 0
            for hit in trans:
                totalAlignLength += int(hit[4]) # total alignment length
            if best == None:
                best = list(trans)
                topScore = int(totalAlignLength)
            elif totalAlignLength > topScore: # new transcript has longest alignment length 
                if (totalAlignLength - topScore) < minExonLength: # if the alignment lengths differ by fewer than 10bp, then one isn't actually better, the same number of exons is aligning. 
                    output = False # so far, the best alignment is the same length as another alignment; therefore the true position is ambiguous. Don't output.
                elif output == False:
                    output = True  # up to this point the best length aligned to two different chromosomes, but now we found a better one. Output this one!
                best = list(trans)
                topScore = totalAlignLength
            elif totalAlignLength < topScore: # not as good as previous best hits; don't update best score                
                if (topScore - totalAlignLength) < minExonLength: # again checking if the alignment lengths are virtually the same; in that case, don't output.
                    output = False
            else: # two transcripts have the exact same aligment length
                if trans == best: # this occurs sometimes becuase we built the same transcript multiple times; we definitely want to keep one of these.
                    pass 
                else: # true position is ambiguous 
                    output = False
                    
    trinity_transcript_length = float(best[0][3])
    if topScore < (trinity_transcript_length * minProportionTranscriptAligned): # make sure enough of the exons aligned
        output = False
    if best:
        if best[0][1] == "Ha0_73Ns": # unaligned/ambiguous contig
            output = False
        






            

#################################################################
# output
################################################################# 
    if output == True:
        print "good_alignment"
        for hit in best:
            print "\t".join(["hit"]+hit)
    else:
        print "not_aligned " + iso 







        
