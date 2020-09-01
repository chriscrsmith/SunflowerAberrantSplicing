
import sys, decimal, numpy
decimal.getcontext().prec = 22 
snpPath = sys.argv[1] # snp bp positions on genome, chrom tab bp
mapPath = sys.argv[2] # genetic map
print "\t".join( ["transcript", "transcript_bp", "chrom", "chrom_bp", "cM", "HA89_allele", "1238_allele" ]) # header






# let's do this by linkage group, to make list-searching quick
for lg in range(1,18)+["cp_gi_88656873", "mt_gi_571031384", "rDNA_gi_563582565"]: # 1-17, plus the organelles
    currentChrom = str(lg)

    # make dictionary of cM positions; we want to exclude redundant bp positions with the same cM position
    map_bins = {}
    map_positions = {}
    with open(mapPath, 'r') as mapFile:
        mapFile.readline()
        for line in mapFile:
            newline = line.strip().split()
            chrom = str(int(newline[0].split("Chr")[1]))
            if chrom == currentChrom:
                bp, cM = newline[1], newline[3]
                if cM not in map_bins:
                    map_bins[cM] = int(bp)
                    map_positions[int(bp)] = cM

    # then go through your list and find positions for each one
    with open(snpPath, 'r') as snpFile:
        for line in snpFile:
            newline = line.strip().split()
            geneID = newline[0]
            oldPos = newline[1]
            chrom = newline[2]
            if newline[2] == "NA":
                pass
            elif "cp" in chrom or "mt" in chrom or "rDNA" in chrom:
                pass
            else:
                chrom = str(int(chrom.split("Chr")[1]))
                if chrom == currentChrom:
                    bp = int(newline[3])
                    mapList = map(int, map_positions.keys())
                    mapList.sort()
                    for map_pos in range(len(mapList)):
                        currentPos = int(mapList[map_pos])
                        if bp >= currentPos: 
                            beforePosition = map_pos

                    # output
                    afterPosition = beforePosition + 1
                    if afterPosition < len(mapList):
                        beforeBP = int(mapList[beforePosition])
                        afterBP = int(mapList[afterPosition])
                        beforecM = float(map_positions[beforeBP])
                        aftercM = float(map_positions[afterBP])
                        refdiff = afterBP - beforeBP
                        snpdiff = int(bp) - beforeBP
                        ratio = decimal.Decimal(snpdiff) / decimal.Decimal(refdiff)
                        aftercM = float(map_positions[afterBP])
                        cMdiff = aftercM - beforecM
                        finalcM = decimal.Decimal(beforecM) + (decimal.Decimal(cMdiff)*decimal.Decimal(ratio))
                        print "\t".join( map(str, [geneID, oldPos, chrom, bp, finalcM, newline[4], newline[5] ]))
                    else:
                        print "\t".join( map(str, [geneID, oldPos, chrom, bp, map_positions[mapList[beforePosition]], newline[4], newline[5] ]))

snpFile.close()
