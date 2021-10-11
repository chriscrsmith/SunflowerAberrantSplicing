

# v3: skipping the new isoform deliniations


import sys

filtered = {}
with open(sys.argv[1]) as infile:
    for line in infile:
        newline = line.strip().split()
        iso = newline[0].split(",")[0]
        if iso not in filtered:
            filtered[iso] = newline

single = {}
with open(sys.argv[2]) as infile:
     for line in infile:
         newline = line.strip().split() 
         iso = newline[1]
         single[iso] = 0

for iso in filtered:
    if iso in single:
        print "\t".join(filtered[iso])
