

DPcutoff = 10 # this is the cutoff we will use for the DP values
              # the median in all six parents, unfiltered, is = 20.
import sys
vcfpath = sys.argv[1]
vcf = open(vcfpath, 'r')

# header
print ' '.join([ 'contig', 'pos', 'HA89version', '1238version'])

nuclist = ['A', 'C', 'T', 'G']
geno = 0 # index for the individual genotype
dp = 2 # index for dp field

for line in vcf:
    if line[0] != '#':
        newline = line.strip().split()
        A1238info = newline[9].split(':')
        b1238info = newline[10].split(':')
        e1238info = newline[11].split(':')
        HA89Ainfo = newline[12].split(':')
        HA89binfo = newline[13].split(':')
        HA89einfo = newline[14].split(':')

        if ( A1238info[geno] == '0/0' and
             b1238info[geno] == '0/0' and
             e1238info[geno] == '0/0' and
             HA89Ainfo[geno] == '1/1' and
             HA89binfo[geno] == '1/1' and
             HA89einfo[geno] == '1/1' ):    
            if ( int(A1238info[dp]) >= DPcutoff and
                 int(b1238info[dp]) >= DPcutoff and
                 int(e1238info[dp]) >= DPcutoff and
                 int(HA89Ainfo[dp]) >= DPcutoff and
                 int(HA89binfo[dp]) >= DPcutoff and
                 int(HA89einfo[dp]) >= DPcutoff ):
                if (newline[3] in nuclist) and (newline[4] in nuclist): # it turns out this is redundant, because you already filtered the "2" allele. whatev.
                    HA89version = newline[4]
                    version1289 = newline[3]
                    print newline[0], newline[1], HA89version, version1289

        if ( A1238info[geno] == '1/1' and
             b1238info[geno] == '1/1' and
             e1238info[geno] == '1/1' and
             HA89Ainfo[geno] == '0/0' and
             HA89binfo[geno] == '0/0' and
             HA89einfo[geno] == '0/0' ):
            if ( int(A1238info[dp]) >= DPcutoff and
                 int(b1238info[dp]) >= DPcutoff and
                 int(e1238info[dp]) >= DPcutoff and
                 int(HA89Ainfo[dp]) >= DPcutoff and
                 int(HA89binfo[dp]) >= DPcutoff and
                 int(HA89einfo[dp]) >= DPcutoff ):
                if (newline[3] in nuclist) and (newline[4] in nuclist):
                    HA89version = newline[3]
                    version1289 = newline[4]
                    print newline[0], newline[1], HA89version, version1289

vcf.close()


