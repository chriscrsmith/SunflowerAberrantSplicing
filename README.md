# SunflowerAberrantSplicing
Code used for sunflower aberrant splicing project





Below is an example workflow.





1.  Blast the assembled transcripts against the reference genome. The specific formatting here is important. 

    ###
    blastn -query Trinity_singleline.fasta -db genome.fasta -outfmt "6 qseqid sseqid pident qlen length qstart\
    	   qend sstart send evalue gaps" -perc_identity 90 > blastOut.txt
    ###





2.  This step "builds" transcript alignments from the individual blast hits (exons), applies some initial filters.
    Heads up: this step uses tons of memory so you need to watch it carefully.

    ###
    python analyzeHAblast_v11.py blastOut.txt 95 > transcript_alignments.txt
    ###

    There is one important parameter that can be played with: Minimum percent identity to include an individual 
    blast hit. This is specified on the command line

    Here’s what the script is doing: An initial filtering step uses a percent identify cutoff to exclude some blast 
    hits. Every possible combination of blast hits is considered, looking for sets of blast hits that align to the 
    same genomic region in a (mostly) non-overlapping arrangement. The “best” alignment is the one with the longest 
    alignment length, i.e. sum of the blast hit lengths. If the second-longest alignment is within 10bp of the best 
    alignment, I called it ambiguous and threw out the gene. If the best alignment length was less than the specified 
    minimum proportion aligned, it was thrown out.

    What the output looks like: Some lines are: good_alignment. This means that the longest isoform for a gene passed 
    all filters. The subsequent lines are the blast hits that made up the alignment (excluding the other blast hits 
    from this isoform). Some lines look like: not_aligned TRINITY_DN68785_c0_g1_i1 Which means that the gene* 
    TRINITY_DN68785_c0_g1 did not align well I could have output the gene ID instead of the isoform ID, but I think 
    I chose the isoform ID to stay consistent with the other output. There is no output for a handful of genes. 
    This happens when: All blast hits are below the minimum % ID, or all blast hits are to scaffolds other than the 
    17 sunflower chromosomes 


    


3.  Next, we want to do a little more filtering to the transcript alignments. Checking if the isoforms identified by 
    Trinity align to the same genomic region. If so, outputting them on the same line. If not, breaking them into 
    multiple "genes". Solo transcripts (no alternative splicing) are excluded at this stage.

    ###    
    python verifySplicing_v5.py transcript_alignments.txt > good_isoforms.txt
    ###





4.  Once we have identified reasonably confident examples of alternative splicing, use the below script to identify 
    transgressive isoforms. The files parent_list.txt an RIL_list.txt are lists of filepaths to gene expression 
    tables output by Trinity. The output will include isoforms that are found in at least one hybrid with >1 TPM, and
    have 0 reads in all parents.

    ###
    python transSplice_v8.py good_isoforms.txt parent_list.txt RIL_list.txt > transgressive_splicing.txt
    ### 





5.  Before proceding, we need to consolidate similar parent transcripts that represent allelic variation. That way, 
    we can sum each samples' expression across all alleles for the same splice isoform.
 
    ###
    python consolidate_isoforms_v2.py transgressive_splicing.txt good_isoforms.txt Trinity_singleline.fasta \
    	   > transgressive_splicing_consolidatedIsos.txt
    ###





6.  Here, we try to understand which splice type is occuring, i.e. intron retention, exon skipping, etc.

    ### 
    python which_splice_type_v2.py transgressive_splicing.txt\
    	   TransgressiveSegregation/transgressive_splicing_consolidatedIsos.txt\
	   good_isoforms.txt parent_list.txt transcript_alignments.txt > transgressive_splicing_spliceType.txt
    ###





7.  For each transgressive isoform, testing for a difference in splicing composition between parent taxa.

    ###
    python parents_diff.py transgressive_splicing.txt good_isoforms.txt transgressive_splicing_consolidatedIsos.txt\
    	   parent_list.txt 
    ###





8.  Code to compare longest open reading frames between a transgressive isoform and parent isoform. The first
    input file is a file containing the transgressive isoform ID. The second and third files are output from
    ORFfinder.

    ###
    python look_at_orfs_v3.py list_of_isoformIDs trans_isoform.orf parent_isoform.orf
    ###





9.  Grabbing the proportion of overall gene expression for each sample, for each transgressive isoform. Will 
    use this for QTL mapping (below).

    ###
    python transPhenoCols_v7.py good_isoforms.txt transgressive_splicing.txt RIL_list.txt > phenotype_table.txt
    ###





10. Taking a step back, processing SNPs resulting from bwa alignments + samtools mpileup. First, identifying 
    SNPs that are fixed between parent genotypes.

    ###
    python parent_fixed_snps_v2.py parents.vcf > fixed_snps.txt
    ###





11. Finding genomic bp position for each snp (the snps were obtained from RNAseq data aligned to HA412 transcriptome). 
    Blast Ha412 trancriptome to Ha412 genome.

    ###
    blastn -query HA412_trinity_noAltSplice_400bpmin.fa -db Han412-HO.fa -outfmt "6 qseqid sseqid pident qlen\
    	   length qstart qend sstart send evalue gaps" -perc_identity 90 > Ha412_transcriptome_x_Ha412_genome_blastOut.txt    
    ###





12. Parse the blast output to identify the genomic position of each transcript.

    ###
    python analyzeHAblast_v9.py HA412_trinity_noAltSplice_400bpmin.fa Ha412_transcriptome_x_Ha412_genome_blastOut.txt\
    97.5 > bpPositions_Ha412_x_Ha412.txt
    ###
  
 



13. Use the transcript alignments to find genomic bp position for each SNP.
    
    ###
    python grab_pos.py bpPositions_Ha412_x_Ha412.txt fixed_snps.txt > fixed_snps_bpPositions.txt
    ###





14. Interpolate cM positions for each fixed snp.

    ###
    python interpolateMapPos_v5.py fixed_snps_bpPositions.txt Ha412HOv2.0-20181130.Nov22k22.geneticmap.txt\
    	   > fixed_snps_cMpositions.txt
    ###





15. Initial SNP filtering in the RILs based on minimum read depth: min depth 10 to include a genotype.

    ###
    python filterSNPs_v3.py fixed_snps_cMpositions.txt RIL_list.txt RIL_VCFs/ 10 > genotype_table_min10reads.txt
    ###





16. Impute missing data.

    ###
    python impute.py genotype_table_min10reads.txt > genotype_table_imputed.txt
    ###





17. Apply final SNP filter: 35% RILs with non-missing data. 

    ###
    python filterSNPs_v4.py genotype_table_imputed.txt > genotype_table.txt 
    ###





18. Epistasis scan. This script has can be run a few different ways.

    ###
    python interaction_scan_v8.py genotype_table.txt phenotype_table.txt <phenotype #> <SNP #> <permute?>
    ### 

    The third command line parameter is the phenotype column # in the phenotype table, minus 1. I.e.
    the number of columns after the RIL ID column.
    The fourth command line parameter is the SNP # in the genotype table, minus 1. This SNP is tested
    against every other SNP (except SNPs on the same chromosome). Alternatively, you may use "all"
    to test every pair of SNPs (except SNPs on the same chromosome).
    The fifth command line parameter specifies whether or not to run permutations. Options include
    "no", "yes", or "stratified". If "yes" is specified, the phenotype column is randomized. If 
    "stratified" is specified, the phenotypes are randomized within each genotype group, as in
    Doerge and Churchill (Doerge RW, Churchill GA. Permutation tests for multiple loci affecting a 
    quantitative character. Genetics. 1996 Jan 1;142(1):285-94.).

    Ex 1. To do the empirical whole genome scan for the first phenotype:

    ###
    python interaction_scan_v8.py genotype_table.txt phenotype_table.txt 1 all no
    ###

    Ex 2. To do stratified permutation analysis for 2nd phenotype, and locus number 10143:

    ###
    python interaction_scan_v8.py genotype_table.txt phenotype_table.txt 2 10143 stratified
    ###





19. Identifying parent isoforms that were missing in hybrids

    ###
    python other_splicing_errors_v4.py parent_list.txt RIL_list.txt
    ###





20. Filter the missing isos for single copy regions

    ###
    python filter_single_v3.py filtered_missing.txt transcript_alignments_singleCopy_50perc.txt > missing_parent_isoforms_singleCopy.txt
    ###







