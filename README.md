# SunflowerTransgressiveSplicing
Code used for sunflower transgressive splicing project

Below is an example workflow using these scripts.


1.  Blast the assembled transcripts against the reference genome. The specific formatting here is important. 

    ###
    blastn -query Trinity_singleline.fasta -db genome.fasta -outfmt "6 qseqid sseqid pident qlen length qstart qend sstart send evalue gaps" -perc_identity 90 > blastOut.txt
    ###




2.  This step "builds" transcript alignments from the individual blast hits (exons), applies some initial filters.
    Heads up: this step uses tons of memory so you need to watch it carefully.

    ###
    python analyzeHAblast_v11.py blastOut.txt 95 > transcript_alignments.txt
    ###

    There is one important parameter that can be played with:
    Minimum percent identity to include an individual blast hit: this is specified on the command line

    Here’s what the script is doing:
    An initial filtering step uses a percent identify cutoff to exclude some blast hits.
    Every possible combination of blast hits is considered, looking for sets of blast hits that align to the same genomic region in a (mostly) non-overlapping arrangement.
    The “best” alignment is the one with the longest alignment length, i.e. sum of the blast hit lengths.
    If the second-longest alignment is within 10bp of the best alignment, I called it ambiguous and threw out the gene.
    If the best alignment length was less than the specified minimum proportion aligned, it was thrown out.

    What the output looks like:
    Some lines are:  good_alignment
    This means that the longest isoform for a gene passed all filters. 
    The subsequent lines are the blast hits that made up the alignment (excluding the other blast hits from this isoform).
    Some lines look like:  not_aligned TRINITY_DN68785_c0_g1_i1
    Which means that the gene* TRINITY_DN68785_c0_g1 did not align well
    I could have output the gene ID instead of the isoform ID, but I think I chose the isoform ID to stay consistent with the other output
    There is no output for a handful of genes. This happens when: 
    All blast hits are below the minimum % ID
    All blast hits are to scaffolds other than the 17 sunflower chromosomes 


    


3.  Next, we want to do a little more filtering to the transcript alignments. 
    Checking if the isoforms identified by Trinity align to the same genomic region. 
    If so, outputting them on the same line. If not, breaking them into multiple "genes". 
    Solo transcripts (no alternative splicing) are excluded at this stage.

    ###    
    python verifySplicing_v5.py transcript_alignments.txt > good_isoforms.txt
    ###









