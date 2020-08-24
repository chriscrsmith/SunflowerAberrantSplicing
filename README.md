# SunflowerTransgressiveSplicing
Code used for sunflower transgressive splicing project

Below is an example workflow using these scripts.


1. Blast the assembled transcripts against the reference genome. The specific formatting here is important. 

     blastn -query Trinity_singleline.fasta -db genome.fasta -outfmt "6 qseqid sseqid pident qlen length qstart qend sstart send evalue gaps" -perc_identity 90 > blastOut.txt



2. This step "builds" transcript alignments from the individual blast hits (exons), applies some initial filters.
   Heads up: this step uses tons of memory so you need to watch it carefully.

     python analyzeHAblast_v11.py blastOut.txt 95 > transcript_alignments.txt

