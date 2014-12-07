NOTE: This is not a finished project and the scripts will be combined for better use. Look at scripts/classify.pl for details.

This project was initiated by Prof. Jae Hyung Lee to look deeper into a lesser-known strain of Prevotella Intermedia and compare it with a well-known strain to recognize its differences and relative functions.

Prevotella Intermedia is known to cause various periodontal diseases, including periodontitis and gingivitis. Despite its important role, this microorganism is not yet heavily studied, and only one of its strains, P. Intermedia 17, has its complete DNA sequences assembled. Meanwhile, a lesser known strain, ATCC 25611, is assembled into 24 DNA contigs, with its coding sequences not fully known. We intend to compare this lesser known strain 25611 to a well known strain 17, thereby further understanding the functional differences between the two strains and ultimately building a global alignment of strain 25611 DNA sequences compared to strain 17.

To achieve this, we used BLAST to compare the proteins of Prevotella Intermedia 17 (PI 17) with open reading frames(ORFs) of Prevotella Intermedia 25611 (PI 25611), and MAUVE to compare the two chromosomes of PI 17 and 24 contains of PI 25611. Once we got the initial alignment data, we classified the protein-ORF pairs according to identity percentages so that we would be able to distinguish more identical pairs from less identical ones. Single Nucleotidal Polymorphisms (SNPs) and indels (insertion/deletions) are also calculated according to this data.


This file lists the files and folders in the original workspace and how they work in the original context.

BiO/apps/data/JM/comparative/PI/
- folder to store Prevotella Intermedia related documents, data, and scripts
adjusted_locations.txt
- a list of PI 25611 shotgun contigs and their relative locations on PI 17 chromosomes
alignment3
- alignment file by MAUVE, required for alignment adjustment
alignment_scores.txt
- a list of scores for each possible location shift candidate
amino.global.aligned.txt
- result of Needleman-Wunsch algorithm with amino acid sequences
amino.mauve.aligned.txt
- result of Needleman-Wunsch algorithm with amino acid anchor-fitting sequences
BLASTp_aligned_matches.txt
- a list of ORF-protein pairs that fit in the alignment locations
BLASTp_all_matches.txt
- a combination of aligned and unaligned ORF-protein pairs
BLASTp_best_matches.txt
- a list of first choice ORF-protein pairs that fit in the alignment locations
BLASTp_consecutive.txt
- a list of ORF-protein pairs that form an ordered block
BLASTp_exclusive.txt
- a list of ORF-protein pairs that are not included in best matches
BLASTp_mauve_best.txt
- a list of ORF-protein pairs that fit in the MAUVE anchor locations with low e values
BLASTp_mauve_matches.txt
- a list of ORF-protein pairs that fit in the MAUVE anchor locations
BLASTp_unaligned_matches.txt
- a list of ORF-protein pairs that do not fit in the alignment locations
BLASTp_parsed.txt
- parsed result of BLASTp analysis
BLASTp_reverse.txt
- parsed result of BLASTp analysis, PI17 protein as query and PI25611 ORF as database
GCA_000439065.1_ASM43906v1_genomic.fna
- shotgun sequences of Prevotella intermedia PI 25611 from NCBI
GCF_000261025.1_ASM26102v1_genomic.fna
- nucleotide sequences of Prevotella intermedia PI17 from NCBI
lastz_result.lav
- lastz result of PI17 chromosome 1 and PI 25611 in lav format
lastz_result.maf
- lastz result of PI17 chromosome 1 and PI 25611 in maf format
mauve_matches_refbase.txt
- a short list of ORF-protein pairs that fit in MAUVE anchors, in PI17 numbering order
mauve_matches.txt
- a short list of ORF-protein pairs that fit in MAUVE anchors
mauve_unique.txt						
- a list of ORF-protein pairs that are unique to MAUVE anchor matches
NC_017860.faa							
- amino-level protein sequences of PI17 chromosome 1 from NCBI
NC_017860.ffn							
- nucleotide-level protein sequences of PI17 chromosome 1 from NCBI
NC_017860.fna							
- nucleotide sequences of PI17 chromosome 1 from NCBI
NC_017860.gff							
- genetic locations and information of PI17 chromosome 1 from NCBI
NC_017860.ptt							
- referential data of PI17 chromosome 1 from NCBI
NC_017861.faa							
- amino-level protein sequences of PI17 chromosome 2 from NCBI
NC_017861.ffn							
- nucleotide-level protein sequences of PI17 chromosome 2 from NCBI
NC_017861.fna							
- nucleotide sequences of PI17 chromosme 2 from NCBI
NC_017861.gff							
- genetic locations and information of PI17 chromosome 2 from NCBI
NC_017861.ptt							
- referential data of PI17 chromosome 2 from NCBI
nucleotide.global.aligned.txt					
- result of Needleman-Wunsch algorithm with nucleotide sequences
PI17.blastp.result.rev.txt					
- BLASTp result with PI25611 as reference and PI17 as query
PI17.blastp.result.txt						
- BLASTp result with PI17 as profile and PI 25611 as query
PI17.db.phr							
- part of PI17 BLAST profile
PI17.db.pin							
- part of PI17 BLAST profile
PI17.db.psq							
- part of PI17 BLAST profile
PI17.faa							
- combined amino-level protein sequence of PI17
PI17.ffn							
- combined nucleotide-level protein sequences of PI17
PI17.mauve.faa							
- combined amino-level anchor-matching protein sequences sorted for global alignment
PI17.mauve.fna							
- combined nucleotide-level anchor-matching protein sequences sorted for global alignment
PI17.ordered.faa						
- combined amino-level protein sequences sorted for global alignment
PI17.ordered.fna						
- combined nucleotide-level protein sequences sorted for global alignment
PI17.ptt							
- combined referential data of PI17
PI25611.contigs.tab						
- a summary of alignment results by MAUVE, required for alignment adjustment
PI25611.mauve.aa.fasta						
- combined amino-level anchor-matching ORF sequences sorted for global alignment
PI25611.mauve.na.fasta						
- combined nucleotide-level anchor-matching ORF sequences sorted for global alignment
PI25611.ordered.aa.fasta					
- combined amino-level ORF sequences sorted for global alignment
PI25611.ordered.na.fasta					
- combined nucleotide-level ORF sequences sorted for global alignment
PI.amino.indels.txt						
- list of indels and mismatches between PI17 and PI25611 on amino acid level
PI.amino.scores.txt						
- list of scores related to PI17 and PI25611 alignment on amino acid level
PI_ATCC_25611_DSM_20706.aa.fasta				
- amino-level open reading frames of PI 25611 from HOMD
PI_ATCC_25611_DSM_20706.aa.fasta.phr				
- part of PI25611 BLAST profile
PI_ATCC_25611_DSM_20706.aa.fasta.pin				
- part of PI25611 BLAST profile
PI_ATCC_25611_DSM_20706.aa.fasta.psq				
- part of PI25611 BLAST profile
PI.mauve.amino.indels.txt					
- list of indels and mismatches between anchor-matching pairs on amino acid level
PI.mauve.amino.scores.txt					
- list of scores related to anchor-matching pairs on amino acid level
PI.mauve.nucleotide.indels.txt					
- list of indels and mismatches between anchor-matching pairs on nucleotide level
PI.mauve.nucleotide.scores.txt					
- list of scores related to anchor-matching pairs on nucleotide level
PI.nucleotide.indels.txt					
- list of indels and mismatches between PI17 and PI25611 on nucleotide level
PI.nucleotide.scores.txt					
- list of scores related to PI17 and PI25611 alignment on nucleotide level
PI.primers.txt							
- text file containing primer sequences for requested PI proteins
Prevotella_Intermedia_ATCC25611_DSM20706.na.fasta		
- nucleotide-level open reading frames of PI 25611 from HOMD

/../PI/GCA/							
- shotgun sequences of PI 25611 parsed into separate files. format = shotgunseg.num.fna

/../PI/GCF/							
- chromosome sequences of PI 17 parsed into separate files. format = refseg.num.fna

/../PI/hmmer/							
- alignment tool hmmer installed folder

/../PI/lastz/							
- alignment tool lastz installed folder

/../PI/needleman_wunsch						
- global alignment tool that implements Needleman-Wunsch algorithm

/../PI/newAlign/						
- alignment results from MAUVE

/../scripts/
alists_indel.pl							
- script to generate an amino acid level indel/mismatch list for aligned pairs
alists_mauve.pl							
- script to generate an amino acid level indel/mismatch list for anchor-matching pairs
anchors_mauve.pl						
- script to use anchors to readjust the results from MAUVE
ascores_indel.pl						
- script to generate an amino acid level list of scores and locations for aligned pairs
ascores_mauve.pl						
- script to generate an amino acid level list of scores and locations for anchor-matching pairs
blast_alignment.pl						
- script to filter the BLASTp results using alignment locations
blast_analysis.pl						
- script to parse results from BLASTp
blast_anchors.pl						
- script to filter the BLASTp results using anchor locations
blastp.pl							
- script to launch BLASTp
chart_calc.pl							
- script to make raw data for drawing scores by location chart, no longer used
classify.pl							
- COMPLETE SCRIPT OF THE ALIGNMENT PIPELINE PROCESS (final script)
comp_aligns.pl							
- script to compare pair numbers between all anchor matches and all aligned matches
comp_baligns.pl							
- script to compare pair numbers between best anchor matches and all aligned matches
comp_best.pl							
- script to compare pair numbers between best anchor matches and best aligned matches
count_reliable.pl						
- script to calculate mean and standard deviation of aligned pair data
match_fasta.pl							
- script to generate ordered fasta files for global alignment
mauve_match_fasta.pl						
- script to generate ordered fasta files of anchor-matching pairs for global alignment
nlists_indel.pl							
- script to generate a nucleotide-level indel/mismatch list for aligned pairs
nlists_mauve.pl							
- script to generate a nucleotide-level indel/mismatch list for anchor-matching pairs
nscores_indel.pl						
- script to generate a nucleotide-level list of scores and locations for aligned pairs
nscores_mauve.pl						
- script to genearte a nucleotide-level list of scores and locations for anchor-matching pairs
order_check.pl							
- script to check how many protein sequences are read in the same/reverse order (no longer used)
separate.pl							
- script to parse sequential data of PI17 and PI25611
sort_blastp.pl							
- script to sort the BLASTp result by order of E value
subst_pg.pl							
- script to generate primer sequences of target PI proteins

