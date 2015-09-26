#Overview
This project was initiated by Prof. Jae Hyung Lee to analyze a lesser known strain of Prevotella Intermedia and compare it with a more well-known strain to recognize its differences and relative functions.
##Prevotella Intermedia 
Prevotella Intermedia is known to cause various periodontal diseases, including periodontitis and gingivitis. Despite its important role, this microorganism is not yet heavily studied, and only one of its strains, P. Intermedia 17, has its complete DNA sequences assembled. Meanwhile, a less known strain, ATCC 25611, is assembled into 24 DNA contigs, with its coding sequences not fully known. We intend to compare this less known strain, 25611, to a well known strain, 17, thereby further understanding the functional differences between the two strains and ultimately building a global alignment of strain 25611 DNA sequences compared to strain 17.

#Methods
To compare a DNA sequence with an imperfect one, we used BLAST to compare the proteins of Prevotella Intermedia 17 (PI 17) with open reading frames(ORFs) of Prevotella Intermedia 25611 (PI 25611), and MAUVE to compare the two chromosomes of PI 17 and 24 contains of PI 25611. Once we retrieved the initial alignment data, we classified the open reading frames according to similarities with already known protein sequences. ORFs with matching protein sequences are more likely to behave in a similar way with already known proteins, and therefore would serve as important information for comparing the two DNA sequences. Single Nucleotidal Polymorphisms (SNPs) and indels (insertion/deletions) are also calculated according to this data.

#More Descriptions
This page does not list the details of every file included in this repository or the original project. For more detailed descriptions on files and folders in the original workspace and how they work in the original context, see readme.txt.

**Note**: Most of the following scripts are merged into *classify.pl*. The bulk of the project is done by classify.pl, while other scripts are supplementary for processing data from other tools, statistical interpretations, or reformatting of data.

Also note that data files, including DNA sequences, parsed data, and protein sequence data are missing in this repository for storage capacity.