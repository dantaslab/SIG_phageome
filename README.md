# SIG_phageome
**Phage ecology drives niche-specific adaptation of *Staphylococcus pseudintermedius* in clinical diagnostic settings**

## Phage analysis pipeline among *Sp* genomes
**Phage identification**
1. Assemblies were piped through Cenote-Taker 3 to identify putative phage contigs with end features as direct terminal repeats (DTRs) indicating circularity, and inverted linear repeats (ITRs) or no features for linear sequences. Run 1a first and then use 1b to compile the results from different samples. (edited)
2. Contigs in each cluster were concatenated and filtered by length and completeness to remove false positives. Specifically, the length limits were 1,000 nt for the detection of circularity, 4,000 nt for ITRs, and 5,000 nt for other linear sequences.
3. The completeness was computed as a ratio between the length of our phage sequence and the length of matched reference genomes by CheckV and the threshold was set to 10.0%.
4. Phage contigs passed these two filters were then run through VIBRANT with “virome” flag to further remove obvious non-viral sequences.

**Phage taxonomy assignment**
* Protein sequences created by CheckV were used as input for vConTACT2 with “DIAMOND” and database “ProkaryoticViralRefSeq207-Merged” to assign taxonomy. For the “unsigned” ones from vConTACT2, we used the tentative taxonomy from Cenote-Taker 3 inferred using BLASTP against a custom database containing Refseq virus and plasmid sequences from GenBank.

**Phage sequencing clustering**
* Based on MIUViG recommended parameters, phages were grouped into populations if they shared ≥95% nucleotide identity across ≥85% of the genome using BLASTN and a CheckV supporting code, anicalc.py (https://bitbucket.org/berkeleylab/checkv/src/master/).

**Gene annotation**
* Phage-encoded metabolic genes will be annotated using VIBRANT. Perform Pharokka for phage gene functional annotation.
