# MinorityReport
Software for generalized analysis of causal genetic variants
written by Jeremy A. Horst

MinorityReport.py is a python script meant to find genetic differences in parent-child diads or strain pairs by comparing genomic sequencing reads aligned to the same reference genome. 

This program takes the genome FASTA file, the corresponding gene model in GFF3 format, and the parent and child SAM files, respectively. FASTQ files should be filtered for quality, aligned to a reference genome with a tool such as BowTie2, and output in SAM format. 

The script MinorityReport-MASTER.py is meant to divide this task across processors on your computer, one per chromosome in the genome. This accelerates the process enormously, but requires one processor per chromosome.


Usage: MinorityReport.py <ref seq FASTA> <ref seq gff> <sam alignment parent> <sam alignment mutant>

e.g.:    ./MinorityReport.py species.fna species.gff species-bug1.bowtie2.sam species-bug2.bowtie2.sam
Options: 
         
         -vp  <minimum_variant_proportion>	    minimum fraction of reads at position supporting variant to accept (default=0.3).
         -wp  <maximum_variant_proportion>	    maximum fraction of reads at position supporting wildtype to accept (default=0.01).
         -vc  <minimum_variant_counts>	        minimum mutant-strain reads covering position to evaluate (default=30).
         -wc  <maximum_wildtype_variant_counts>	maximum parent-strain reads with variant to report (default=0).
         -wtc <minimum_wildtype_total_counts>	minimum total parent-strain reads covering position to report (default=0).

         -cnv     analyze Copy Number Variants. Also available through CNV_caller.py
         -median  calculate median instead of mean CNV for tile range of positions (default: mean).
         -window_size       <sliding window length>  cnv: size of window to check for each position  (default=3000).
         -window_increment  <sliding window spacing> cnv: increment for progressing through chromosome(s) (default=3000).
         -report_frequency  <report window spacing>  cnv: increment for reporting - must be multiple of window_increment (default=3000).
         -o_cnv   <output_file> specify output file for CNV analysis (default: print to stdout).

         -pe  paired read matching alignments ONLY (default=False).
         -gene_type  <GFF gene type>		descriptor to find correct entry lines in GFF file.
         -cds evaluate CDS for start and stop codons in each CDS, i.e. don't trust GFF file (default=False).
         -all report all variants whether nonsynonymous or only nucleotide changes (default=False).
         -position_read_report <output_file> output a file with the number of counts that map to each position in the genome for parent & mutant (default=False).
         -o   <output_file> specify output file (default: print to stdout).

ABSTRACT

	Background: The widespread availability of next generation genome sequencing
	technologies has enabled a wide range of variant detection applications,
	especially in cancer and inborn genetic disorders. For model systems and
	microorganisms, the same technology may be used to discover the causative
	mutations for any phenotype, including those generated in response to chemical
	perturbation. In the case of pathogenic organisms, these approaches have allowed
	the determination of drug targets by means of resistance selection followed by
	genome sequencing.
	
	Results: Here, we present open source software written in python,
	MinorityReport, to facilitate the comparison of any two sets of genome
	alignments for the purpose of rapidly identifying the spectrum of nonsynonymous
	changes, insertions or deletions, and copy number variations in a presumed
	mutant relative to its parent. Specifically, MinorityReport relates mapped
	sequence reads in SAM format output from any alignment tool for both the mutant
	and parent genome, relative to a reference genome, and produces the set of
	variants that distinguishes the mutant from the parent, all presented in an
	intuitive, straightforward report format. MinorityReport features tunable
	parameters for evaluating evidence and a scoring system that prioritizes
	reported variants based on relative proportions of read counts supporting the
	variant in the mutant versus parent data sets. We demonstrate the utility of
	MinorityReport using publicly available data sets that we previously published
	to find the determinants of resistance for novel anti-malarial drugs.
	
	Conclusions: MinorityReport is readily available (github: xxxxxxx) to identify
	the genetic mechanisms of drug resistance in plasmodium, genotype-phenotype
	relationships in human diads, or genomic variations between any two related
	organisms.

MinorityReport.py
20150826 Jeremy Horst
20161128 last update (added SoftClip option in CIGAR string)

	INPUT:	FASTA of reference sequence
			GFF gene annotation file (SAME CHROMOSOME NAMES AS FASTA!!!)
			parent strain SAM alignment file
			child/mutant strain SAM alignment file
	
	nonsynonymous variants
	PROCESS: *populate hash table for each position from SAM files, with evidence = counts
			  handle deletions as #-'s. 
			  handle insertions as sequence position = A[CTGCTGCTG], where A is the reference sequence position left of the insertion
			  e.g. sequence[2456821] = {'A':3,'---':426,'C':5,'T':3,'G':51, 'GTCGTACGTAGCTAGC':20}
			  check for variants over % threshold
			  check if variant is in protein coding region
			  check if variant changes the amino acid
			  compare parent to child
	OUTPUT:   report all variants over n% abundant in reads mapped to this position, unique in child w.r.t. parent
	
	copy number variants
	PROCESS:  build gene model & data structure from FASTA & GFF
			 *read in each SAM file, with the start & end of each read
			  examine tiles / sliding windows (depending on settings) across each chromosome
			 *check overlap of reads to each tile
			  take ratio of mutant to parent, in context of overall counts  
			  do statistics on distribution for read ratio for each tile
			 * = slow steps / bottlenecks
	OUTPUT:	report statistics of read ratios for each tile
