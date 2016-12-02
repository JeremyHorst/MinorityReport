#test-readme.txt

1. download test files from here (127MB):
	<https://drive.google.com/open?id=0B7G5rgobx2cDNjNaVUxESnFTY00>

2. in a command line shell, unpack the test tarball as follows: 
	% tar -xzf test.tgz 
	(this will take a minute)

3. cd to "test" directory
	also find the path to the MinorityReport.py script file. Swap this out for [path] in the next step.

4. run test as follows:
	% python [path]/MinorityReport.py PlasmoDB-3D7_9.3-chr12.fasta Pfgenes-3D7_9.3-chr12.gff 3D7Parent-100Q20-chr12.sam s83-100Q20-chr12.sam
	
	use standard linux notation to save the output through standard output and standard error redirects.
	(This takes 5-6 minutes on a 2011 MacBook Pro while lots of other processes are running.)
	
	The standard output should look like:
	
	#score	chromosome	genomic_mutation	gene_mutation	protein_mutation	mutation_counts	mutation_proportion	mutant_position_countsparent_mutation_counts	parent_mutation_proportion	parent_position_counts	|	gene_id	strand	protein_length	gene_description	||	A-gene_variant	protein_variant	variant_counts	variant_proportion	|	parent_variant_counts	parent_variant_proportion	||	C-gene_variant	protein_variant	variant_counts	variant_proportion	|	parent_variant_counts	parent_variant_proportion	||	G-gene_variantprotein_variant	variant_counts	variant_proportion	|	parent_variant_counts	parent_variant_proportion	||	T-gene_variant	protein_variant	variant_counts	variant_proportion	|	parent_variant_counts	parent_variant_proportion	||	nonACGT-gene_variant	protein_variant	variant_counts	variant_proportion	|	parent_variant_counts	parent_variant_proportion	||
	4.295	chr12	G531566T	C1233A	P412T	72	0.973	74	|	1	0.01	104.0	|	PF3D7_1211900-1	-	1264	non-SERCA-type Ca2  -transporting P-ATPase (ATP4)	||	C1233T	P412S	0	(0.0)	|	1	(0.01)	||	C1233G	P412A	0	(0.0)	(0.0)	||	C1233C	-	1	(0.014)	|	102	(0.981)	||	C1233A	P412T	72	(0.973)	|	1	(0.01)	||

