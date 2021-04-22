SPLICER: Generalised Splice-site prediction program.

Dependency: pysam (Splicer requires pysam as a dependency, please install pysam)

Splicer will predict splice sites for query sequences using a Position Weight Matrix (PWM) generated from organism specific alignments.

Before reading the detailed description of how to run the program I will outline the basic steps below:

1) Create a bedfile containing known donor and acceptor locations: these sites will be used to generate the organism-specific PWM.
2) Run Splicer by passing in the bedfile from step 1, including the file containing the scaffolds (fasta file) to generate the alignment and matrix files (the matrix file will be used to scan query sequences)
3) Run Splicer and pass in query sequences (fasta file)

STEP 1:

In order for Splicer to generate an organism specific PWM it will require alignments of known donor and acceptor sites from your organism of study. You can do this by aligning transcripts from transcriptome data to the genome and recording the start and ends locations of exons. This information will be passed into Splicer as a bedfile: a tab-delimited file with the extension .bed.

Quick reminder:

---------------<<<< EXON 1 >>>>-----------------------------<<<< EXON 2 >>>>-------------------
															|                             |
												  donor site                    acceptor site										                               

Example:

chrX1	1000	1100	Exon1	0	+
chrX1	1150	1170	Exon2	0	+
chrX2	1000	1120	Exon3	0	-
chrX2	1200	1300	Exon2	0	-
chrX2	1200	1300	Exon2	0	-
chrX2	1400	1450	Exon1	0	-

Column 1 contains the scaffold ID of your gene with known donor and acceptor sites, in this example, the first two exons are located on the scaffold chrX1. Column 1 contains the start of the exon and columns 2 contains end of the exon (reverse if they ar on the "-" strand). The next two columns are not important to Splicer so you can put in any value you'd like. The final column contains the direction of the strand.

Important Notes: Just make sure that the locations are in <ascending order>, regardless of the strand. If you do not do this  Splicer will incorrectly swap donor and acceptor sites.

One final note is that there must be an even amount of lines in the bedfile for Splicer to correctly identify the known donor and acceptor sites. Notice that the gene on the "-" strand has 3 exons. In order to keep the lines even I have copied the details of Exon2 twice, this is because Splicer will look at only one column for each line to identify a donor or acceptor site like so:

chrX2	1000	1120	Exon3	0	-    <---- column 2: acceptor site
chrX2	1200	1300	Exon2	0	-    <---- column 1: donor site
chrX2	1200	1300	Exon2	0	-    <---- column 2: acceptor site
chrX2	1400	1450	Exon1	0	-    <---- column 1: donor site


STEP 2:

Splice will also need the genome fasta file, alternatively, you can create a fasta file containing all the scaffolds listed in your bedfile (just like you would using samtools faidx). The scaffolds specified in the bedfile from step 1 will be used to search the genome/scaffold files so it important that the scaffold IDs in the bedfile and fasta file are <identical>.

Now you can run Splicer by passing in the file: bedfile and fasta genome/scaffold file. There are multiple flags which can be used but we will stick with -g for now.

Run Splicer with the following commandline arguments:

$python3 splicer.py <flag> <file_containing_scaffolds> <bed_file> <output_matrix_file> <d1> <d2> <a1> <a2>

$python3 splicer.py -g scaffolds.fasta bed_file.bed output-matrix.txt 5 3 5 3

The commandline arguments d1, d2, a1, and a2 refer to how many bases before and after the donor and acceptor sites you'd like to inspect, for example:

---------------<<<< EXON 1 >>>>-----------------------------<<<< EXON 2 >>>>-------------------
															|                             |
												  donor site                    acceptor site
													d1      d2                    a1         a2

by assigned d1 = 5, d2 = 3, a1 = 5, and a2 = 12 you are saying that you'd like Splicer to inspect 5 bases before the donor site (d1) and 3 bases after the donor site (d2). The same logic applies to the acceptor site. d1, d2, a1, a2 must all have values greater than zero. Please remember the values as the same values will need to be specified when scanning your query sequences.

If you know the frequency of seeing each nucleotide in your genome (background models) you can specify them as well (as shown below). Splicer will assume an equal frequency of each nucleotide (0.25) if you do not explicitly specify them:

$python $python3 splicer.py <flag> <file_containing_scaffolds> <bed_file> <output_matrix_file> <d1> <d2> <a1> <a2> <freq_a> <feq_c> <freq_g> <freq_t>

Two files will be generated: An alignment file and the PWM text file. The alignment file is there in case you'd like to inspect the splice site sequences pulled from the scaffolds. Like a bedfile, it is tab-delimited:

AGGGAAAA	TAATTTTT
AGAAGAAA	TTTCTTTT

The first column contains the known donor splice sites and the second column contains the known acceptor sites pulled from the scaffold file. Since we specified that we'd like to inspect 5 bases before and 3 bases after the splice sites, the donor and acceptors sites each have a total length of 8. Notice that all sequences in column 1 are equal in length to each other. The same applies for column 2. However, sequences between the two columns do not need to equal in length. If we had specified d1=5, d2=3, a1=2, a2=2 then we'd have:

AGGGAAAA	TTTT
AGAAGAAA	CTTT

The other file generated is the PWM, this is a text file containing the values used to calculate scores for your query sequences. Please do not modify this file in any way.

STEP 3:

You can now run Splicer to scan your query sequences by passing in the PWM generated in the previous step:

$python3 splicer.py <flag> <query_sequences> <matrix> <output_name> <d1> <d2> <a1> <a2>

$python3 splicer.py -s query-sequences.fasta output-matrix.txt results.txt 5 3 5 3

Note: You must use the same d1, d2, a1, and a2 parameters used when you generated the matrix file.

OUTPUT:

You will now have the results as a text file.

Donor sites: For each sequence, Splicer will output the query sequence from the start of the sequence up until it reaches a putative donor splice site. You will be able to see this as the sequence will end with ">". Each query sequence will also have the relative score (%), this just tells the you how close to the donor site score got to the maximum score possible. Thus, the higher the percentage the closer it resembles the known donor splice sites from your organism.

Acceptor sites: For each sequence, Splicer will output the query sequence from first putative acceptor splice site until the end of the query sequence. Thus, putative acceptor sites will start with "<"

OTHER FLAGS:

You can skip steps 1-2 if you create your own alignment file in the same format specified in step 2. For example, you might already have known donor and acceptor sequences, in this case you can create an alignment text file (tab-delimited) and use the flag: -a

$python3 splicer.py <flag> <alignments> <output_matrix> <d1> <d2> <a1> <a2>

Like the -g flag, you can specify the background models for each of the nucleotides. Splicer will assume equal probability (0.25) for each nucleotide if you don't.

$python3 splicer.py <flag> <alignments> <output_matrix> <d1> <d2> <a1> <a2> <freq_a> <feq_c> <freq_g> <freq_t>
