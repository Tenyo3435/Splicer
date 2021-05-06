import sys
from file_handler import file_op
from sequence_operator import seq_operator
#COMMAND LINE ARGUMENTS
# <scaffolds> lists the file which contains the scaffolds to be searched, can be fasta genome or samtools generated file.
# <bed> is a bed file which contains all the putative exons and their start/end locations. This should have been generated
#by aligning transcriptome data to the genome in order to identity putative exons and introns.
# <output> is the name of the output file you specify
# <donor_1> is the number of letters before the donor splice site you'd like to isolate.
# <donor_2> is the number of letters after the donor splice site you'd like to isolate.
# <acceptor_1> is the number of letter before the acceptor splice site you'd like to isolate.
# <acceptor_2> is the number of letters after the acceptor splice site you'd like to isolate.
try:
    flag = sys.argv[1]
except:
    print("\nERROR: flag does not exist: -g or -a or -s\n")

#############################################################
# FLAG -g: GENERATE PWM FROM SCAFFOLD AND .BED FILE##########
#############################################################
if (flag == "-g"):
    scaffolds = sys.argv[2]
    bed = sys.argv[3]
    output = sys.argv[4]
    donor_1 = int(sys.argv[5])
    donor_2 = int(sys.argv[6])
    acceptor_1 = int(sys.argv[7])
    acceptor_2 = int(sys.argv[8])

    if(donor_1 < 0 or donor_2 < 0 or acceptor_1 < 0 or acceptor_2 < 0):
        print("\nERROR: d1, d2, a1, a2 must be >0")
        sys.exit(0)

    #Initialise file operator object to work with fasta and bed files
    file_operator = file_op()

    #Create intron objects from bedfile using file_operator
    introns = file_operator.introns_list_init(bed)

    #Create dictionary that contains all the scaffolds from the scaffold file. Where key = scaffold ID and value = corresponding sequence
    scaffold_dict = file_operator.scaffold_dict_init(scaffolds, introns)

    #Instantiate seq_operator to create alignments and pwm
    operator = seq_operator(donor_1, donor_2, acceptor_1, acceptor_2)

    #Create alignments
    alignments = operator.alignments_init(scaffold_dict, introns)

    #Write alignments to file for later reference
    file_operator.write_file(alignments, "alignments.txt")

    #Return list containing donor and acceptor sites from alignments generated
    site_alignments = operator.return_sites(alignments)
    size = len(site_alignments)

    #Calculate PFM
    pfm = operator.generate_pfm(site_alignments, donor_1 + donor_2, acceptor_1 + acceptor_2)

    pwm = None
    #Calculate default PWM
    if (len(sys.argv) == 9):

    #Calculate PWM with custom background probabilties
        pwm = operator.generate_pwm(size, pfm, 0.25, 0.25, 0.25, 0.25)
    elif (len(sys.argv) == 13):
        try:
            a_freq = float(sys.argv[9])
            c_freq = float(sys.argv[10])
            g_freq = float(sys.argv[11])
            t_freq = float(sys.argv[12])
        except:
            print("\nERROR: background probabilities must be decimals numbers <= 1.0\n")


        pwm = operator.generate_pwm(size, pfm, a_freq, c_freq, g_freq, t_freq)
    else:
        print("\nERROR: you must specify 4 background probabilities in the following order as decimals A C G T\n")
        sys.exit(0)

    file_operator.output_matrix(pwm, output)

#############################################################
# FLAG -a: GENERATE PWM FROM ALIGNMENT FILE #################
#############################################################

elif (flag == "-a"):
    #Alignment file
    alignment_file = sys.argv[2]
    output = sys.argv[3]
    donor_1 = int(sys.argv[4])
    donor_2 = int(sys.argv[5])
    acceptor_1 = int(sys.argv[6])
    acceptor_2 = int(sys.argv[7])

    if(donor_1 < 0 or donor_2 < 0 or acceptor_1 < 0 or acceptor_2 < 0):
        print("\nERROR: d1, d2, a1, a2 must be >0")
        sys.exit(0)

    #Initialise file operator object to work with fasta and bed files
    file_operator = file_op()

    #Initialise sequence operator object to do sequence calculations
    operator = seq_operator(donor_1, donor_2, acceptor_1, acceptor_2)

    #Generate alignments
    alignments = file_operator.alignments_init(alignment_file)
    size = len(alignments)
    #Calculate PFM
    pfm = operator.generate_pfm(alignments, donor_1 + donor_2, acceptor_1 + acceptor_2)

    #Calcualte PWM
    pwm = None

    #Calculate PWM with default values
    if (len(sys.argv) == 8):

        pwm = operator.generate_pwm(size, pfm, 0.25, 0.25, 0.25, 0.25)
    #Calculate PWM with custom background model values
    elif (len(sys.argv) == 12):
        try:
            a_freq = float(sys.argv[9])
            c_freq = float(sys.argv[10])
            g_freq = float(sys.argv[11])
            t_freq = float(sys.argv[12])
        except:
            print("\nERROR: background probabilities must be decimals numbers <= 1.0\n")

        pwm = operator.generate_pwm(size, pfm, a_freq, c_freq, g_freq, t_freq)

    file_operator.output_matrix(pwm, output)

#############################################################
# FLAG -s: SCAN QUERY SEQUENCES FOR SPLICE SITES ############
#############################################################

elif (flag == "-s"):
    sequences = sys.argv[2]
    matrix_file = sys.argv[3]
    output = sys.argv[4]
    donor_1 = int(sys.argv[5])
    donor_2 = int(sys.argv[6])
    acceptor_1 = int(sys.argv[7])
    acceptor_2 = int(sys.argv[8])

    if(donor_1 < 0 or donor_2 < 0 or acceptor_1 < 0 or acceptor_2 < 0):
        print("\nERROR: d1, d2, a1, a2 must be >0")
        sys.exit(0)

    #Instantiate file operator
    file_operator = file_op()

    #Open sequences and create dictionary
    sequence_dict = file_operator.sequence_dict_init(sequences)

    #Generate matrices from reading matrix file
    pwm = file_operator.read_matrix(matrix_file)

    #Instantiate sequence operator
    operator = seq_operator(donor_1, donor_2, acceptor_1, acceptor_2)

    #Calculate max and min scores of pwm
    max_donor_score = operator.matrix_donor_max(pwm)
    min_donor_score = operator.matrix_donor_min(pwm)
    #Scan sequences and generate donor results
    donor_results = operator.scan_donor_sites(sequence_dict, pwm, max_donor_score, min_donor_score)

    #Calculate max and min scores of pwm
    max_acceptor_score = operator.matrix_acceptor_max(pwm)
    min_acceptor_score = operator.matrix_acceptor_min(pwm)

    #Scan sequences and generate acceptor results
    acceptor_results = operator.scan_acceptor_sites(sequence_dict, pwm, max_acceptor_score, min_acceptor_score)

    #Output results from scanning sequences
    file_operator.output_results(donor_results, acceptor_results, output)

else:
    print("\nERROR: Please check that you are using the correct flag and command line arguments\n")
