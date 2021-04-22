from pysam import FastaFile
from intron_object import intron
import sys

class file_op:
    #Read file and create dictionary with scaffold IDs and corresponding sequences
    #All new line characters have been removed in order to get accurate length of sequences
    def scaffold_dict_init(self, file, introns):
        sequences_object = FastaFile(file)
        dict = {}
        for i in introns:
            if (not(i.scaffold_id in dict)):
                dict[i.scaffold_id] = sequences_object.fetch(i.scaffold_id)
        return dict
    
    def sequence_dict_init(self, file):
        dict = {}
        key = None
        try:
            with open(file) as seq_file:
                for line in seq_file:
                    if (line[0] == ">"):
                        line = line.split("\n")
                        line = line[0].split(">")
                        key = line[1]
                        dict[key] = ""
                    else:
                        line = line.split("\n")
                        dict[key] += line[0]
        except: 
            print("\nERROR: file does not exist or is not fasta format\n")
        return dict
    
    def alignments_init(self, file):
        alignments = []
        acceptors = []
        try:
            with open(file) as f:
                for line in f:
                    line = line.split("\n")
                    line = line[0].split("\t")
                    alignments.append(line[0])
                    acceptors.append(line[1])
                for a in acceptors:
                    alignments.append(a)
        except:
            print("\nERROR: Alignment file does not exist or is not correctly formatted\n")
        return alignments
                

    #initialise and return list of intron objects 
    def introns_list_init(self, bed_file):
        introns = []
        try:
            bed_file = open(bed_file, "r")
        except:
            print("\nERROR: bed file does not exist\n")
            sys.exit(0)
        i = 0
        start = 0
        end = 0
        try:
            for line in bed_file:
                bed_line = line.split("\n")
                bed_line = bed_line[0].split("\t")
                if (bed_line[5] == "+"):
                    if(i % 2 == 0):
                        start = int(bed_line[2])
                    else:
                        end = int(bed_line[1])
                #Can cut this out, wait to see what everyone votes
                elif (bed_line[5] == "-"):
                    if(i % 2 == 0):
                        end = int(bed_line[2])
                    else:
                        start = int(bed_line[1])
                if(i % 2 != 0):
                    new_intron = intron(bed_line[0], start, end, bed_line[5])
                    introns.append(new_intron)
                i += 1
        except:
            print("\nERROR: bedfile is not tab-delimited\n")
        return introns

    def write_file(self, list, output):
        file = open(output, "w")
        for line in list:
            file.write(line + "\n")
        file.close()
    
    def output_matrix(self, matrix, output):
        output_string = ""
        i = 0
        j = 0
        while i < len(matrix):
            j = 0
            while j < len(matrix[i]):
                if(j == len(matrix[i]) - 1):
                    output_string += str(matrix[i][j]) + "\n"
                else:
                    output_string += str(matrix[i][j]) + " "
                j += 1
            i += 1
        file = open(output, "w")
        file.write(output_string)
        file.close()

    def read_matrix(self, matrix_file):
        try:
            file = open(matrix_file, "r")
        except:
            print("\nERROR: matrix file does not exist\n")
        matrix = []
        try:
            for line in file:
                line = line.split("\n")
                line = line[0].split(" ")
                row = []
                for n in line:
                    row.append(float(n))
                matrix.append(row)
        except:
            print("\nERROR: Matrix file is not correctly formatted, please seperate donor and acceptor matrices with tabs and values by \s\n")
        return matrix
    
    def output_results(self, donor_results, acceptor_results, output):
        output_file = open(output, "w")
        output_string = "Splicer Results v1.0\n\n                           =================================\n                           ========== DONOR SITES ==========\n                           =================================\n\n"
        for splice_site in donor_results:
            output_string += "\nID:" + splice_site.scaffold_id + " Score:{:.2f}".format(splice_site.score) + " Relative Score (%):{:.2f}".format(splice_site.relative_score) + " Index:" +\
                str(splice_site.index) + "\n" + "SEQUENCE: \n" + splice_site.sequence[0 : splice_site.index] + ">\n"
        output_string += "\n\n\n\n                           =================================\n                           ======== ACCEPTOR SITES =========\n                           =================================\n\n"

        for splice_site in acceptor_results:
            output_string += "\nID:" + splice_site.scaffold_id + " Score:{:.2f}".format(splice_site.score) + " Relative Score (%):{:.2f}".format(splice_site.relative_score) + " Index:" +\
                str(splice_site.index) + "\n" + "SEQUENCE: \n" + "<" + splice_site.sequence[splice_site.index:] + "\n"
        output_file.write(output_string)

