import math
import sys
from splice_object import splice_site

class seq_operator:
    def __init__(self, d1, d2, a1, a2):
        self.d1 = d1
        self.d2 = d2
        self.a1 = a1
        self.a2 = a2

    #Create list of alignments, where donor and acceptor sites are separated by tabs
    def alignments_init(self, scaffold_dict, introns):
        alignments = []
        for intron in introns:
            scaffold = scaffold_dict[intron.scaffold_id]
            scaffold_len = len(scaffold)
            splice_sites = ""
            #Find splice sites for positive strands
            if (intron.strand == "+"):
                if (((intron.donor - self.d1 + 1) < 0) or ((intron.acceptor + self.a2 + 1) > scaffold_len)):
                    print("\nERROR: please specify smaller donor/acceptor values\n")
                    sys.exit(0)
                splice_sites = scaffold[intron.donor - self.d1 : intron.donor + self.d2] \
                    + "\t" + scaffold[intron.acceptor - self.a1 : intron.acceptor + self.a2]
            #Find splice sites for negative strands
            else:
                if ((intron.donor - self.d2 < 0) or (intron.donor - self.a2 < 0) or (intron.acceptor + self.a1  > scaffold_len) or (intron.acceptor + self.d1 > scaffold_len)):
                    print("\nERROR: please specify smaller donor/acceptor values\n")
                    sys.exit(0)
                start = scaffold[intron.donor - self.d2 - 1 : intron.donor + self.d1 - 1]
                end = scaffold[intron.acceptor - self.a2 - 1: intron.acceptor + self.a1 - 1]
                #reverse transcribe
                splice_sites = self.reverse_transcribe(start, end)
            alignments.append(splice_sites)
        return alignments

    def reverse_transcribe(self, start, end):
        nucleotides = {"a" : "T", "c" : "G", "g" : "C", "t" : "A"}
        transcribed_sequence = ""
        new_sequence = start[::-1] + "\t" + end[::-1]
        for nucl in new_sequence:
            if(nucl == "\t"):
                transcribed_sequence += "\t"
            else:
                transcribed_sequence += nucleotides[nucl.lower()]
        return transcribed_sequence

    #CHECK
    def return_sites(self, alignments):
        sites = []
        lines = None
        for seq in alignments:
            lines = seq.split("\n")
            lines = lines[0].split("\t")
            sites.append(lines[0])

        for seq in alignments:
            lines = seq.split("\n")
            lines = lines[0].split("\t")
            sites.append(lines[1])
        return sites

    def generate_pfm(self, alignments, len_donors, len_acceptors):
        total_seq_length = len_donors + len_acceptors
        n_seq = len(alignments)
        n_alignments = int(len(alignments)/2)
        #Initialise 2D array
        pfm = []
        for r in range(4):
            row = []
            for c in range(total_seq_length):
                row.append(0)
            pfm.append(row)
        #Calculate PPM
        i = 0
        j = 0
        while i < n_seq:
            j = 0
            if (i < n_alignments):
                k = 0
            else:
                k = len_donors
            seq_length = len(alignments[i])
            while j < seq_length:
                nucl = alignments[i][j].lower()
                #Scoring donor sites
                if (nucl == 'a'):
                    pfm[0][k] += 1/n_alignments
                elif (nucl == 'c'):
                    pfm[1][k] += 1/n_alignments
                elif (nucl == 'g'):
                    pfm[2][k] += 1/n_alignments
                elif (nucl == 't'):
                    pfm[3][k] += 1/n_alignments
                j += 1
                k += 1
            i += 1
        print(pfm)
        return pfm

    def generate_pwm(self, size, pfm, a_freq, c_freq, g_freq, t_freq):
        #background model and corresponding probabilities: A, C, G, T
        b_k = [a_freq, c_freq, g_freq, t_freq]
        i = 0
        j = 0
        b = 0
        while i < len(pfm):
            j = 0
            while j < len(pfm[i]):
                m = math.log((pfm[i][j] + (1/size) * (1/size)) / b_k[b], 2)
                pfm[i][j] = m
                j += 1
            b += 1
            i += 1
        return pfm

    def scan_donor_sites(self, sequence_dict, pwm, max_score, min_score):
        #max score possible
        max_score = max_score

        #min score possible
        min_score = min_score

        donor_sites = []
        #Calculate donor sequence with highest score
        for keys in sequence_dict:
            i = self.d1
            current_score = 0
            previous_score = 0
            donor_index = 0
            current_sequence = sequence_dict[keys]
            length = len(current_sequence)
            while i < length - self.d2:
                current_score = 0
                j = 0
                donor_index = i
                for c in current_sequence[i - self.d1 : i + self.d2]:
                    nucl = c.lower()
                    if nucl == "a":
                        current_score += pwm[0][j]
                    elif nucl == "c":
                        current_score += pwm[1][j]
                    elif nucl == "g":
                        current_score += pwm[2][j]
                    elif nucl == "t":
                        current_score += pwm[3][j]
                    else:
                        print("\n ERROR: sequence has unknown nucleotide - sequence must only contains A, C, G or T\n")
                        sys.exit(0)
                    j += 1
                rel_score = ((current_score - min_score)/(max_score - min_score)) * 100
                if (rel_score >= 75):
                    donor_site = splice_site(keys, sequence_dict[keys], donor_index, current_score, rel_score)
                    donor_sites.append(donor_site)
                i += 1
        #donor_sites.reverse()
        return donor_sites

    #CHECK
    def scan_acceptor_sites(self, sequence_dict, pwm, max_score, min_score):
        #max score possible
        max_score = max_score
        #min score possible
        min_score = min_score
        acceptor_sites = []
        #Calculate donor sequence with highest score
        for keys in sequence_dict:
            i = self.a1
            donor_length = self.d1 + self.d2
            current_score = 0
            previous_score = 0
            acceptor_index = 0
            current_sequence = sequence_dict[keys]
            length = len(current_sequence)
            while i < length - self.a2:
                current_score = 0
                j = donor_length
                acceptor_index = i
                for c in current_sequence[i - self.a1 : i + self.a2]:
                    nucl = c.lower()
                    if nucl == "a":
                        current_score += pwm[0][j]
                    elif nucl == "c":
                        current_score += pwm[1][j]
                    elif nucl == "g":
                        current_score += pwm[2][j]
                    elif nucl == "t":
                        current_score += pwm[3][j]
                    else:
                        print("\n ERROR: sequence has unknown nucleotide - sequence must only contains A, C, G or T\n")
                        sys.exit(0)
                    j += 1
                rel_score = ((current_score - min_score)/(max_score - min_score)) * 100
                if (current_score >= 75):
                    acceptor_site = splice_site(keys, sequence_dict[keys], acceptor_index, current_score, rel_score)
                    acceptor_sites.append(acceptor_site)
                i += 1
        return acceptor_sites


    def matrix_donor_max(self, pwm):
        i = 0
        max_score = 0
        length = self.d1 + self.d2
        while i < length:
            max = 0
            j = 0
            while j < 4:
                current = pwm[j][i]
                if current > max:
                    max = current
                j += 1
            max_score += max
            i += 1
        return max_score

    def matrix_donor_min(self, pwm):
        i = 0
        min_score = 0
        length = self.d1 + self.d2
        while i < length:
            j = 0
            min = pwm[j][i]
            while j < 4:
                current = pwm[j][i]
                if current < min:
                    min = current
                j += 1
            min_score += min
            i += 1
        return min_score

    def matrix_acceptor_max(self, pwm):
        max_score = 0
        i = self.d1 + self.d2
        total_length = len(pwm[0])
        while i < total_length:
            j = 0
            max = pwm[j][i]
            while j < 4:
                current = pwm[j][i]
                if current > max:
                    max = current
                j += 1
            max_score += max
            i += 1
        return max_score

    def matrix_acceptor_min(self, pwm):
        min_score = 0
        i = self.d1 + self.d2
        total_length = len(pwm[0])
        while i < total_length:
            j = 0
            min = pwm[j][i]
            while j < 4:
                current = pwm[j][i]
                if current < min:
                    min = current
                j += 1
            min_score += min
            i += 1
        return min_score
