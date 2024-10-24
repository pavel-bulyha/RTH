import numpy as np
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

def palindrome_finder(sequence, score=None, calc="no"):
    n_rut = 0
    list_scores = []
    sequence = Seq(sequence)
    sequence_r = str(sequence.reverse_complement())
    _pattern_pause_site1 = r"GG\D{8}[C,T]G"
    _pattern_pause_site2 = r"C[G,A]\D{8}CC"
    for k in range(4, 8):
        for i in range(150):
            _x1 = i
            _y1 = i + k
            if _y1 > len(sequence):
                break
            else:
                window = sequence[_x1:_y1]
                window = str(window)
                gc_content = 100 * (window.count("C") + window.count("G")) / len(window)
                score_p = score
                sub_sequence_r = sequence_r
                _start = 0
                while True:
                    if re.search(window, sub_sequence_r):
                        _research = re.search(window, sub_sequence_r)
                        _positions = _research.span()
                        _x2 = _positions[0]
                        _y2 = _positions[1]
                        sub_sequence_r = sub_sequence_r[_y2:]
                        _x2 += _start
                        _y2 += _start
                        _start += _y2
                        if 4 <= len(sequence) - _y2 - _y1 + 1 <= 8:
                            if calc == "yes":
                                loop = len(sequence) - _y2 - _y1
                                score_p += 3
                                if gc_content > 20:
                                    score_p += 2
                                elif gc_content > 10:
                                    score_p += 1
                                if len(window) > 4:
                                    score_p += 1
                                if loop < 6:
                                    score_p += 1
                                if score_p > 0:
                                    list_scores.append(score_p)
                            else:
                                n_rut += 1
                    else:
                        break
    if calc == "yes":
        return np.max(list_scores) if list_scores else 0
    elif calc == "no":
        return n_rut

sequences_to_analyze = {}

print("RhoTermPredict is a genome-wide predictor of transcription Rho-dependent terminators "
      "in bacterial genomes. It analyzes both the strands.\n\n")
print("Input: Genome sequences file")
print("Output: a GenBank file containing Rho-dependent terminators annotations")
print("Genome file must be in genbank format")
file_genome = input("Enter the input genome file name: ")

try:
    print("Starting the search for Rho-dependent terminators...")
    for seq_record in SeqIO.parse(file_genome, "genbank"):
        genome = str(seq_record.seq)
        genome = genome.upper()
        gc_whole_genome = 100 * (genome.count("G") + genome.count("C")) / len(genome)
        pattern1 = r"C\D{11,13}C\D{11,13}C\D{11,13}C\D{11,13}C\D{11,13}C"
        pattern2 = r"G\D{11,13}G\D{11,13}G\D{11,13}G\D{11,13}G\D{11,13}G"
        pattern_pause_site1 = r"GG\D{8}[C,T]G"
        pattern_pause_site2 = r"C[G,A]\D{8}CC"
        found_terminators = []  # To store found terminators

        # positive strand
        scale = 0
        for j in range(len(genome)):
            x1 = scale + j
            x2 = scale + j + 78
            if x2 > len(genome):
                break
            else:
                w = genome[x1:x2]
                if w.count("G") > 0:
                    numG = w.count("G")
                    c_over_g = w.count("C") / numG
                else:
                    numG = 1
                    c_over_g = (w.count("C") + 1) / numG
                if c_over_g > 1:
                    if re.search(pattern1, w):
                        data = []
                        for g in range(50):
                            if x2 + g <= len(genome):
                                w = genome[x1 + g:x2 + g]
                                if w.count("G") > 0:
                                    numG = w.count("G")
                                    c_over_g = w.count("C") / numG
                                else:
                                    numG = 1
                                    c_over_g = (w.count("C") + 1) / numG
                                if c_over_g > 1:
                                    if re.search(pattern1, w):
                                        data.append(c_over_g)
                                    else:
                                        data.append(0)
                                else:
                                    data.append(0)
                            else:
                                break
                        maxP = np.argmax(data)
                        c_over_g = np.max(data)
                        x1 = x1 + maxP
                        x2 = x2 + maxP
                        scale = x2 - j - 1
                        s = genome[x2:x2 + 150]
                        score = 3
                        ctrl = palindrome_finder(s)
                        if ctrl > 0 or re.search(pattern_pause_site1, s):
                            found_terminators.append((int(x1), int(x2)))  # Store found terminators

        # negative strand
        scale = 0
        for j in range(len(genome)):
            x1 = scale + j
            x2 = scale + j + 78
            if x2 > len(genome):
                break
            else:
                w = genome[x1:x2]
                if w.count("C") > 0:
                    numC = w.count("C")
                    c_over_g = w.count("G") / numC
                else:
                    numC = 1
                    c_over_g = (w.count("G") + 1) / numC
                if c_over_g > 1:
                    if re.search(pattern2, w):
                        data = []
                        for g in range(50):
                            if x2 + g <= len(genome):
                                w = genome[x1 + g:x2 + g]
                                if w.count("C") > 0:
                                    numC = w.count("C")
                                    c_over_g = w.count("G") / numC
                                else:
                                    numC = 1
                                    c_over_g = (w.count("G") + 1) / numC
                                if c_over_g > 1:
                                    if re.search(pattern2, w):
                                        data.append(c_over_g)
                                    else:
                                        data.append(0)
                                else:
                                    data.append(0)
                            else:
                                break
                        maxP = np.argmax(data)
                        c_over_g = np.max(data)
                        x1 = x1 + maxP
                        x2 = x2 + maxP
                        scale = x2 - j - 1
                        s = genome[x1 - 150:x1]
                        score = 3
                        ctrl = palindrome_finder(s)
                        if ctrl > 0 or re.search(pattern_pause_site2, s):
                            found_terminators.append((int(x1), int(x2)))  # Store found terminators

        # Add annotations to the sequence record
        for x1, x2 in found_terminators:
            feature = SeqFeature(FeatureLocation(start=int(x1), end=int(x2)), type="RhoTerm", qualifiers={
                "note": [f"Predicted Rho-dependent terminator at {x1}-{x2}"]
            })
            seq_record.features.append(feature)

    # Save the updated GenBank file back to the original file
    SeqIO.write(seq_record, file_genome, "genbank")

    print("Work finished, the original file has been updated.")

except IOError:
    print(f"File {file_genome} not existent in the current directory!")