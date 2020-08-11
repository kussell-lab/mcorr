#!/usr/bin/env python3
"""
This script PRECISELY chooses some {number} strains from the gene-by-gene alignments.
Created by Asher, based on script by Mingzhi Lin.
"""
import sys
import random
def parse_description_line(line):
    """
    Parse description line of a FASTA record.
    And return the strain and alignment names.
    """
    line = line.rstrip()[1:] # ignore '>'
    words = line.split(' ')
    alignment_name = words[len(words)-1]
    strain_name = words[0].split(':')[0]
    return (strain_name, alignment_name)

def read_fasta(filename):
    """
    Read a FASTA or XMFA file.
    And return a list of (strain, alignment, sequence).
    """
    sequence_list = []
    strain_name = None
    alignment_name = None
    sequence = ""
    with open(filename, 'rU') as reader:
        for line in reader:
            line = line.strip()
            if line.startswith('>'): # description line
                if strain_name is not None:
                    sequence_list.append((strain_name, alignment_name, sequence))
                strain_name, alignment_name = parse_description_line(line)
                sequence = ""
            elif line.startswith('=') or len(line) == 0:
                continue
            else:
                sequence = sequence + line
    if strain_name is not None and sequence != "":
        sequence_list.append((strain_name, alignment_name, sequence))
    return sequence_list

def filter_gapped_sequences(sequences):
    """
    Filter sequences with too many gaps
    """
    results = []
    for (strain, sequence) in sequences:
        num_gap = 0
        for c in sequence:
            if c == '-':
                num_gap = num_gap + 1
        if float(num_gap) < float(len(sequence)) * 0.02:
                results.append((strain, sequence))
    return results

def main():
    """chooses some number of strains from the gene-by-gene alignment file"""
    align_file = sys.argv[1]
    out_file = sys.argv[2]
    num = int(sys.argv[3])
    sequence_list = read_fasta(align_file)
    strain_set = set()
    for strain, _, _ in sequence_list:
        strain_set.add(strain)
    strains = list(strain_set)
    choices = set(strains[0:(num)])
    #choices = set(strains[0])
    alignments = {}
    for (strain, alignment_name, sequence) in sequence_list:
        #if strain in choices:
        if strain == '661_HUP_B14':
            sequences = alignments.get(alignment_name, [])
            sequences.append((strain, sequence))
            alignments[alignment_name] = sequences
    with open(out_file, 'w') as writer:
        for alignment_name, sequences in alignments.items():
            sequences = filter_gapped_sequences(sequences)
            # if len(sequences) < 2:
            #     continue
            for (strain, sequence) in sequences:
                writer.write(">%s %s\n" % (alignment_name, strain))
                writer.write(sequence + "\n")
            writer.write("=\n")

if __name__ == "__main__":
    main()
