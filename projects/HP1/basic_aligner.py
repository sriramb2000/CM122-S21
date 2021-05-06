import sys
import argparse
import time
import zipfile
from collections import defaultdict
from collections import Counter


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need
"""
def bwt(text):
    rotated_texts = [text[len(text)-i:] + text[:len(text)-i] for i in range(len(text))]
    sorted_rotated_texts = sorted(rotated_texts)
    output = ""
    for i in range(len(sorted_rotated_texts)):
        output += sorted_rotated_texts[i][-1]
    return output

def better_bwt_matching(bwt_text, pattern, first_occurrence, suffix_array):
    top = 0
    bottom = len(bwt_text) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if symbol in bwt_text[top:bottom+1]:
                top = first_occurrence[symbol] + bwt_text[:top].count(symbol)
                bottom = first_occurrence[symbol] + bwt_text[:bottom+1].count(symbol) - 1
            else:
                return []
        else:
            return [suffix_array[i] for i in range(top, bottom+1)]

def kmer_matching(text, pattern, d, bwt_text, first_occurrence, suffix_array):
    k_size = len(pattern) // (d + 1)
    sol = set()
    for i in range(0, len(pattern), k_size):
        indices = better_bwt_matching(bwt_text, pattern[i:i + k_size], first_occurrence, suffix_array)
        for j in sorted(indices):
            start = j - i
            if start < 0: continue
            mismatches = 0
            for k in range(start, start+len(pattern)):
                if mismatches > d or k >= len(text): 
                    mismatches = d + 1
                    break
                if text[k] != pattern[k-start]:
                    mismatches += 1
            if mismatches <= d:
                sol.add(start)
    return list(sol)

def get_first_occurrence_map(first_col):
    first_occurrence = {}
    for i in range(len(first_col)):
        if first_col[i] not in first_occurrence:
            first_occurrence[first_col[i]] = i
    return first_occurrence

def get_suffix_array(bwt_text):
    suffixes = []
    for i in range(len(bwt_text)):
        suffixes.append((reference[i:], i))
    suffixes.sort()
    return [x[1] for x in suffixes]

def get_votes_by_reference_position(reference, read1, read1_index, read2, read2_index, limit, vote_dict):
    for i in range(read1_index, read1_index + len(read1)):
        if reference[i] != read1[i-read1_index]:
            vote_dict[i][read1[i-read1_index]] += 1
        else:
            vote_dict[i]['original'] += 1
    
    read2_bwd = read2[::-1]
    for i in range(read2_index, read2_index + len(read2)):
        if reference[i] != read2_bwd[i-read2_index]:
            vote_dict[i][read2_bwd[i-read2_index]] += 1
        else:
            vote_dict[i]['original'] += 1

    return vote_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    """
        TODO: Call functions to do the actual read alignment here
        
    """
    reference += '$'
    bwt_text = bwt(reference)
    first_col = sorted(bwt_text)
    first_occurrence = get_first_occurrence_map(first_col)
    suffix_array = get_suffix_array(bwt_text)
    snps = set()
    vote_dict = defaultdict(Counter)
    d = 1
    for [read1, read2] in input_reads:
        # This basically finds start indices of length 50 of a read with <= 1 mismatch
        read1_indices = \
            kmer_matching(reference, read1, d, bwt_text, first_occurrence, suffix_array) \
        # + kmer_matching(reference, read1[::-1], d, bwt_text, first_occurrence, suffix_array)   
        read2_indices = \
            kmer_matching(reference, read2[::-1], d, bwt_text, first_occurrence, suffix_array) \
            # + kmer_matching(reference, read2, d, bwt_text, first_occurrence, suffix_array) 
        for i in read1_indices:
            for j in read2_indices:
                if j - (i + 50) < 111 and j - (i + 50) > 89:
                    get_votes_by_reference_position(reference, read1, i, read2, j, d, vote_dict)

    for pos,votes in vote_dict.items():
        curr_max = (0, None)
        for key,count in votes.items():
            if count > curr_max[0]:
                curr_max = (count, key)
        if curr_max[1] != 'original' and curr_max[1] != None:
            snps.add((reference[pos], curr_max[1], pos))

    # snps = sorted(list(snps), key=lambda x: x[2])
    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)