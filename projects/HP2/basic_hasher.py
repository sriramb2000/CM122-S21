import sys
import argparse
import time
import zipfile
from itertools import zip_longest, islice
from collections import defaultdict
from collections import Counter
from pydivsufsort import divsufsort
import numpy as np

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads

    HINT: This might not work well if the number of reads is too large to handle in memory
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

def auxiliary_from_suffix_array(suffix_array, ref_text):
    bwt = ""
    for i in range(len(ref_text)):
        bwt += ref_text[suffix_array[i]-1] 

    counts = defaultdict(Counter)
    for i in range(len(bwt)):
        counts[i] = counts[i-1].copy()
        counts[i][bwt[i]] += 1
    first_col = sorted(bwt)
    first_occurrences = {}
    for i in range(len(first_col)):
        if first_col[i] not in first_occurrences:
            first_occurrences[first_col[i]] = i

    return bwt, first_occurrences, counts

def better_bwt_matching(bwt_text_len, pattern, first_occurrence, suffix_array, symbol_counts):
    top = 0
    bottom = bwt_text_len - 1 
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if symbol_counts[top-1][symbol] < symbol_counts[bottom][symbol]: 
                top = first_occurrence[symbol] + symbol_counts[top-1][symbol]
                bottom = first_occurrence[symbol] + symbol_counts[bottom][symbol] - 1 
            else:
                return []
        else:
            return [suffix_array[i] for i in range(top, bottom+1)]

def kmer_matching(text, pattern, d, bwt_text_len, first_occurrence, suffix_array, symbol_counts):
    k_size = len(pattern) // (d + 1)
    sol = set()
    for i in range(0, len(pattern), k_size):
        indices = better_bwt_matching(bwt_text_len, pattern[i:i + k_size], first_occurrence, suffix_array, symbol_counts)
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

def accumulate_ref_pos_votes(reference, read1, read1_index, read2, read2_index, ref_pos_votes):
    for i in range(read1_index, read1_index + len(read1)):
        if reference[i] != read1[i-read1_index]:
            ref_pos_votes[i][read1[i-read1_index]] += 1
        else:
            ref_pos_votes[i]['match'] += 1
    
    for i in range(read2_index, read2_index + len(read2)):
        if reference[i] != read2[i-read2_index]:
            ref_pos_votes[i][read2[i-read2_index]] += 1
        else:
            ref_pos_votes[i]['match'] += 1

    return ref_pos_votes

# GRAPH LEVELS
BOTTOM = 0
MIDDLE = 1
TOP = 2

def insert_del(word, i):
    return word[:i] + '-' + word[i:]

def fit_align_with_affine_gap(v, w, sigma, epsilon):
    sB = [[0 for _ in range(len(w)+1)] for _ in range(len(v)+1)]
    sM = [[0 for _ in range(len(w)+1)] for _ in range(len(v)+1)]
    sT = [[0 for _ in range(len(w)+1)] for _ in range(len(v)+1)]
    backtrack = [[0 for _ in range(len(w)+1)] for _ in range(len(v)+1)]

    max_score = -1

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            sB[i][j] = max([sB[i-1][j] - epsilon, sM[i-1][j] - sigma])
            sT[i][j] = max([sT[i][j-1] - epsilon, sM[i][j-1] - sigma])
            match_val = 1 if v[i-1] == w[j-1] else -1
            middle_scores = [sB[i][j], sM[i-1][j-1] + match_val, sT[i][j]]
            sM[i][j] = max(middle_scores)
            backtrack[i][j] = middle_scores.index(sM[i][j])

            if sM[i][j] > max_score:
                max_score = sM[i][j]

    j = len(w)
    i = max(enumerate([sM[row][j] for row in range(len(w), len(v))]),key=lambda x: x[1])[0] + len(w)

    v_aligned, w_aligned = v[:i], w[:j]

    while i != 0 and j != 0:
        if backtrack[i][j] == BOTTOM:
            i -= 1
            w_aligned = insert_del(w_aligned, j)
        elif backtrack[i][j] == MIDDLE:
            i -= 1
            j -= 1
            # Don't need to worry about mismatches, only indels
        elif backtrack[i][j] == TOP:
            j -= 1
            v_aligned = insert_del(v_aligned, i)

    # Cut the strings at the ending point of the backtrack.
    v_aligned = v_aligned[i:]

    return max_score, i, v_aligned, w_aligned

def indels_from_align_ref(final_ref, final_read, read_i):
    insertions, deletions = set(), set()
    curr_ins, curr_del = ["", None], ["", None]
    for i in range(len(final_ref)):
        if final_ref[i] == '-':
            if curr_ins[1]:
                curr_ins[0] += final_read[i]
            else:
                curr_ins = [final_read[i], read_i + i]
        elif final_read[i] == '-':
            if curr_del[1]:
                curr_del[0] += final_ref[i]
            else:
                curr_del = [final_ref[i], read_i + i]
        else:
            if curr_ins[1]:
                insertions.add((curr_ins[0], curr_ins[1]))
                curr_ins = ["", None]
            if curr_del[1]:
                deletions.add((curr_del[0], curr_del[1]))
                curr_del = ["", None]
    return insertions, deletions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    # Calling SNPS
    print("Calling SNPs")
    reference += '$'

        # use efficient way to build suffix array and build bwt from suffix array
    suffix_array = divsufsort(np.unique(np.array(list(reference)), return_inverse=True)[1])
    print("Initialized suffix array")
    bwt_text, first_occurrence, symbol_counts = auxiliary_from_suffix_array(suffix_array, reference)
    # build hash table
    referece_dic_25 = defaultdict(list)
    for i in range(len(reference) - 25):
        referece_dic_25[reference[i:i+25]].append(i)

    snps = set()
    insertions = set()
    deletions = set()
    ref_pos_votes = defaultdict(Counter)
    d = 1
    sigma = 11
    epsilon = 2
    print("Running algorithm...")
    ccc = 1
    for [read1, read2] in input_reads:
        print("Processing Read {}".format(ccc))
        ccc += 1
        # SNPS
        maxMatchCount = -1
        max_read1_indices, max_read2_indices, max_read1, max_read2 = None, None, None, None
        for flip1 in [1, -1]:
            for flip2 in [1, -1]:
                if flip1 == -1 and flip2 == -1:
                    continue
                read1_, read2_ = read1[::flip1], read2[::flip2]
                read1_indices = kmer_matching(reference, read1_, d, len(bwt_text), first_occurrence, suffix_array, symbol_counts)
                read2_indices = kmer_matching(reference, read2_, d, len(bwt_text), first_occurrence, suffix_array, symbol_counts)
                matchCount = 0
                for i in read1_indices:
                    for j in read2_indices:
                        if j - (i + 50) < 111 and j - (i + 50) > 89:
                            matchCount += 1
                if matchCount > maxMatchCount:
                    max_read1_indices, max_read2_indices, max_read1, max_read2 = read1_indices, read2_indices, read1_, read2_
                    maxMatchCount = matchCount

        snpFound = False
        for i in max_read1_indices:
            for j in max_read2_indices:
                if j - (i + 50) < 111 and j - (i + 50) > 89:
                    accumulate_ref_pos_votes(reference, max_read1, i, max_read2, j, ref_pos_votes)
                    snpFound = True
        
        if snpFound: continue

        # INDELS
        maxAlignCount = -1
        max_read1_occurrences, max_read2_occurrences, max_read1, max_read2 = None, None, None, None

        for flip1 in [1, -1]:
            for flip2 in [1, -1]:
                read1_ = read1[::flip1]
                read2_ = read2[::flip2]
                read1_occurrences, read2_occurrences = [], []

                for i in [0, 25]:
                    read1_occurrences.extend([x - i for x in referece_dic_25[read1_[i:i+25]]])
                    read2_occurrences.extend([x - i for x in referece_dic_25[read2_[i:i+25]]])
                alignCount = 0
                for i in read1_occurrences:
                    for j in read2_occurrences:
                        if j - i <= 300 and j - i >= 50 and i > 100 and j > 100:
                            alignCount += 1

                if alignCount > maxAlignCount:
                    max_read1_occurrences, max_read2_occurrences, max_read1, max_read2 = read1_occurrences, read2_occurrences, read1_, read2_
                    maxAlignCount = alignCount

        for i in max_read1_occurrences:
            for j in max_read2_occurrences:
                if j - i <= 300 and j - i >= 50 and i > 100 and j > 100:
                    read1_score, read1_i, read1_ref, final_read1 = \
                        fit_align_with_affine_gap(\
                            reference[i-100:i+100], max_read1, sigma, epsilon)
                    read2_score, read2_i, read2_ref, final_read2 = \
                        fit_align_with_affine_gap(\
                            reference[j-100:j+100], max_read2, sigma, epsilon)

                    if read1_score + read2_score > 0:
                        ins1, dels1 = indels_from_align_ref(read1_ref, final_read1, i - 100 + read1_i)
                        ins2, dels2 = indels_from_align_ref(read2_ref, final_read2, j - 100 + read2_i)

                        insertions.update(ins1)
                        deletions.update(dels1)
                        insertions.update(ins2)
                        deletions.update(dels2)

    # Majority vote for SNPS locations
    for pos,votes in ref_pos_votes.items():
        curr_max = max(votes.items(), key=lambda x: x[1])
        if curr_max[0] != 'match' and curr_max[0] != None:
            snps.add((reference[pos], curr_max[0], pos))

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
