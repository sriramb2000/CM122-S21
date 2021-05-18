import sys
import argparse
import time
import zipfile
from itertools import zip_longest, islice
from collections import defaultdict
from collections import Counter

#  From https://louisabraham.github.io/notebooks/suffix_arrays.html

def to_int_keys_best(l):
    """
    l: iterable of keys
    returns: a list with integer keys
    """
    seen = set()
    ls = []
    for e in l:
        if not e in seen:
            ls.append(e)
            seen.add(e)
    ls.sort()
    index = {v: i for i, v in enumerate(ls)}
    return [index[v] for v in l]

def suffix_array_best(s):
    """
    suffix array of s
    O(n * log(n)^2)
    """
    n = len(s)
    k = 1
    line = to_int_keys_best(s)
    while max(line) < n - 1:
        line = to_int_keys_best(
            [a * (n + 1) + b + 1
             for (a, b) in
             zip_longest(line, islice(line, k, None),
                         fillvalue=-1)])
        k <<= 1
    return line

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


"""
    TODO: Use this space to implement any additional functions you might need

"""
# # DEPRECATED: use bwt_from_suffix_array in junction with build_suffix_array instead
# def bwt(text):
#     rotated_texts = [text[len(text)-i:] + text[:len(text)-i] for i in range(len(text))]
#     sorted_rotated_texts = sorted(rotated_texts)
#     output = ""
#     for i in range(len(sorted_rotated_texts)):
#         output += sorted_rotated_texts[i][-1]
#     return output

def bwt_from_suffix_array(suffix_array, ref_text):
    bwt = ""
    for i in range(len(ref_text)):
        bwt += ref_text[suffix_array[i]-1] # use BWT <-> SA correspondence
    return bwt

def better_bwt_matching(bwt_text_len, pattern, first_occurrence, suffix_array, symbol_counts):
    top = 0
    bottom = bwt_text_len - 1 # len(bwt_text) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if symbol_counts[top-1][symbol] < symbol_counts[bottom][symbol]: # symbol in bwt_text[top:bottom+1]:
                top = first_occurrence[symbol] + symbol_counts[top-1][symbol] # bwt_text[:top].count(symbol)
                bottom = first_occurrence[symbol] + symbol_counts[bottom][symbol] - 1 # bwt_text[:bottom+1].count(symbol) - 1
            else:
                return []
        else:
            return [suffix_array[i] for i in range(top, bottom+1)]

def kmer_matching(text, pattern, d, bwt_text_len, first_occurrence, suffix_array, symbol_counts):
    k_size = len(pattern) // (d + 1)
    if k_size == 0: return []
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

def get_first_occurrence_map(first_col):
    first_occurrence = {}
    for i in range(len(first_col)):
        if first_col[i] not in first_occurrence:
            first_occurrence[first_col[i]] = i
    return first_occurrence

# """
# More efficient suffix array construction algorithm
# Source: https://www.geeksforgeeks.org/suffix-array-set-2-a-nlognlogn-algorithm/
# """
# def build_suffix_array(ref_text):
#     suffixes = [{'index': i, 'rank': [ord(ref_text[i]) - ord('a'), ord(ref_text[i+1]) - ord('a') if i < len(ref_text)-1 else -1]} for i in range(len(ref_text))]
#     suffixes = sorted(suffixes, key = lambda x: (x['rank'][0], x['rank'][1]))
#     ind = [0] * len(ref_text)
#     k = 4
#     while (k < 2*len(ref_text)):
#         rank, prev_rank = 0, suffixes[0]['rank'][0]
#         suffixes[0]['rank'][0] = rank
#         ind[suffixes[0]['index']] = 0
#         for i in range(1, len(ref_text)):
#             if (suffixes[i]['rank'][0] == prev_rank and suffixes[i]['rank'][1] == suffixes[i - 1]['rank'][1]):
#                 prev_rank = suffixes[i]['rank'][0]
#                 suffixes[i]['rank'][0] = rank  
#             else: 
#                 prev_rank = suffixes[i]['rank'][0]
#                 rank += 1
#                 suffixes[i]['rank'][0] = rank
#             ind[suffixes[i]['index']] = i

#         for i in range(len(ref_text)):
#             nextindex = suffixes[i]['index'] + k // 2
#             suffixes[i]['rank'][1] = suffixes[ind[nextindex]]['rank'][0] if (nextindex < len(ref_text)) else -1

#         suffixes = sorted(suffixes, key = lambda x: (x['rank'][0], x['rank'][1]))
#         k *= 2

#     return [suffixes[i]['index'] for i in range(len(suffixes))]

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

def fitting_alignment_affine_gap_penalty(v, w, scoring_matrix, sigma, epsilon):
    '''
    Returns the score and local alignment substrings for strings v, w with the
    given scoring matrix, gap opening penalty sigma, and gap extension penalty epsilon.
    '''

    # Initialize the matrices.
    S_lower = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    S_middle = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    S_upper = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    backtrack = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]

    # Initialize the maximum score below the lowest possible score.
    max_score = -1

    # Fill in the Score and Backtrack matrices.
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            S_lower[i][j] = max([S_lower[i-1][j] - epsilon, S_middle[i-1][j] - sigma])
            S_upper[i][j] = max([S_upper[i][j-1] - epsilon, S_middle[i][j-1] - sigma])
            middle_scores = [S_lower[i][j], S_middle[i-1][j-1] + scoring_matrix[v[i-1]][w[j-1]], S_upper[i][j]]
            S_middle[i][j] = max(middle_scores)
            backtrack[i][j] = middle_scores.index(S_middle[i][j])

            if S_middle[i][j] > max_score:
                max_score = S_middle[i][j]

    # Initialize the indices to start at the position of the high score.
    j = len(w)
    i = max(enumerate([S_middle[row][j] for row in range(len(w), len(v))]),key=lambda x: x[1])[0] + len(w)

    # Initialize the aligned strings as the input strings up to the position of the high score.
    v_aligned, w_aligned = v[:i], w[:j]

    insert_indel = lambda word, i: word[:i] + '-' + word[i:]

    # Backtrack to start of the local alignment starting at the highest scoring cell.
    # Note: the solution format specifically asks for substrings, so no indel insertion necessary.
    while i*j != 0:
        if backtrack[i][j] == 0:
            i -= 1
            w_aligned = insert_indel(w_aligned, j)
        elif backtrack[i][j] == 1:
            i -= 1
            j -= 1
        elif backtrack[i][j] == 2:
            j -= 1
            v_aligned = insert_indel(v_aligned, i)

    # Cut the strings at the ending point of the backtrack.
    v_aligned = v_aligned[i:]

    return max_score, i, v_aligned, w_aligned

def get_indels_from_final_read_ref(final_ref, final_read, read_i):
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
                # print("ins", curr_ins)
                curr_ins = ["", None]
            if curr_del[1]:
                deletions.add((curr_del[0], curr_del[1]))
                # print("del", curr_del)
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

    """
        TODO: Call functions to do the actual read alignment here

    """
    # Calling SNPS
    print("Calling SNPs")
    print("Initializing support data structures...")
    reference += '$'

    # build hash table
    lookup_table_25_bp = defaultdict(list)
    for i in range(len(reference) - 25):
        lookup_table_25_bp[reference[i:i+25]].append(i)

    # use efficient way to build suffix array and build bwt from suffix array
    # bwt_text, suffix_array = build_suffix_array(reference)
    # bwt_text = bwt_from_suffix_array(suffix_array, reference)
    suffix_array = suffix_array_best(reference)
    bwt_text = bwt_from_suffix_array(suffix_array, reference)

    # use bwt to get count of each allele at each position (from beginning)
    symbol_counts = defaultdict(Counter)
    for i in range(len(bwt_text)):
        symbol_counts[i] = symbol_counts[i-1].copy()
        symbol_counts[i][bwt_text[i]] += 1
    first_col = sorted(bwt_text)
    first_occurrence = get_first_occurrence_map(first_col)
    snps = set()
    insertions = set()
    deletions = set()
    vote_dict = defaultdict(Counter)
    indel_mismatch_scoring_matrix = defaultdict(Counter)
    for i in ['A', 'C', 'T', 'G']:
        for j in ['A', 'C', 'T', 'G']:
            if i == j: indel_mismatch_scoring_matrix[i][j] = 1
            else: indel_mismatch_scoring_matrix[i][j] = -1
    d = 1
    print("Running matching algorithm...")
    for [read1, read2] in input_reads:
        # SNPS
        # This basically finds start indices of length 50 of a read with <= 1 mismatch
        read1_indices = \
            kmer_matching(reference, read1, d, len(bwt_text), first_occurrence, suffix_array, symbol_counts)
        read2_indices = \
            kmer_matching(reference, read2[::-1], d, len(bwt_text), first_occurrence, suffix_array, symbol_counts)

        snp_found_do_not_check_indels = False
        for i in read1_indices:
            for j in read2_indices:
                if j - (i + 50) < 111 and j - (i + 50) > 89:
                    get_votes_by_reference_position(reference, read1, i, read2, j, d, vote_dict)
                    snp_found_do_not_check_indels = True
        
        if snp_found_do_not_check_indels: continue

        # INDELS
        read1_occurrences, read2_occurrences = [], []
        read2_bwd = read2[::-1]
        for i in [0, 25]:
            read1_occurrences.extend([x - i for x in lookup_table_25_bp[read1[i:i+25]]])
            read2_occurrences.extend([x - i for x in lookup_table_25_bp[read2_bwd[i:i+25]]])

        for i in read1_occurrences:
            for j in read2_occurrences:
                if j - i <= 300 and j - i >= 50 and i > 100 and j > 100:
                    read1_score, read1_i, final_ref1, final_read1 = \
                        fitting_alignment_affine_gap_penalty(\
                            reference[i-100:i+100], read1, indel_mismatch_scoring_matrix, 11, 2)
                        # (v, w, scoring_matrix, sigma, epsilon):
                    read2_score, read2_i, final_ref2, final_read2 = \
                        fitting_alignment_affine_gap_penalty(\
                            reference[j-100:j+100], read2_bwd, indel_mismatch_scoring_matrix, 11, 2)

                    # if read1_score > 0: print(read1_score, read1_i, final_ref1, final_read1)
                    # if read2_score > 0: print(read2_score, read2_i, final_ref2, final_read2)

                    if read1_score + read2_score > 0:
                        ins, dels = get_indels_from_final_read_ref(final_ref1, final_read1, i - 100 + read1_i)
                        insertions.update(ins)
                        deletions.update(dels)

                        ins, dels = get_indels_from_final_read_ref(final_ref2, final_read2, j - 100 + read2_i)
                        insertions.update(ins)
                        deletions.update(dels)

    # Majority vote for SNPS locations
    for pos,votes in vote_dict.items():
        curr_max = (0, None)
        for key,count in votes.items():
            if count > curr_max[0]:
                curr_max = (count, key)
        if curr_max[1] != 'original' and curr_max[1] != None:
            snps.add((reference[pos], curr_max[1], pos))

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
