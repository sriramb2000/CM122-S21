import sys
import argparse
import time
import zipfile
from collections import defaultdict
from collections import Counter

def cmp(a, b):
    return (a>b)-(a<b)

# Linear time BWT/Suffix Array Construction from https://github.com/michaelting/BWT/blob/master/bwt.py

def process_FASTA(infile):
    """
    Generator that reads a FASTA file and outputs line by line
    """
    name, seq = None, []
    for line in infile:
        if line.startswith(">"):            # sequence header found
            if name:                        # next sequence in file found
                yield (name, ''.join(seq))  # output the sequence
            name, seq = line.strip(), []    # create a new sequence item
        else:
            seq.append(line.strip())        # continuing portion of sequence appended
    if name: 
        yield (name, ''.join(seq))	

def prep_DC3(charseq):
    """
    Prepare for DC3 to run on numerical lists since our input is a sequence
    of characters
    Input:
        charseq - our sequence of interest from our FASTA file
    Output:
        numlst  - the numerical representation of our characters
    """
    # stored as [0,1,]
    numlst = []
    # convert chars to ascii numerical values
    for i in range(len(charseq)):
        char = charseq[i]
        num = int(ord(char))    # ensure value stored as an int
        numlst.append(num)
    return numlst

def DC3(numlst):
    """
    DC3/Skew Algorithm
    Computes the suffix array S for the string seq
    Input:
        numlst - the sequence of interest represented by a numerical list
    Output:
        suffixarray  - the suffix array of the input sequence numlst
    """
    B0, B1, B2, C = DC3_sample(numlst) # {C} = {B1}U{B2}
    
    # extract suffixes and indices of suffix start positions
    samsuftups = DC3_num_to_suffix(numlst, C)
    nonsamsuftups = DC3_num_to_suffix(numlst, B0)
    
    # sort sample suffixes
    sortedsamples, indexrankmap = DC3_samsort(samsuftups)
    
    # pad rank(S_n+1) and rank(S_n+2) with 0's
    n = len(numlst)
    indexrankmap[n]     = 0
    indexrankmap[n+1]   = 0
    
    # pad the ends
    numlst.extend([0,0])
    
    # sort the nonsample suffixes using the indexrankmap
    sortednonsamples = DC3_nonsamsort(nonsamsuftups, numlst, indexrankmap)
    
    # merge the sets of sorted suffixes into one big sorted list
    allsortedsuffixes, suffixarray = DC3_mergesamples(numlst, sortednonsamples, sortedsamples, indexrankmap)
    
    return suffixarray

    
def DC3_sample(numlst):
    """
    Sort the indices and return the sample indices for B0, B1, B2, and B1UB2
    Input:
        numlst - the sequence of interest in numerical list representation
    Output:
        B0  - the set of indices i for which i = 0mod3
        B1  - the set of indices i for which i = 1mod3
        B2  - the set of indices i for which i = 2mod3
        C   - the set of indices that is the union of {B1} and {B2}
    """
    B0, B1, B2 = [], [], []
    # sort our indices modulo 3 into our samples
    for i in range(len(numlst)):
        if i % 3 == 1:
            B1.append(i)
        elif i % 3 == 2:
            B2.append(i)
        elif i % 3 == 0:
            B0.append(i)
        else:
            raise Exception("Invalid index found: %d" % i)
    # join our difference cover sample positions
    C = B1[:]       # copy list by value, not by reference
    C.extend(B2)    # join B1 and B2 by value
    return B0, B1, B2, C
    
def DC3_num_to_suffix(numR, samind):
    """
    Convert numerical sequence representation to suffixes for recursions
    Input:
        numR    - our new string R represented numerically as a list of ints
        samind  - sample indices corresponding to R1'UR2'
    Output:
        suffixes     - a list of tuples, with each tuple appearing:
                    ([65,2,1],1), ([34,0,0],4)
                    So the list looks like
                    [([65,2,1],1), ([34,0,0],4), ...]
        Tuples contain (suffix, index) where index corresponds to the position
        of the beginning of the suffix in the original number list.
    """
    # suffixes looks like [[1,3,2],[4,5,0],...]
    suffixes = []
    for index in samind:
        sufnumlst = numR[index:index+3]    # extract suffix of length 2 starting at index
        while len(sufnumlst) < 3:   # pad short suffixes with 0's
            sufnumlst.append(0)
        # each tuple stores the suffix and the index it came from in the original numlist
        sufnumtuple = (sufnumlst, index)
        suffixes.append(sufnumtuple)
    return suffixes
    
def DC3_samsort(suftuplst):
    """
    Sort the samples from C = B1UB2
    Input:  
        suftuplst   - a list of tuples of (suffix, origindex)
    Output:
        ordered     - list of ordered suffixes
        indexrankmap    - a mapping of original string indices i to the suffix rank, rank(S_i)
    """

    ITER = 0
    
    # lexicographically sorted suffix tuples (suffix, index)
    ordered = DC3_radixsort(suftuplst, ITER)
    
    # determine whether we need to recurse on R' (flag)
    # get a mapping of original suffix index to its rank {index:rank}
    # checked contains exactly one copy of a suffix (elim repeats)
    flag, indexrankmap, checked = DC3_check_radix(ordered)

    # create new string R' using mappings from radix sort on original suffixes
    newR = []
    for suftup in suftuplst:
        suf = suftup[0]     # suffix
        sufind = suftup[1]  # index
        newR.append(indexrankmap[sufind])
    
    # need to recurse on newR since a duplicate suffix is found in the suffixes
    while flag:
        # DC3 should return a list of sorted suffixes
        suffixarray = DC3(newR)  #suffixarray = rank:index
        
        invsuffixarray = dict((v,k) for k,v in suffixarray.items()) # index: rank
        # sort by ranks
        indices = sorted(invsuffixarray, cmp=lambda x,y:cmp(x,y))
        # get ordering of suffixes using suffix array information
        # relabel our suffixes using our ranks, which are now unique from recursion
        arraysortedsamples = []
        for index in indices:
            rank = invsuffixarray[index]
            oldtup = suftuplst[index]
            newtup = ([rank], oldtup[1])  # rank becomes the new suffix name
            arraysortedsamples.append(newtup)
            
        ordered = DC3_radixsort(arraysortedsamples, ITER)
        flag, indexrankmap, checked = DC3_check_radix(ordered)
    
    return ordered, indexrankmap
    
def DC3_radixsort(suftuplst, iternum):
    """
    Sort the suffixes in linear time using radix sort
    Input:
        samlst  - a list of samples (triplets of indices)
        iternum - the iteration number, valid from 0 to 2, base case at 3
    Output:
        ordered - a list of tuples of numerically ordered samples [ ([1,2,3],1), ([2,4,5],5),...]
    """    
    # stop at iternum = 3 since we are in mod 3 and suffixes are max length 3
    if iternum > 2:
        return suftuplst
    
    # buckets have form {1:[1,2,3],[1,7,2],...]
    #                    3:[3,2,5],[3,18,20],...}
    buckets = {}
    # process all sample suffixes using radix sort
    for suftuple in suftuplst:
        # sort suffixes into buckets by character at index iternum {0, 1, 2}
        sufnumlst = suftuple[0]
        label = sufnumlst[iternum]
        # create bucket if it does not yet exist
        if label not in buckets:
            buckets[label] = []
        # buckets will look like [ ([1,2,3],1), ([1,2,4],3), ...]]
        buckets[label].append(suftuple)
    # sort bucket keys numerically increasing, 1,2,3,...
    keysortedbuckets = sorted(buckets, cmp=lambda x,y:cmp(x,y)) # sort by numerical value
    # order the suffixes
    ordered = []
    for key in keysortedbuckets:
        # get list of tuples in the bucket
        samples = buckets[key]
        # recurse on buckets with multiple tuples
        if len(buckets[key]) > 1:   # more than one tuple in a bucket
            samples = DC3_radixsort(samples, iternum+1) # use next indexed num in suffix
        ordered.extend(samples)
    return ordered

def DC3_check_radix(sortedtuplst):
    """
    Determines whether or not we need to recurse with R'
    Input:
        sortedtuplst    - result from running DC3_radixsort
    Output:
        flag    - True if we need to recurse with DC3 on R' (duplicates found)
                - False if all characters are different, no need to recurse
sortedsamranks  - a list of ranks mapping on the same indices to sortedsamlst
        indexrankmap    - a dictionary mapping original suffix index to its rank
        checked - a dictionary of samnumlist:rank, {(1,2,3):1,(2,3,4):2,...}
    """
    checked = {}
    indexrankmap = {}
    rank = 1
    flag = False
    # rank samples from the radix-sorted list of tuples
    for samnumtuple in sortedtuplst:
        samnumlst = samnumtuple[0]  # extract numlst from tuple
        sufindex = samnumtuple[1]   # index of suffix from original string
        samnumlst = tuple(samnumlst)    # need immutable keys       
        # sample is unique
        if samnumlst not in checked:
            checked[samnumlst] = rank
            #sortedsamranks.append(rank)
            indexrankmap[sufindex] = rank
            rank += 1
        # sample is a repeat
        else:
            flag = True
            #sortedsamranks.append(checked[samnumlst])
            indexrankmap[sufindex] = checked[samnumlst] # rank tie, use from table
    return flag, indexrankmap, checked

def DC3_nonsamsort(nonsortedsamtups, numlst, indexrankmap):
    """
    Sort the non-sample suffixes using information from the sorted sample suffixes
    Input:
        nonsortedsamtups    - a list of tuples (suffix, origindex) from B0
        numlst      - the original list of numbers representing our sequence
        indexrankmap    - a mapping from indices i in our original numlst to suffix rank, rank(S_i)
    Output:
        A non-sample indexrankmap dictionary
        OR a sorted list of tuples of non-sample suffixes and their indices
    """
    
    nonsampairtuplst = []
    for samtup in nonsortedsamtups:
        samlst = samtup[0]
        samind = samtup[1]
        # create the pair (t_i, rank(S_i+1))
        t_i = samlst[0] # the first num in our suffix        
        trank = indexrankmap[samind+1]  # rank(S_i+1)
        nonsampair = (t_i, trank)
        nonsampairtup = (nonsampair, samind)    # keep track of original index
        nonsampairtuplst.append(nonsampairtup)
    # radix sort the non-sample suffix tuples
    sortednonsampairtuplst = DC3_nonsamradixsort(nonsampairtuplst, 0)    
    return sortednonsampairtuplst
    
def DC3_nonsamradixsort(nonsampairtuplst, iternum):
    """
    Sort the non-sample suffix tuples (t_i, trank) where trank = rank(S_i+1)
    Input:
        nonsampairtuplst   - a list of non-sample suffix tuples [((t_i, rank(S_i+1)), origindex), ...]
        iternum     - the iteration number, which is in {0, 1} since we only have 2 values for radix sort
    Output:
        ordered - the radix-sorted list of non-sample suffix tuples with their original indices
    """

    # stop at iternum = 2 since we only have 2 values to compare
    if iternum > 1:
        return nonsampairtuplst

    # buckets have form {1:[ (1,1), (1,4), ...]
    #                    3:[ (3,2), (3,20), ...}
    buckets = {}
    # process all sample suffixes using radix sort
    for nonsampairtup in nonsampairtuplst:
        nonsampair = nonsampairtup[0]   # extract the (t_i, rank(S_i+1)) tuple
        # need to sort by trank since chars t are equivalent
        if iternum > 0:
            label = nonsampair[1]   # sorting by trank
            if label not in buckets:
                buckets[label] = []
            buckets[label].append(nonsampairtup)
        # first sort by chars t
        else:
            # sort suffixes into buckets
            label = nonsampair[0]   # sorting into buckets by t_i
            # create bucket if it does not yet exist
            if label not in buckets:
                buckets[label] = []
            # buckets will look like [ (1,1), (2,3), ...]]
            buckets[label].append(nonsampairtup)
    # sort bucket keys numerically increasing, 1,2,3,...
    keysortedbuckets = sorted(buckets, cmp=lambda x,y:cmp(x,y)) # sort by numerical value
    # order the suffixes
    ordered = []
    for key in keysortedbuckets:
        # get list of tuples in the bucket
        samples = buckets[key]
        # recurse on buckets with multiple tuples
        if len(buckets[key]) > 1:   # more than one tuple in a bucket
            samples = DC3_nonsamradixsort(samples, iternum+1) # use next indexed num in suffix
        ordered.extend(samples)
    return ordered

def DC3_mergesamples(numlst, sortedB0, sortedC, indexrankmap):
    """
    Merge the sample sets C and B0
    
    B0 form:    [((0, 0), 12), ((1, 5), 6), ((1, 7), 9), ((2, 2), 3), ((3, 1), 0)]
    C=B1UB2 form: [([1, 2, 3], 1), ([2, 3, 4], 2), ([4, 5, 6], 4), ([5, 6, 0], 5)]
    
    S_i corresponds to C=B1UB2
    S_j corresponds to B0
    
    """
    TUPIND = 0
    SUFFIX = 0
    INDEX = 1    
    
    allsorted = []
    suffixarray = {}  # rank:index
    # loop through every suffix
    k = 0
    while k < len(numlst):
        # no more pairs in either list
        if (not sortedC) and (not sortedB0):
            break
        # no more pairs in C=B1UB2
        if not sortedC:
            suffixarray[k] = sortedB0[TUPIND][INDEX]
            allsorted.append(sortedB0.pop(0)) # remove first item from S_j
            k += 1
            continue
        if not sortedB0:
            suffixarray[k] = sortedC[TUPIND][INDEX]
            allsorted.append(sortedC.pop(0))    # remove first item from S_i
            k += 1
            continue
            
        # extract information from sorted lists
        S_i = sortedC[TUPIND]   # ((0, 0), 12)      (suffix, index)
        S_i_suffix = S_i[SUFFIX]
        S_i_index = S_i[INDEX]
        S_j = sortedB0[TUPIND]    # ([1, 2, 3], 1)    (suffix, index)
        S_j_suffix = S_j[SUFFIX]
        S_j_index = S_j[INDEX]
        
        if S_i_index % 3 == 1:
            # handle B1
            # iff ((t_i, rank(S_i+1)) <= (t_j, rank(S_j+1))), then S_i <= S_j
            # compare:
            ipair = (numlst[S_i_index], indexrankmap[S_i_index+1])
            jpair = (numlst[S_j_index], indexrankmap[S_j_index+1])
        elif S_i_index % 3 == 2:
            # handle B2
            # iff (t_i, t_i+1, rank(S_i+2)) <= (t_j, t_j+1, rank(S_j+2)), then S_i <= S_j
            ipair = (numlst[S_i_index], numlst[S_i_index+1], indexrankmap[S_i_index+2])
            jpair = (numlst[S_j_index], numlst[S_j_index+1], indexrankmap[S_j_index+2])
        else:
            raise Exception("Invalid suffix index found in C=B1UB2")

        comparison = compare(ipair, jpair)
        # ipair > jpair
        if comparison > 0:
            suffixarray[k] = S_j_index
            allsorted.append(sortedB0.pop(0)) # remove first item from S_j
        # ipair < jpair
        elif comparison < 0:
            suffixarray[k] = S_i_index
            allsorted.append(sortedC.pop(0))  # remove first item from S_i
        elif comparison == 0:
            raise Exception("Pairs found to be equal!")
        
        k += 1
    
    return allsorted, suffixarray
            
def compare(x, y):
    """
    Returns:
    1 for x > y
    0 for x = y
    -1 for x < y
    """
    # check lengths are the same
    if len(x) != len(y):
        raise Exception("Lengths of merging pairs to compare not equal!")
        
    result = 0
    # check ordering
    for i in range(len(x)):
        if x[i] < y[i]:
            result = -1
            break
        elif x[i] > y[i]:
            result = 1
            break
        
    return result    

def encode_BWT(seq):
    """
    Constructs the Burrows-Wheeler Transform of the input sequence
    Input:
        seq - the input sequence (assumes $ at end of seq)
    Output:
        BWT - the Burrows-Wheeler Transform of seq (which is L in pi_sorted)
    """
    
    suffixarray = DC3(prep_DC3(seq))
    L = []
    for rank in range(len(seq)):
        L.append(seq[suffixarray[rank]-1])
    BWT = ''.join(L)
    return BWT, suffixarray

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

# def bwt_from_suffix_array(suffix_array, ref_text):
#     bwt = ""
#     for i in range(len(ref_text)):
#         bwt += ref_text[suffix_array[i]-1] # use BWT <-> SA correspondence
#     return bwt

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
    bwt_text, suffix_array = encode_BWT(reference)

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
