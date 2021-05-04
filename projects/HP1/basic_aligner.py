import sys
import argparse
import time
import zipfile
from typing import List
from collections import defaultdict

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

def rotate(strg, n):
    return strg[n:] + strg[:n]

def getBwt(text: str):
    rotated = [(i, rotate(text, i)) for i in range(len(text))]
    rotated.sort(key = lambda x: x[1])
    bwt = [p[1][-1] for p in rotated]
    suffix = [p[0] for p in rotated]

    firstCol = sorted(text)
    bands = {}
    for s in set(text):
        bands[s] = firstCol.index(s)
    return bwt, suffix

def approximateMatching(text: str, patterns: List[str], d: int):
    bwt, suffix = getBwt(text)

    def getCountSymbol(s: str):
        d = {}
        res = [{}]
        for c in s:
            d[c] = d.get(c, 0) + 1
            res.append(d.copy())
        return res
    count_symbol = getCountSymbol(bwt)

    def getFirstOccurrence(s: str):
        ordered = sorted(s)
        d = {}
        for c in set(s):
            d[c] = ordered.index(c)
        return d
    first_occurrence = getFirstOccurrence(bwt)

    def genSeeds(p: str):
        l = len(p)
        assert l>d
        minsize = l//(d+1)
        cut = list(range(0,l-minsize+1,minsize))
        cut.append(l)
        seeds = [(p[cut[i-1]:cut[i]],cut[i-1]) for i in range(1,len(cut))]
        return seeds
    
    def findSeed(seed: str):
        top = 0
        bottom = len(text) - 1
        while top <= bottom:
            if seed:
                symbol = seed[-1]
                seed = seed[:-1]
                if count_symbol[bottom+1].get(symbol, 0) > count_symbol[top].get(symbol, 0):
                    top = first_occurrence[symbol] + count_symbol[top].get(symbol,0)
                    bottom = first_occurrence[symbol] + count_symbol[bottom+1].get(symbol, 0) - 1
                else:
                    break
            else:
                return [suffix[i] for i in range(top, bottom+1)]
        return []
    
    def kMismatch(offset: int, p: str, k: int):
        snps = []
        for i, c in enumerate(p):
            if (c != text[offset+i]):
                k -= 1
                if k < 0:
                    return (False, snps)
                # original, actual, index
                snps.append((text[offset+i], c, offset+i))
        return (True, snps)
    
    def approxPatternPositions(p: str):
        positions = set()
        seedOffsetPairs = genSeeds(p)
        snpz = []
        for (seed, offset) in seedOffsetPairs:
            candidate_poz = findSeed(seed)
            for pos in candidate_poz:
                pattern_pos = pos - offset
                if pattern_pos < 0 or (pattern_pos + len(p) > len(text)):
                    continue
                kMatch, snps = kMismatch(pattern_pos, p, d)
                if kMatch:
                    positions.add(pattern_pos)
                    snpz += snps
        return (list(positions), snpz)
    
    res = []
    all_snps = []
    for pattern in patterns:
        r, s = approxPatternPositions(pattern)
        res += r
        all_snps += s
    return res, all_snps

def flatten_read_list(reads):
    temp = []
    for r in reads:
        temp += r
    return temp

def convert_to_freq_dic(all_snps):
    freq_dic = defaultdict(lambda: defaultdict(int))
    for snp in all_snps:
        original, actual, index = snp
        freq_dic[index][actual] += 1
    
    return freq_dic


def filter_read_errs(all_snps):
    filtered_dic = {}
    for index in all_snps:
        if len(all_snps[index]) == 1:
            only_val_list = list(all_snps[index].keys())
            if all_snps[index][only_val_list[0]] == 1:
                continue
        filtered_dic[index] = all_snps[index].copy()
    
    return filtered_dic

def majority_snps(ref_string, all_snps):
    final_snps = []
    for index in all_snps:
        if len(all_snps[index]) == 1:
            # Only one type of variation at this index, include it
            for actual in all_snps[index]:
                final_snps.append((ref_string[index], actual, index))
        else:
            # Find the majority if possible and include it
            max_val = -1
            max_val_actual = None
            total_count = 0
            for actual in all_snps[index]:
                actual_val = all_snps[index][actual]
                total_count += actual_val
                if actual_val > max_val:
                    max_val = actual_val
                    max_val_actual = actual
            # If there is a majority, include it
            if (max_val - (total_count // 2)) > 0:
                final_snps.append((ref_string[index], max_val_actual, index))
    
    return final_snps


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
    # Get all snps from reads (currently ignores paired ordering)
    _, temp_snps = approximateMatching(reference, flatten_read_list(input_reads), 1) 
    # Filter out apparent read errors
    read_errs_filtered = filter_read_errs(convert_to_freq_dic(temp_snps))

    # Filter out any non-majority snps
    snps = majority_snps(reference, read_errs_filtered)

    # snps = [['A', 'G', 3425]]

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
