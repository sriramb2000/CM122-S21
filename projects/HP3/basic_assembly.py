from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))


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

def generate_contigs_from_reads(graph):
    degrees = graph_degrees(graph)
    contigs = []

    for v in graph.keys():
        if degrees[v] == [1, 1]:
            continue
        for u in graph[v]:
            contig = v
            w = u
            while True:
                contig += w[-1]
                w_degree = degrees[w]
                if w_degree == [1, 1]:
                    w = graph[w][0]
                else:
                    break
            contigs.append(contig)
    return sorted(contigs)

def debrujin_graph_from_reads(reads, k):
    de_bruijn_counter = defaultdict(Counter)
    for read in reads:
        # Cut the read into k-mers
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])

    # Eemoves the nodes from the DeBruijn Graph that we have not seen enough.
    de_bruijn_graph = defaultdict(list)
    for key, val in de_bruijn_counter.items():
        if not val:
            continue
        for end, count in val.items():
            if count > 3: # could be 5
                de_bruijn_graph[key].append(end)
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}
    return de_bruijn_graph

def graph_degrees(graph):
    degrees = {}
    for i in graph.keys():
        neighbors = graph[i]
        out_degree = len(neighbors)

        if i in degrees:
            degrees[i][1] = out_degree
        else:
            degrees[i] = [0, out_degree]

        for j in neighbors:
            if j in degrees:
                degrees[j][0] += 1
            else:
                degrees[j] = [1, 0]

    return degrees


"""
    TODO: Use this space to implement any additional functions you might need

"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    """
            TODO: Call functions to do the actual assembly here

    """
    all_reads = [item for sublist in input_reads for item in sublist]
    contigs = generate_contigs_from_reads(debrujin_graph_from_reads(all_reads, 30))

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
