from typing import Optional, List
from collections import defaultdict
import sys

def constructGraph(kmers: List[str]):
    suffixToPrefix = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffixToPrefix[prefix].append(kmer)
    
    for kmer in kmers:
        suffix = kmer[1:]
        if suffix in suffixToPrefix:
            print("{} -> {}".format(kmer, ",".join(suffixToPrefix[suffix])))

    return

def part1(content: List[str]) -> None:
    print(content[0] + "".join([s[-1] for s in content[1:]]))

def part2(content: List[str]) -> None:
    constructGraph(content)

def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: python sol.py 1|2 <filename>")
        return
    
    filename = sys.argv[2]
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content] 

    if sys.argv[1] == "1":
        part1(content)
    elif sys.argv[1] == "2":
        part2(content)

if __name__ == "__main__":
    main()