import sys
from typing import List
from collections import defaultdict

def kMerize(text: str, k: int) -> List[str]:
    return sorted([text[i:i + k] for i in range(len(text) - k + 1)])
    
def constructGraph(kmers: List[str]):
    suffixToPrefix = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        suffixToPrefix[prefix].append(suffix)
    
    for suffix in suffixToPrefix:
        if suffix in suffixToPrefix:
            print("{} -> {}".format(suffix, ",".join(suffixToPrefix[suffix])))

    return

def part1(content: List[str]) -> None:
    k = int(content[0])
    dna = str(content[1])

    constructGraph(kMerize(dna,k))

def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python sol.py <filename>")
        return
    
    filename = sys.argv[1]
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]

    part1(content)

if __name__ == "__main__":
    main()