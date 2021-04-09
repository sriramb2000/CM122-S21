import sys
from typing import List

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
        mismatches = 0
        for i, c in enumerate(p):
            if (c != text[offset+i]):
                k -= 1
                if k < 0:
                    return False
        return True
    
    def approxPatternPositions(p: str):
        positions = set()
        seedOffsetPairs = genSeeds(p)
        for (seed, offset) in seedOffsetPairs:
            candidate_poz = findSeed(seed)
            for pos in candidate_poz:
                pattern_pos = pos - offset
                if pattern_pos < 0 or (pattern_pos + len(p) > len(text)):
                    continue
                elif kMismatch(pattern_pos, p, d):
                    positions.add(pattern_pos)
        return list(positions)
    
    res = []
    for pattern in patterns:
        res += approxPatternPositions(pattern)
    res.sort()
    return res
    
def part1(content: List[str]) -> None:
    text = content[0]
    queries = content[1].split(" ")
    mismatches = int(content[2])
    
    positions = approximateMatching(text, queries, mismatches)
    print(" ".join([str(p) for p in positions]))

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