import sys
from typing import List
from collections import defaultdict

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

def bwtQuery(text, query) -> int:
    temp = [i[0] for i in sorted(enumerate(list(text)), key=lambda x:x[1])]
    l2f = {}
    for f, l in enumerate(temp):
        l2f[l] = f
    first = 0
    last = len(text) - 1
    while query:
        cur = query[-1]
        query = query[:-1]
        tf, tl = text.find(cur, first, last+1), text.rfind(cur, first, last+1)
        if tf == -1 or tl == -1:
            return 0
        first, last = l2f[tf], l2f[tl]
        
    return last - first + 1

def betterBwtQuery(text, query) -> int:
    first_col = ''.join(sorted(text))
    first = 0
    last = len(text) - 1
    while first <= last:
        if query:
            cur = query[-1]
            query = query[:-1]
            if text.find(cur, first, last+1) >= 0:
                first = first_col.find(cur) + text.count(cur, 0, first)
                last = first_col.find(cur) + text.count(cur, 0, last + 1) - 1
            else:
                return 0
        else:
            return last - first + 1

def findPatternOccurrences(text: str, patterns: List[str]):
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

    def findPattern(p: str):
        top = 0
        bottom = len(text) - 1
        while top <= bottom:
            if p:
                symbol = p[-1]
                p = p[:-1]
                if count_symbol[bottom+1].get(symbol, 0) > count_symbol[top].get(symbol, 0):
                    top = first_occurrence[symbol] + count_symbol[top].get(symbol,0)
                    bottom = first_occurrence[symbol] + count_symbol[bottom+1].get(symbol, 0) - 1
                else:
                    break
            else:
                return [suffix[i] for i in range(top, bottom+1)]
        return []
    
    res = []
    for pattern in patterns:
        res += findPattern(pattern)
    res.sort()
    return res

def countSymbol(symbol: str, text: List, i: int) -> int:
    return text[:i].count(symbol)

def getFirstOccurences(bwtString: str):
    first = sorted(list(bwtString))
    first_occurrences = {}
    for i in range(len(first)):
        if first[i] not in first_occurrences:
            first_occurrences[first[i]] = i
    return first_occurrences

def part1(content: List[str]) -> None:
    text = content[0]
    queryStrings = content[1].split(" ")
    print(" ".join([str(bwtQuery(text, q)) for q in queryStrings]))

def part2(content: List[str]) -> None:
    text = content[0]
    queryStrings = content[1:]
    print(" ".join([str(betterBwtQuery(text, q)) for q in queryStrings]))

def part3(content: List[str]) -> None:
    text = content[0]
    queryStrings = content[1:]
    print(" ".join([str(x) for x in findPatternOccurrences(text, queryStrings)]))

def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: python sol.py 1|2|3 <filename>")
        return
    
    filename = sys.argv[2]
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content] 

    if sys.argv[1] == "1":
        part1(content)
    elif sys.argv[1] == "2":
        part2(content)
    elif sys.argv[1] == "3":
        part3(content)

if __name__ == "__main__":
    main()