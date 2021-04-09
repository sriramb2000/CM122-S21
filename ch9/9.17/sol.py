import sys
from typing import List

def partialSuffix(text, k) -> None:
    suffixs = []
    for length in range(1,len(text)+1):
        suffix = text[-length:]
        suffixs.append((suffix, len(text) - length))
    suffixs.sort(key = lambda x: x[0])

    res = []
    for i in range(len(suffixs)):
        if suffixs[i][1] % k == 0:
            res.append((i, suffixs[i][1]))
    return res

def part1(content: List[str]) -> None:
    text = content[0]
    k = int(content[1])
    
    partial = partialSuffix(text, k)
    for p in partial:
        print("{},{}".format(p[0], p[1]))


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