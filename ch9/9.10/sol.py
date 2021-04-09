import sys
from typing import List

def bwtQuery(text, query) -> int:
    temp = [i[0] for i in sorted(enumerate(list(text)), key=lambda x:x[1])]
    l2f = {}
    for f, l in enumerate(temp):
        l2f[l] = f
    first = 0
    last = len(text) - 1
    while query:
        cur = query[-1]
        tf, tl = text.find(cur, first, last+1), text.rfind(cur, first, last+1)
        if tf == -1 or tl == -1:
            return 0
        first, last = l2f[tf], l2f[tl]
        query = query[:-1]

    return last - first + 1

def part1(content: List[str]) -> None:
    text = content[0]
    queryStrings = content[1].split(" ")
    print(" ".join([str(bwtQuery(text, q)) for q in queryStrings]))

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