import sys
from typing import List

def rotate(strg, n):
    return strg[n:] + strg[:n]

def part1(content: List[str]) -> None:
    text = content[0]
    rotated = sorted([rotate(text, i) for i in range(len(text))])
    print("".join([s[-1] for s in rotated]))

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