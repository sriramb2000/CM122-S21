import sys
from typing import List

def part1(content: List[str]) -> None:
    text = content[0]
    suffixs = []
    for length in range(1,len(text)+1):
        suffix = text[-length:]
        suffixs.append((suffix, len(text) - length))
    suffixs.sort(key = lambda x: x[0])
    sortedSuffix = [str(x[1]) for x in suffixs]
    print(", ".join(sortedSuffix))

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