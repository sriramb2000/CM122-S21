import sys
from typing import List

def decryptBWT(s):
    tokens = ["" for i in range(len(s))]
    for length in range(len(s)):
        tokens = [s[i]+tokens[i] for i in range(len(s))]
        tokens.sort()
        print(tokens)
    for s in tokens:
        if s[-1] == "$":
            return s
    

def part1(content: List[str]) -> None:
    text = content[0]
    print(decryptBWT(text))

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