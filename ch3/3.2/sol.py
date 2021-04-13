import sys
from typing import List

def kMerize(text: str, k: int) -> List[str]:
    return sorted([text[i:i+k] for i in range(len(text)-k+1)])

def part1(content: List[str]) -> None:
    k = int(content[0])
    dna = str(content[1])

    print("\n".join(kMerize(dna, k)))

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