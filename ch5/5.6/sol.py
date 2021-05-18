import sys
from typing import List
from collections import defaultdict

def manhattanTourist(n: int, m: int, Down: List[List[int]], Right: List[List[int]]) -> int:
    dp = {}
    dp[(0,0)] = 0
    for i in range(1, n+1):
        dp[(i, 0)] = dp[(i-1, 0)] + Down[i-1][0]
    for j in range(1, m+1):
        dp[(0, j)] = dp[(0, j-1)] + Right[0][j-1]
    for i in range(1, n+1):
        for j in range(1, m+1):
            dp[(i,j)] = max(dp[(i-1, j)] + Down[i-1][j], dp[(i, j-1)] + Right[i][j-1])
    return dp[(n,m)]

def part1(content: List[str]):
    n, m = [int(x) for x in content[0].split(' ')]
    splitDex = content.index('-')
    down = []
    right = []
    for i in range(1, splitDex):
        down.append([int(x) for x in content[i].split(' ')])
    for j in range(splitDex+1, len(content)):
        right.append([int(x) for x in content[j].split(' ')])

    print(manhattanTourist(n, m, down, right))

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