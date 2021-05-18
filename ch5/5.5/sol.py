import sys
from typing import List
from collections import defaultdict

def dpChange(target: int, denoms: List[int]) -> int:
    dp = defaultdict(lambda: float('inf'))
    dp[0] = 0
    for m in range(target+1):
        for coin in denoms:
            dp[m] = (dp[m - coin] + 1) if (m >= coin and ((dp[m - coin] + 1) < dp[m])) else dp[m]

    return dp[target]

def part1(content: List[str]):
    target = int(content[0])
    denoms = [int(d) for d in content[1].split(',')]

    print(dpChange(target, denoms))

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