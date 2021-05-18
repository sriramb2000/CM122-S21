import sys
from typing import List
from collections import defaultdict

DOWN = 0
RIGHT = 1
DIAG = 2

def lcsBacktrack(v: str, w: str):
    s = defaultdict(int)
    backtrack = {}
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            match = 0
            if v[i-1] == w[j-1]:
                match = 1
            s[(i, j)] = max(s[(i-1, j)], s[(i, j-1)], s[(i-1,j-1)] + match)
            if s[(i,j)] == s[(i-1, j)]:
                backtrack[(i,j)] = DOWN
            elif s[(i,j)] == s[(i, j-1)]:
                backtrack[(i,j)] = RIGHT
            else:
                backtrack[(i,j)] = DIAG

    return backtrack

def outputLcs(backtrack, v: str, i: int, j: int):
    res = ""
    while (i != 0 and j != 0):
        if backtrack[(i,j)] == DOWN:
            i -= 1
        elif backtrack[(i,j)] == RIGHT:
            j -= 1
        else:
            res += v[i-1]
            i -= 1
            j -= 1
    
    return res[::-1]

def buildPredecessorGraph(edges: List[str]):
    predecessors = defaultdict(lambda: defaultdict(int))
    successors = defaultdict(lambda: defaultdict(int))
    for edge in edges:
        e, weight = edge.split(":")
        weight = int(weight)
        parent, child = e.split("->")
        successors[parent][child] = weight
        predecessors[child][parent] = weight
    
    return (predecessors, successors)

def longestPath(source, sink, successors, predecessors):
    dp = defaultdict(lambda: (0, None))
    dp[source] = (0, None)
    queue = list(successors[source].keys())
    while queue:
        curNode = queue.pop(0)
        # Get all parent weights, find the max
        maxParent = max([(dp[parent][0] + predecessors[curNode][parent], parent) for parent in list(predecessors[curNode].keys())], key = lambda i : i[0])
        dp[curNode] = (maxParent[0], maxParent[1])
        # push children, if sink, push None
        if curNode == sink:
            continue
        queue += list(successors[curNode].keys())

    path = []
    cur = sink
    while cur:
        path.append(cur)
        cur = dp[cur][1]
    path.reverse()

    return (dp[sink][0], path)


def part1(content: List[str]):
    v = content[0]
    w = content[1]
    print(outputLcs(lcsBacktrack(v, w), v, len(v), len(w)))

def part2(content: List[str]):
    source = content[0]
    sink = content[1]
    p, s = buildPredecessorGraph(content[2:])
    length, path = longestPath(source, sink, s, p)
    print(length)
    print("->".join(path))
    
def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: python sol.py 1|2 <filename>")
        return
    
    filename = sys.argv[2]
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]

    if sys.argv[1] == "1":
        part1(content)
    elif sys.argv[1] == "2":
        part2(content)

if __name__ == "__main__":
    main()