from typing import Optional, List
from collections import defaultdict
from itertools import product
import random
import sys

def constructGraph(kmers: List[str]):
    suffixToPrefix = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffixToPrefix[prefix].append(kmer)
    
    for kmer in kmers:
        suffix = kmer[1:]
        if suffix in suffixToPrefix:
            print("{} -> {}".format(kmer, ",".join(suffixToPrefix[suffix])))

    return

def createGraph(content: List[str]):
    edges = set()
    nodeMap = defaultdict(set)
    appearanceAsChild = defaultdict(int)
    for line in content:
        parent, children = [item.strip() for item in line.split('->')]
        for child in [c.strip() for c in children.split(',')]:
            edge = (parent, child)
            nodeMap[parent].add(child)
            edges.add(edge)
            appearanceAsChild[child] += 1
    return (edges, nodeMap, appearanceAsChild)

def findStartAndEnd(nodeMap: defaultdict(set), appearanceAsChild: defaultdict(int)):
    start, end = None, None
    keys = set(nodeMap.keys()) | set(appearanceAsChild.keys())
    for key in keys:
        numChildren = len(nodeMap[key])
        numParents = appearanceAsChild[key]
        if numChildren > numParents:
            start = key
        elif numParents > numChildren:
            end = key

    return (start, end)

def findCycle(edges: set, nodeMap: defaultdict(set), seed: Optional[str]=None):
    parent, child = None, None
    if seed:
        parent = seed
        child = nodeMap[parent].pop()
    else:
        parent,child = edges.pop()
        nodeMap[parent].remove(child)

    if len(nodeMap[parent]) == 0:
        del nodeMap[parent]

    cycle = [parent]
    while len(edges) != 0:
        try:
            cycle.append(child)
            # Remove edge from graph
            edges.remove((parent,child))
            if len(nodeMap[parent]) == 0:
                del nodeMap[parent]

            parent = child
            child = nodeMap[parent].pop()
        except:
            if len(nodeMap[parent]) == 0:
                del nodeMap[parent]
            break

    return (cycle, edges, nodeMap)

def mergeCycles(oldCycle: List[str], newCycle: List[str]):
    merge_index = oldCycle.index(newCycle[0])
    return oldCycle[:merge_index]+newCycle+oldCycle[merge_index+1:]

def findEulerianCycle(edges: set, nodeMap: defaultdict(set), seed: Optional[str] = None):
    cycle, unvisitedEdges, unvisitedNodeMap = findCycle(edges.copy(), nodeMap.copy(), seed)
    while len(unvisitedEdges) != 0:
        potential = [item for item in cycle if item in unvisitedNodeMap]
        seed = random.choice(potential)
        newCycle, unvisitedEdges, unvisitedNodeMap = findCycle(unvisitedEdges, unvisitedNodeMap, seed)
        cycle = mergeCycles(cycle,newCycle)
    return cycle

def createEulerianPathFromCycle(cycle: List[str], start: str, end: str):
    if cycle[0] == start and cycle[-1] == end:
        return cycle
    cycleString = '->'.join(cycle)
    endStartEdge = "{}->{}".format(end, start)
    edgeIndexStart = cycleString.rfind(endStartEdge)
    edgeIndexEnd = edgeIndexStart + len(endStartEdge)

    newCycleStart = cycleString[edgeIndexEnd - len(start):-len(start)]
    newCycleEnd = cycleString[:edgeIndexStart + len(end)]
    newCycle = "{}{}".format(newCycleStart, newCycleEnd)

    return newCycle.split("->")

def findEulerianPath(edges: List[str], nodeMap: defaultdict(set), appearanceAsChild: defaultdict(int)):
    start, end = findStartAndEnd(nodeMap, appearanceAsChild)
    # Add an edge to make it cycle
    edges.add((end, start))
    nodeMap[end].add(start)

    # Find the cycle
    ec = findEulerianCycle(edges, nodeMap, start)
    return createEulerianPathFromCycle(ec, start, end)

def constructKmerAdjacencyList(kmers: List[str]) -> List[str]:
    suffixToPrefix = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        suffixToPrefix[prefix].append(suffix)
    
    adjList = []

    for suffix in suffixToPrefix:
        if suffix in suffixToPrefix:
            adjList.append("{} -> {}".format(suffix, ",".join(suffixToPrefix[suffix])))

    return adjList

def reconstructGenome(kmers: List[str]) -> str:
    adjList = constructKmerAdjacencyList(kmers)
    edges, nodeMap, appearanceAsChild = createGraph(adjList)
    ep = findEulerianPath(edges, nodeMap, appearanceAsChild)

    return pathToGenome(ep)

def pathToGenome(eulerianPath: List[str]) -> str:
    return eulerianPath[0] + "".join([s[-1] for s in eulerianPath[1:]])

def prefixSuffixMerge(suffixPrefixList: List[str]) -> List[str]:
    res = []
    for i in range(len(suffixPrefixList) - 1):
        res.append(suffixPrefixList[i] + suffixPrefixList[i+1][-1])
    return suffixPrefixList

def generateKBinaryStrings(k: int) -> List[str]:
    return ["".join(w) for w in product(['0','1'], repeat=k)]

def generateKCircularString(k: int) -> str:
    adjList = constructKmerAdjacencyList(generateKBinaryStrings(k))
    edges, nodeMap, appearanceAsChild = createGraph(adjList)
    ec = findEulerianCycle(edges, nodeMap)

    pathString = ''.join([item[-1] for item in ec[1:]])
    return pathString

def part1(content: List[str]) -> None:
    edges, nodeMap, _ = createGraph(content)
    eulerianCycle = findEulerianCycle(edges, nodeMap)
    print("->".join(eulerianCycle))

def part2(content: List[str]) -> None:
    edges, nodeMap, appearanceAsChild = createGraph(content)
    ep = findEulerianPath(edges, nodeMap, appearanceAsChild)

    print("->".join(ep))

def part3(content: List[str]) -> None:
    k = int(content[0])
    print(reconstructGenome(content[1:]))

def part4(content: List[str]) -> None:
    k = int(content[0])
    print(generateKCircularString(k))

def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: python sol.py 1|2|3|4 <filename>")
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
    elif sys.argv[1] == "4":
        part4(content)

if __name__ == "__main__":
    main()