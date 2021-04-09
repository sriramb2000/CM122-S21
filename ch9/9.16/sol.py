import sys
from typing import List

class Node(object):
    def __init__(self, children):
        self.children = children
        self.color = "gray"

def colorNodes(nodeDic) -> None:
    def nodeColor(node):
        if node.color == "gray":
            childColors = set([nodeColor(nodeDic[c]) for c in node.children])
            if len(childColors) > 1:
                node.color = "purple"
            else:
                node.color = childColors.pop()
        return node.color
    
    for node in nodeDic:
        nodeColor(nodeDic[node])

def part1(content: List[str]) -> None:
    splitdex = content.index("-")
    nodeDic = {}
    for row in content[:splitdex]:
        tokens = [x.strip() for x in row.split("->")]
        children = []
        if tokens[1] != "{}":
            children = tokens[1].split(",")
        nodeDic[tokens[0]] = Node(children)
    for row in content[splitdex+1:]:
        tokens = [x.strip() for x in row.split(":")]
        nodeDic[tokens[0]].color = tokens[1]

    colorNodes(nodeDic)

    keylist = list(nodeDic.keys())
    keylist.sort()
    for node in keylist:
        print("{}: {}".format(node, nodeDic[node].color))


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