from typing import List, Optional, Tuple
import sys

class SuffixTree(object):

    class Node(object):
        def __init__(self, lab):
            self.lab = lab # label on path leading to this node
            self.out = {}  # outgoing edges; maps characters to nodes
            self.index = -1

    def __init__(self, s):
        """ Make suffix tree, without suffix links, from s in quadratic time
            and linear space """
        suffix=[]
        self.suffix=suffix
        self.original_string = s
        self.root = self.Node(None)
        self.root.out[s[0]] = self.Node(s) # trie for just longest suf
        # add the rest of the suffixes, from longest to shortest
        for i in range(1, len(s)):
            # start at root; we’ll walk down as far as we can go
            cur = self.root
            j = i
            while j < len(s):
                if s[j] in cur.out:

                    child = cur.out[s[j]]
                    lab = child.lab
                    # Walk along edge until we exhaust edge label or
                    # until we mismatch
                    k = j+1 
                    while k-j < len(lab) and s[k] == lab[k-j]:
                        k += 1
                    if k-j == len(lab):
                        cur = child # we exhausted the edge 
                        suffix.append(child.lab)
                        j = k
                    else:
                        # we fell off in middle of edge
                        cExist, cNew = lab[k-j], s[k]
                        # create “mid”: new node bisecting edge
                        mid = self.Node(lab[:k-j])
                        mid.out[cNew] = self.Node(s[k:])
                        # original child becomes mid’s child
                        mid.out[cExist] = child
                        # original child’s label is curtailed
                        child.lab = lab[k-j:]

                        # mid becomes new child of original parent
                        cur.out[s[j]] = mid
                else:
                    # Fell off tree at a node: make new edge hanging off it
                    cur.out[s[j]] = self.Node(s[j:])

    def followPath(self, s):
        """ Follow path given by s.  If we fall off tree, return None.  If we
            finish mid-edge, return (node, offset) where 'node' is child and
            'offset' is label offset.  If we finish on a node, return (node,
            None). """
        cur = self.root
        i = 0
        while i < len(s):
            c = s[i]
            if c not in cur.out:
                return (None, None) # fell off at a node
            child = cur.out[s[i]]
            lab = child.lab
            j = i+1
            while j-i < len(lab) and j < len(s) and s[j] == lab[j-i]:
                j += 1
            if j-i == len(lab):
                cur = child # exhausted edge
                i = j
            elif j == len(s):
                return (child, j-i) # exhausted query string in middle of edge
            else:
                return (None, None) # fell off in the middle of the edge
        return (cur, None) # exhausted query string at internal node

    def hasSubstring(self, s):
        """ Return true iff s appears as a substring """
        node, off = self.followPath(s)
        return node

    def hasSuffix(self, s):
        """ Return true iff s is a suffix """
        node, off = self.followPath(s)
        if node is None:
            return False # fell off the tree
        if off is None:
            # finished on top of a node
            return '$' in node.out
        else:
            # finished at offset 'off' within an edge leading to 'node'
            return node.lab[off] == '$'
    
    # def markNodes(self, node):
    #     ret = -1
    #     if len(node.out) > 1:
    #         for key in node.out:
    #             next = node.out[key]
    #             ret = self.markNodes(next)
    #             if (node.index == -1):
    #                 node.index = ret
    #             elif ((node.index == -2 and ret == -3) or (node.index == -3 and ret == -2) or (node.index == -4)):
    #                 node.index = -4
    #     elif (node.index > -1 and node.index < self.splitdex):
    #         return -2
    #     elif (node.index >= self.splitdex):
    #         return -3

    #     return node.index

def deepestInternal(node, builder):
    if len(node.out) > 1:
        temps = []
        for key in node.out:
            next = node.out[key]
            temps.append(deepestInternal(next, builder + (node.lab if node.lab else "")))
        temps.append(builder)
        builder = max(temps, key = len)
    return builder



def part1(content: List[str]) -> None:
    search_str = content[0]
    tree = SuffixTree(search_str)
    tree.printEdges()

def part2(content: List[str]) -> None:
    search_str = content[0]
    tree = SuffixTree(search_str)
    print(deepestInternal(tree.root, ""))

def part3(content: List[str]) -> None:
    s1, s2 = content[0], content[1]
    tree1, tree2 = SuffixTree(s1 + "$"), SuffixTree(s2 + "$")
    test_str = s1 if (len(s1) < len(s2)) else s2
    substrs = [test_str[j: j+i] for i in range(len(test_str), 0, -1)
          for j in range(0, len(test_str) - i + 1)]
    for s in substrs:
        if tree1.hasSubstring(s) and tree2.hasSubstring(s):
            print(s)
            break

def part4(content: List[str]) -> None:
    s1, s2 = content[0], content[1]
    tree1, tree2 = SuffixTree(s1 + "$"), SuffixTree(s2 + "$")
    for length in range(max(len(s1), len(s2))):
        if length <= len(s1):
            for j in range(0, len(s1) - length + 1):
                sub = s1[j:j+length]
                if (tree1.hasSubstring(sub) and not tree2.hasSubstring(sub)) or (tree1.hasSubstring(sub) and not tree2.hasSubstring(sub)):
                    print(sub)
                    return
        if length <= len(s2):
            for j in range(0, len(s2) - length + 1):
                sub = s2[j:j+length]
                if (tree1.hasSubstring(sub) and not tree2.hasSubstring(sub)) or (tree1.hasSubstring(sub) and not tree2.hasSubstring(sub)):
                    print(sub)
                    return

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