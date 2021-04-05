from typing import Optional, List
import sys

class TrieNode(object):
    id = 0
    def __init__(self, char: Optional[str]):
        self.char = char
        self.children = {}
        self.id = TrieNode.id
        TrieNode.id += 1

    def addWord(self, word: str):
        if len(word) > 0:
            letter = word[0]
            if letter not in self.children:
                self.children[letter] = TrieNode(letter)
            self.children[letter].addWord(word[1:])
        return

def TrieConstruction(patterns: List[str]) -> TrieNode:
    root = TrieNode(None)
    for pattern in patterns:
        root.addWord(pattern)
    return root

def PrintTrie(root: TrieNode) -> None:
    stack = [root]
    while stack:
        current = stack.pop(0)
        for letter in current.children:
            child = current.children[letter]
            print('{:d}->{:d}:{:s}'.format(current.id, child.id, letter))
            stack.append(child)
    return

def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python sol.py <filename>")
        return
    
    filename = sys.argv[1]
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content] 

    PrintTrie(TrieConstruction(content))

if __name__ == "__main__":
    main()
    
