from typing import Optional, List
import sys

class TrieNode(object):
    id = 0
    def __init__(self, char: Optional[str]):
        self.char = char
        self.children = {}
        self.id = TrieNode.id
        self.isEnd = False
        TrieNode.id += 1

    def addWord(self, word: str) -> None:
        if len(word) > 0:
            letter = word[0]
            if letter not in self.children:
                self.children[letter] = TrieNode(letter)
            self.children[letter].addWord(word[1:])
        else:
            self.isEnd = True
    
    def lookupWord(self, word: str) -> bool:
        if self.isEnd:
            return True
        if len(word) == 0:
            return False
        letter = word[0]
        return (letter in self.children) and (self.children[letter].lookupWord(word[1:]))

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

def part1(content: List[str]) -> None:
    PrintTrie(TrieConstruction(content))

def part2(content: List[str]) -> None:
    text = content[0]
    patterns = content[1:]

    trie = TrieConstruction(patterns)
    i = 0
    while text:
        if trie.lookupWord(text):
            print(i, end=" ")
        text = text[1:]
        i += 1
    print()

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
    
