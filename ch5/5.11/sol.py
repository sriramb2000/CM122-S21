import sys
from typing import List
from collections import defaultdict

DOWN = 1
RIGHT = 2
DIAG = 3
SOURCE = 0

def editDistance(v: str, w: str):
    s = {}
    backtrack = {}
    s[0,0] = 0
    for i in range(1, len(v) + 1):
        s[(i, 0)] = i
        backtrack[(i, 0)] = DOWN
    for j in range(1, len(w) + 1):
        s[(0, j)] = j
        backtrack[(0, j)] = RIGHT

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s[(i, j)] = min(s[(i-1, j)], s[(i, j-1)], s[(i-1,j-1)]) + 1
            if v[i-1] == w[j-1]:
                s[i,j] = s[i-1, j-1]
            if s[(i,j)] == (s[(i-1, j)]):
                backtrack[(i,j)] = DOWN
            elif s[(i,j)] == (s[(i, j-1)]):
                backtrack[(i,j)] = RIGHT
            else:
                backtrack[(i,j)] = DIAG
    return s, backtrack

def fittingAlignment(v: str, w: str):
    indel_penalty = -1
    s = {}
    backtrack = {}
    for i in range(len(v) + 1):
        # Free taxi from source to all points in v
        s[(i, 0)] = 0
        backtrack[(i, 0)] = SOURCE
    for j in range(1, len(w) + 1):
        s[(0, j)] = s[0, j-1] + indel_penalty
        backtrack[(0, j)] = RIGHT
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            score = (1 if v[i-1] == w[j-1] else indel_penalty)
            s[(i, j)] = max(s[(i-1, j)] + indel_penalty, s[(i, j-1)] + indel_penalty, s[(i-1,j-1)] + score)
            if s[(i,j)] == (s[(i-1, j)] + indel_penalty):
                backtrack[(i,j)] = DOWN
            elif s[(i,j)] == (s[(i, j-1)] + indel_penalty):
                backtrack[(i,j)] = RIGHT
            else:
                backtrack[(i,j)] = DIAG
            # print("s[{}][{}] = {}, backtrack = {}".format(i, j, s[i,j], backtrack[i,j]))
    maxScore = s[0, len(w)]
    maxRow = 0
    for i in range(1, len(v) + 1):
        if s[i, len(w)] > maxScore:
            maxScore = s[i, len(w)]
            maxRow = i
    return s[maxRow,len(w)], backtrack, maxRow, len(w)

def overlapAlignment(v: str, w: str):
    indel_penalty = -2
    s = {}
    backtrack = {}
    for i in range(len(v) + 1):
        # Free taxi from source to all points in v
        s[(i, 0)] = 0
        backtrack[(i, 0)] = SOURCE
    for j in range(1, len(w) + 1):
        s[(0, j)] = s[0, j-1] + indel_penalty
        backtrack[(0, j)] = RIGHT
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            score = (1 if v[i-1] == w[j-1] else indel_penalty)
            s[(i, j)] = max(s[(i-1, j)] + indel_penalty, s[(i, j-1)] + indel_penalty, s[(i-1,j-1)] + score)
            if s[(i,j)] == (s[(i-1, j)] + indel_penalty):
                backtrack[(i,j)] = DOWN
            elif s[(i,j)] == (s[(i, j-1)] + indel_penalty):
                backtrack[(i,j)] = RIGHT
            else:
                backtrack[(i,j)] = DIAG
            # print("s[{}][{}] = {}, backtrack = {}".format(i, j, s[i,j], backtrack[i,j]))
    maxCol = 0
    maxScore = s[len(v), maxCol]
    
    for j in range(1, len(w) + 1):
        if s[len(v), j] > maxScore:
            maxScore = s[len(v), j]
            maxCol = j
    # print(len(v), maxCol)
    return s[len(v), maxCol], backtrack, len(v), maxCol


def outputLcs(backtrack, v: str, w: str, i: int, j: int):
    vRes, wRes = "", ""
    while (i != 0 or j != 0):
        if backtrack[(i,j)] == DOWN:
            vRes += v[i-1]
            wRes += "-"
            i -= 1
        elif backtrack[(i,j)] == RIGHT:
            vRes += "-"
            wRes += w[j-1]
            j -= 1
        elif backtrack[i,j] == SOURCE:
            break
        else:
            vRes += v[i-1]
            wRes += w[j-1]
            i -= 1
            j -= 1
    
    return (vRes[::-1], wRes[::-1])

def part1(content: List[str]):
    v = content[0]
    w = content[1]
    scores, backtrack = editDistance(v, w)
    print(scores[len(v), len(w)])

def part2(content: List[str]):
    v = content[0]
    w = content[1]
    maxScore, backtrack, maxV, maxW = fittingAlignment(v, w)
    v_, w_ = outputLcs(backtrack, v, w, maxV, maxW)
    print(maxScore)
    print(v_)
    print(w_)

def part3(content: List[str]):
    v = content[0]
    w = content[1]
    maxScore, backtrack, maxV, maxW = overlapAlignment(v, w)
    v_, w_ = outputLcs(backtrack, v, w, maxV, maxW)
    print(maxScore)
    print(v_)
    print(w_)
    
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
    elif sys.argv[1] == "3":
        part3(content)

if __name__ == "__main__":
    main()