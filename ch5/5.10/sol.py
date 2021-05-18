import sys
from typing import List
from collections import defaultdict

DOWN = 0
RIGHT = 1
DIAG = 2
SOURCE = 3

INDEL_PENALTY = -5

def lcsBacktrack(v: str, w: str, score):
    s = {}
    backtrack = {}
    for i in range(len(v) + 1):
        s[(i, 0)] = i * INDEL_PENALTY
        backtrack[(i, 0)] = DOWN
    for j in range(len(w) + 1):
        s[(0, j)] = j * INDEL_PENALTY
        backtrack[(0, j)] = RIGHT
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s[(i, j)] = max(s[(i-1, j)] + INDEL_PENALTY, s[(i, j-1)] + INDEL_PENALTY, s[(i-1,j-1)] + score[v[i-1] + w[j-1]])
            if s[(i,j)] == (s[(i-1, j)] + INDEL_PENALTY):
                backtrack[(i,j)] = DOWN
            elif s[(i,j)] == (s[(i, j-1)] + INDEL_PENALTY):
                backtrack[(i,j)] = RIGHT
            else:
                backtrack[(i,j)] = DIAG
    return s, backtrack

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
        else:
            vRes += v[i-1]
            wRes += w[j-1]
            i -= 1
            j -= 1
    
    return (vRes[::-1], wRes[::-1])

def lcsBacktrack2(v: str, w: str, score):
    s = {}
    backtrack = {}
    for i in range(len(v) + 1):
        s[(i, 0)] = 0
        backtrack[(i, 0)] = DOWN
    for j in range(len(w) + 1):
        s[(0, j)] = 0
        backtrack[(0, j)] = DOWN
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s[(i, j)] = max(s[(i-1, j)] + INDEL_PENALTY, s[(i, j-1)] + INDEL_PENALTY, s[(i-1,j-1)] + score[v[i-1] + w[j-1]], 0)
            if s[(i,j)] == (s[(i-1, j)] + INDEL_PENALTY):
                backtrack[(i,j)] = DOWN
            elif s[(i,j)] == (s[(i, j-1)] + INDEL_PENALTY):
                backtrack[(i,j)] = RIGHT
            elif 0 == s[(i,j)]:
                backtrack[(i,j)] = SOURCE
            else:
                backtrack[(i,j)] = DIAG
    return s, backtrack

def outputLcs2(backtrack, v: str, w: str, i: int, j: int):
    vRes, wRes = "", ""
    while (i != 0 or j != 0):
        if i == 0:
            vRes += '-'*j
            wRes += w[:j]
            break
        elif j == 0:
            wRes += '-'*i
            vRes += v[:i]
            break
        elif backtrack[(i,j)] == DOWN:
            vRes += v[i-1]
            wRes += "-"
            i -= 1
        elif backtrack[(i,j)] == RIGHT:
            vRes += "-"
            wRes += w[j-1]
            j -= 1
        elif backtrack[(i,j)] == DIAG:
            vRes += v[i-1]
            wRes += w[j-1]
            i -= 1
            j -= 1
        else:
            break
    
    return (vRes[::-1], wRes[::-1])

blosum62="""
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
"""

pam250 ="""
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  2 -2  0  0 -3  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3
C -2 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0
D  0 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4
E  0 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4
F -3 -4 -6 -5  9 -5 -2  1 -5  2  0 -3 -5 -5 -4 -3 -3 -1  0  7
G  1 -3  1  0 -5  5 -2 -3 -2 -4 -3  0  0 -1 -3  1  0 -1 -7 -5
H -1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0
I -1 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1
K -1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4
L -2 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1
M -1 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2
N  0 -4  2  1 -3  0  2 -2  1 -3 -2  2  0  1  0  1  0 -2 -4 -2
P  1 -3 -1 -1 -5  0  0 -2 -1 -3 -2  0  6  0  0  1  0 -1 -6 -5
Q  0 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4
R -2 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4
S  1  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3
T  1 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3
V  0 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2
W -6 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0
Y -3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10
"""

def getScoringMatrix(scores):
    score = [line for line in scores.splitlines() if line]
    Amino_acids = score[0].split()
    score_dict = {}
    for i in range(len(Amino_acids)):
        for j in range(len(Amino_acids)):
            score_dict[Amino_acids[i]+Amino_acids[j]]=int(score[i+1].split()[j+1])

    return score_dict

def part1(content: List[str]):
    v = content[0]
    w = content[1]
    scores,backtrack=lcsBacktrack(v,w, getScoringMatrix(blosum62))
    v_, w_ = outputLcs(backtrack, v, w, len(v), len(w))
    print(scores[len(v),len(w)])
    print(v_)
    print(w_)

def part2(content: List[str]):
    v = content[0]
    w = content[1]
    scores, backtrack = lcsBacktrack2(v, w, getScoringMatrix(pam250))
    max_score = max(list(scores.values()))
    max_index = None
    for key in scores:
        if max_score == scores[key]:
            max_index = key
    v_, w_ = outputLcs2(backtrack, v, w, max_index[0], max_index[1])
    print(max_score)
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

if __name__ == "__main__":
    main()