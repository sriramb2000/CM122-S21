import sys
from typing import List
from collections import defaultdict

# GRAPH LEVELS
BOTTOM = 0
MIDDLE = 1
TOP = 2

# BACKTRACK DIRECTIONS
DOWN = 1
RIGHT = 2
DIAG = 3
SOURCE = 0

def gapAlignment(v: str, w: str, score, sigma = 11, epsilon = 1):
    sB = defaultdict(lambda: float('-inf'))
    sM = defaultdict(lambda: float('-inf'))
    sU = defaultdict(lambda: float('-inf'))
    s = [sB, sM, sU]
    backtrackB = defaultdict(int)
    backtrackM = defaultdict(int)
    backtrackU = defaultdict(int)
    backtrack = [backtrackB, backtrackM, backtrackU]

    s[BOTTOM][0,0] = 0
    s[MIDDLE][0,0] = 0
    s[TOP][0,0] = 0
    s[BOTTOM][0,1] = -sigma
    s[MIDDLE][0,1] = -sigma
    for i in range(2, len(v) + 1):
        s[BOTTOM][i, 0] = s[BOTTOM][i-1, 0] - epsilon
        s[MIDDLE][i, 0] = s[0][i, 0]
    s[TOP][0,1] = -sigma
    s[MIDDLE][0,1] = -sigma
    for j in range(2, len(w) + 1):
        s[TOP][0, j] = s[TOP][0, j-1] - epsilon
        s[MIDDLE][0, j] = s[TOP][0, j]
    
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s[BOTTOM][i, j] = max(s[BOTTOM][i-1, j] - epsilon, s[MIDDLE][i-1, j] - sigma)
            if s[BOTTOM][i,j] == s[MIDDLE][i-1, j] - sigma:
                backtrack[BOTTOM][i, j] = MIDDLE

            s[TOP][i, j] = max(s[TOP][i, j-1] - epsilon, s[MIDDLE][i, j-1] - sigma)
            if s[TOP][i, j] == s[MIDDLE][i, j-1] - sigma:
                backtrack[TOP][i, j] = MIDDLE

            s[MIDDLE][i, j] = max(s[BOTTOM][i,j], s[TOP][i,j], s[MIDDLE][i-1, j-1] + score[v[i-1] + w[j-1]])
            if s[MIDDLE][i, j] == s[MIDDLE][i-1, j-1] + score[v[i-1] + w[j-1]]:
                backtrack[MIDDLE][i, j] = MIDDLE
            elif s[MIDDLE][i, j] == s[TOP][i, j]:
                backtrack[MIDDLE][i, j] = TOP
            
    return s[MIDDLE][len(v), len(w)], backtrack, len(v), len(w)

def outputLcs(backtrack, v: str, w: str, i: int, j: int):
    vRes, wRes = "", ""
    lvl = MIDDLE
    while (i != 0 or j != 0):
        if lvl == BOTTOM:
            if backtrack[lvl][i,j] == MIDDLE:
                lvl = MIDDLE
            vRes += v[i-1]
            wRes += "-"
            i -= 1
        elif lvl == TOP:
            if backtrack[lvl][i,j] == MIDDLE:
                lvl = MIDDLE
            vRes += "-"
            wRes += w[j-1]
            j -= 1
            pass
        else:
            if backtrack[lvl][i,j] == MIDDLE:
                vRes += v[i-1]
                wRes += w[j-1]
                i -= 1
                j -= 1
            else:
                lvl = backtrack[lvl][i,j]
            pass
    
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

def getScoringMatrix():
    score = [line for line in blosum62.splitlines() if line]
    Amino_acids = score[0].split()
    score_dict = {}
    for i in range(len(Amino_acids)):
        for j in range(len(Amino_acids)):
            score_dict[Amino_acids[i]+Amino_acids[j]]=int(score[i+1].split()[j+1])

    return score_dict

def part1(content: List[str]):
    v = content[0]
    w = content[1]
    maxScore, backtrack, maxV, maxW = gapAlignment(v, w, getScoringMatrix())
    v_, w_ = outputLcs(backtrack, v, w, maxV, maxW)
    print(maxScore)
    print(v_)
    print(w_)

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