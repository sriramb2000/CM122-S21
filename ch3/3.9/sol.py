import sys
from typing import List
import collections
import itertools

# I did not write this: https://algorithmcode.wordpress.com/2014/11/18/string-reconstruction/

def part1(content: List[str]) -> None:
    line = 1
    in_pkmer=[]
    in_skmer=[]
    in_kmer={}
    bkmer=[]
    in_distance = 0
    okmer=[]
    i=0
    _, in_distance = [int(token) for token in content[0].split(" ")]
    for ldata in content[1:]:
        ls = ldata.split("|")
        in_kmer[i]=[ls[0],ls[1]]
        i+=1

        line+=1

    bkmer=[False]*len(in_kmer)

    i=0
    pos = len(in_kmer[0][0])-1
    while(bkmer.count(False)>1):
        if(not bkmer[i]):
            ikmer = in_kmer[i][0]
            ikmer1 = in_kmer[i][1]
            ek = ikmer[len(ikmer)-pos:len(ikmer)]
            ek1 = ikmer1[len(ikmer1)-pos:len(ikmer1)]
            for j in range(0,len(in_kmer)):
                if(not bkmer[j]):
                    jkmer = in_kmer[j][0]
                    jkmer1 = in_kmer[j][1]
                    sk = jkmer[0:pos]
                    sk1 = jkmer1[0:pos]
                    if(ek==sk and ek1==sk1):
                        npref =ikmer+jkmer[pos:len(jkmer)]
                        nsuff =in_kmer[i][1]+in_kmer[j][1][pos:len(jkmer)]
                        in_kmer[i][0]=npref
                        in_kmer[i][1]=nsuff
                        bkmer[j]=True
                        i=-1
                        break
        i+=1

    bpos = bkmer.index(False)
    pre = in_kmer[bpos][0]
    suf = in_kmer[bpos][1]
    res = pre+suf[len(suf)-(pos+1+in_distance):len(suf)]
    print(res)

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