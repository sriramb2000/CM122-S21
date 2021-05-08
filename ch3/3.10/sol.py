import sys
from typing import List
import collections   
import itertools

# I did not write this: https://algorithmcode.wordpress.com/2014/11/18/string-reconstruction/

def part1(content: List[str]) -> None:
    in_kmer=content
    bkmer=[]

    bkmer=[False]*len(in_kmer)
    nkmer=[]
    pos = len(in_kmer[0])-1
    inout=[]
    for i,im in enumerate(in_kmer):
        sim = im[0:pos]
        eim = im[len(im)-pos:len(im)]
        ii = 0
        io = 0
        for j,jm in enumerate(in_kmer):
            sjm = jm[0:pos]
            ejm = jm[len(jm)-pos:len(jm)]
            if(sim==sjm):
                ii+=1
            if(sim==ejm):
                io+=1
        if(ii==1 and io==1):
            inout.append(im)
        else:
            if(not bkmer[i]):
                bkmer[i]=True
                nkmer.append(im)

    i=0
    while(len(inout)>0):
        iv=inout[i]
        liv = iv[0:pos]
        # print iv, liv
        for x in nkmer:
            if liv in x:
                nkmer.append(x+iv[pos:len(iv)])
                nkmer.remove(x)
                inout.remove(iv)
                i=-1
                break
        i+=1
            
        
    print(" ".join([im for i, im in enumerate(sorted(nkmer))]))

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