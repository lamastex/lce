import sys
import fileinput
from decimal import Decimal 

'''takes in from stdin the branch-lengths in newick-tree from ms via newickutils and outputs epoch-times
   example useage and output:
   
$ ms 9 3 -T  | tail -n +4 | grep -v // | nw_distance -s a -t - | python time2epoch.py -
0.001	0.109	0.007	0.038	0.165	0.059	0.047	0.192
0.005	0.013	0.006	0.066	0.026	0.09	0.04	0.652
0.001	0.014	0.03	0.029	0.051	0.05	0.136	0.042

'''
#import math

n = 0
maxnsam=1000 # make max sample size 1000
m = [0.0 for i in range(maxnsam)]
nsam = 0 # unknown now

for l in fileinput.input():
    a = l.strip().split('\t')
    print("~~~~~~~~~~"); 
    print(len(a), a)
    a = set(a)
    print(len(a), a)
    a = sorted([float(i) for i in a])
    a.reverse()
    print(len(a), a)
    b = [s - t for s, t in zip(a, a[1:])]; # epoch times when there are n, n-1, ..., 2 lineages
    if (n==0): 
        nsam=len(b)
    print(n,nsam)
    print("....")
    print(nsam,b)
    assert(nsam==len(b))
    n = n + 1
    for i in range(len(b)):
        m[i] = m[i] + (b[i] - m[i]) / n

    print(b)
    #print('\t'.join([str(i) for i in b]))
    print(m[0:nsam])

print('\t'.join([str(i) for i in m[0:nsam]]))
