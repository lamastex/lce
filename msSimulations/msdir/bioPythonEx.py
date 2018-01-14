from Bio import Phylo
from cStringIO import StringIO
import numpy as np

'''to obtain various statistics of the tree output from ms
example:
./ms 15 1000 -T  | tail -n +4 | grep -v //  > out && python bioPythonEx.py 
Colles Index - mean and std :
0.137291208791 0.04952766606
mean of relative lengths of lineages leading to singleton, doubleton, ..., (n-1)-ton ineages and total tree size
[ 0.03551693  0.05677844  0.07127058  0.09012895  0.10335559  0.09314382
  0.08704998  0.08034704  0.0745306   0.08074013  0.06361253  0.05812927
  0.04279431  0.06260186]
1.43528098563
'''

# see http://biopython.org/DIST/docs/tutorial/Tutorial.html

trees = Phylo.parse("out", "newick")
collesIndices=[]
lineageLengths=[]
for tree in trees:
    n = tree.root.count_terminals()
    #print(tree)
    collesIndex=0.0;
    #Phylo.draw_ascii(tree)

    for clade in tree.get_nonterminals():
        aa=[]
        for cladeChild in clade: aa.append(cladeChild.count_terminals())
        collesIndex = collesIndex + abs(aa[0]-aa[1])
        #print aa
    collesIndices.append(collesIndex/(n*n - 3.0*n + 2.0))

    linL=np.zeros(n)
    for clade in tree.find_clades():
        if(clade.branch_length) : 
            #print clade.branch_length, clade.count_terminals()
            linL[clade.count_terminals()]=clade.branch_length
            #print linL 
    lineageLengths.append(linL)    

cI = np.array(collesIndices);
print "Colles Index - mean and std :"
print np.mean(cI), np.std(cI)
lL = np.array(lineageLengths);
#print lL
print "mean of relative lengths of lineages leading to singleton, doubleton, ..., (n-1)-ton ineages and total tree size"
meanlL = np.mean(lL,0)[1:n]
sumMeanlL = sum(meanlL)
print meanlL/sumMeanlL
print sumMeanlL
#print np.std(lL,0)[1:n]
