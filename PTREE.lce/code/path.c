/*                       path.c                    */
#include "tree.h"
/* Return
    0 if label not in nodelabel's path to the root,
    1 if a node in the path,
      Assume the root has a label 0.
*/    

int inpath(int label,int nodelabel) {
    if (label == nodelabel) return 1;
    if (nodelabel == 0) return 0;
    return inpath(label,parent(nodelabel));
}
int interior(TREE *tree,int leaf) {
    int i,test = 0;
    for(i = 0;tree->multiplicity[i] != 0 & test == 0;i++)
        if (tree->leaves[i] != leaf)
            test = inpath(leaf,tree->leaves[i]);
    return test;
}        

int parent(int label) { return nodelist[label].ancestor ; }
