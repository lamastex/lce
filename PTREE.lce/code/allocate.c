/*                      allocate.c                          */
#include "tree.h"
TREE *tallocate(int leafnumber) {
   TREE *tree;
   tree = (TREE *)getmem(sizeof(TREE));
   tree->leaves = (int *)getmem(leafnumber*sizeof(int));
   tree->multiplicity = (int *)getmem((leafnumber+1)*sizeof(int));
   return tree;
}
void tdispose(TREE *tree) {
    free(tree->multiplicity);
    free(tree->leaves);
    free(tree);
}
void copy(TREE *tr,TREE *tree,int leafnumber) {
    memcpy(tr->leaves,tree->leaves,leafnumber*sizeof(int));
    memcpy(tr->multiplicity,tree->multiplicity,(leafnumber+1)*sizeof(int));
}
void ridnodes(void) {
    int i = 0;
    while(nodelist[i].sibs != NULL) 
        free(nodelist[i++].sibs);
    free(nodelist);
}                
