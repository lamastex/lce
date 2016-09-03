/*                    utility.c                       */
#include "tree.h"
int getnsample(TREE *tree) {    /* Sample number. */
    int count = 0;
    register int i;
    for(i = 0;tree->multiplicity[i];i++)
        count += tree->multiplicity[i];
    return count;
}    
int getnleaves(TREE *tree) {    /* Number of distinct sequences. */
    register int i;
    for(i = 0;tree->multiplicity[i];i++) ;
    return i;
}
int getnsites(NODE *nodelist) { /* Number of segregating sites. */
    register int i;
    for(i = 1;nodelist[i].sibs != NULL;i++);
    return --i;
}
int agecheck(TREE *tree) {
    int i,last,next;
    for(i = 0;tree->multiplicity[i];i++) {
        next = last = tree->leaves[i];
        while(next) {
            next = parent(last);
            if(next > last) 
                return 0;
            last = next;    
        }
    }
    return 1;
}
