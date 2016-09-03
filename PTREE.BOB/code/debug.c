/*                     debug.c                 */
#include "tree.h"
void printtree(TREE *tree) {
    int leaf;
    printf("\n\tleaf\tmultiplicity");
    for(leaf = 0;tree->multiplicity[leaf];leaf++)
        printf("\n\t%d\t%d",tree->leaves[leaf],tree->multiplicity[leaf]);
    printf("\n\n");    
}
void debugnodes(NODE *nodelist,int d) {
    int i,j,printblock;
    for (i = 0;i < d;i++)  {
        printf("\n\tNode %d",i);
        printf("\n\tancestor = %d",nodelist[i].ancestor);
#if SCALE        
	printblock = nodelist[i].scale == -1 ? -1 : block[i];	
        printf("\n\tBlock = %d, Scale = %d",printblock,nodelist[i].scale);
#endif        
        printf("\n\tSiblist ");
        if (!nodelist[i].sibnumber) printf("Null");
        else for (j = 0;j < nodelist[i].sibnumber;j++)
            printf("%d  ",nodelist[i].sibs[j]);
        printf("\n");            
    }
    printf("\n\n");
}
/* d = sample size + number of segregating sites + 1 */            
