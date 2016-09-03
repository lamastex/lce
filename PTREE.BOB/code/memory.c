/*                  memory.c                         */
#include <stdio.h>
#include <stdlib.h>
void *getmem(int bytes) {
    void *tmp;
    if(bytes <= 0) return NULL;
    if((tmp = malloc(bytes)) == NULL) {
        printf("\n\tMemory allocation error\n");
        printf("\tCannot allocate %d bytes\n",bytes);
        exit(1);
    }
    return tmp;
}
void *getmemc(int nobj,int siz) {
void *tmp;
    if(siz <= 0) return NULL; 
     if((tmp = calloc(nobj,siz)) == NULL) {
        printf("\nMemory allocation error\n");
        printf("\tCannot allocate %d bytes\n",nobj*siz);        
        exit(1);
    }
    return tmp;
}   
