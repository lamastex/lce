/*                    pattern.h              */
/* The distinct binary site patterns of the rows of incidence[]
    and their multiplicity.Use the skeleton tree. 
*/    
typedef struct {
    int *pattern;
    int multiplicity; /* Null terminated. */ 
} PATTERN;

extern PATTERN* getpattern(int*,int,int);
extern int inpath(int,int);
extern int nextlexper(int*,int,int*,int);
extern int permtest(int*,int*,PATTERN*,int,int);

