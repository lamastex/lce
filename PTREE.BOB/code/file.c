/*                  file.c                          */
#include "tree.h"

void mismatch(void) {
    printf("\n\tComments must be enclosed by  /*    */\n");
    exit(1);
}
/* Strip C comments from a file, and read in the non-comments
    into a buffer returned by strip. *plength is the length of
    the returned (null-terminated) string.
*/    
char *strip(char *fname) {
    FILE *fp;
    int i = 0,lastc = 0,c = 0,comment = 0,filelength = 0;
    char *buffer;
    if((fp = fopen(fname,"r")) == NULL) {
        printf("\n\tNo data input file %s\n",fname);
        exit(1);
    }
    while(c != EOF) {
        c = getc(fp);
        filelength++;
    }
    buffer = (char *)getmem(filelength);
    rewind(fp); 
    c = 0;
    while(c != EOF) {
        c = getc(fp);
        if(c == '*' & lastc == '/') {
            if(comment)
              mismatch();
            else 
                comment = 1;
        }
        if(c == '/' & lastc == '*') {
            if(!comment) 
               mismatch();
            else 
                comment = 0;
        }
        if(!comment & c != '/') {
            if((isdigit(c)|isspace(c)))
                buffer[i++] = c;
            else if(c != EOF) 
                mismatch();
        }
        lastc = c;
    }
    if(fclose(fp) == EOF) {
        printf("\n\tFile close error\n");
        exit(1);
    }
    buffer[i] = 0;
    return buffer;
}    
    
/* Get stream from a null terminated buffer */

int *getstream(char *buf) {
    int i = 0,in,test = 1;
    int *stream;
    char *tbuf;
    tbuf = buf;
    while(isspace(*tbuf)) tbuf++;
    stream = (int *)getmem((strlen(buf)+2)*sizeof(int));
    while ((test = sscanf(tbuf,"%d",&in)) != EOF) {
            printf("\n\tData input  = \t %d",in);
        if(!test) {
            printf("\n\tData input error.\n");
            exit(1);
        }
        else  {
            stream[i++] = in;
			while(!isspace(*tbuf))
				tbuf++;
            while(isspace(*tbuf))
                tbuf++;
        }
    }
        stream[i] = -1;
        free(buf);       
        if(!i) {
            printf("\n\tNo input data\n");
            exit(1);
        }
        return stream;
}
