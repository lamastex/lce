This code of Bob's doesn't compile under:


i686-apple-darwin8-gcc-4.0.1 (GCC) 4.0.1 (Apple Computer, Inc. build 5367)
Copyright (C) 2005 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

gives error:

cc -Wall -W -ansi -pedantic    -o simulate SIMULATE.c 
SIMULATE.c:31: error: conflicting types for 'random'
/usr/include/stdlib.h:222: error: previous declaration of 'random' was here
SIMULATE.c: In function 'randtree':
SIMULATE.c:64: warning: implicit declaration of function 'memcpy'
SIMULATE.c:64: warning: incompatible implicit declaration of built-in function 'memcpy'
SIMULATE.c: In function 'main':
SIMULATE.c:88: warning: suggest parentheses around comparison in operand of |
make: *** [all] Error 1


SIMULATE.EXE file came with Bob's original package and it compiles under windows, I guess.

Also, I am not sure what SIMULATE does.  Don't remember ever using it.

Someday, look at the source SIMULATE.C and start commenting it heavily and try to make it compile cleanly.

RS
Tue Feb  3 22:42:25 NZDT 2009

