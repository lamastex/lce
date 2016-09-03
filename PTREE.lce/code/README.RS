Try to make Bob's code compile cleaner and warning free.  
Currently under:
--------------------------------------------------------------------------
$ cc --version
i686-apple-darwin8-gcc-4.0.1 (GCC) 4.0.1 (Apple Computer, Inc. build 5367)
Copyright (C) 2005 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
--------------------------------------------------------------------------
We get these warnings.  Perhaps easily fixable someday -- rz.
--------------------------------------------------------------------------
$ make
cc -Wall -W -ansi -pedantic    -o ptree allocate.c comb.c debug.c file.c main.c memory.c path.c pattern.c prob.c sample.c setup.c store.c utility.c 
comb.c: In function 'multinomial':
comb.c:35: warning: suggest parentheses around comparison in operand of |
comb.c: In function 'nextper':
comb.c:56: warning: suggest parentheses around comparison in operand of &
comb.c:57: warning: suggest parentheses around comparison in operand of &
comb.c:57: warning: suggest parentheses around comparison in operand of &
comb.c:57: warning: suggest parentheses around comparison in operand of &
comb.c:58: warning: suggest parentheses around comparison in operand of &
comb.c:58: warning: suggest parentheses around comparison in operand of &
comb.c:59: warning: suggest parentheses around comparison in operand of |
comb.c:66: warning: suggest parentheses around comparison in operand of &
comb.c:73: warning: suggest parentheses around comparison in operand of &
comb.c:76: warning: suggest parentheses around comparison in operand of &
comb.c:79: warning: suggest parentheses around comparison in operand of &
comb.c:82: warning: suggest parentheses around comparison in operand of &
comb.c: In function 'nextlex':
comb.c:101: warning: suggest parentheses around comparison in operand of &
comb.c: In function 'nextlexper':
comb.c:119: warning: suggest parentheses around comparison in operand of &
file.c: In function 'strip':
file.c:29: warning: suggest parentheses around comparison in operand of &
file.c:35: warning: suggest parentheses around comparison in operand of &
file.c:41: warning: suggest parentheses around comparison in operand of &
main.c: In function 'main':
main.c:75: warning: suggest parentheses around comparison in operand of |
main.c:129: warning: suggest parentheses around comparison in operand of &
main.c:142: warning: ISO C90 does not support the '%lf' printf format
main.c:145: warning: ISO C90 does not support the '%lf' printf format
main.c:150: warning: control reaches end of non-void function
path.c: In function 'interior':
path.c:16: warning: suggest parentheses around comparison in operand of &
pattern.c: In function 'factor':
pattern.c:44: warning: unused variable 'j'
prob.c: In function 'P':
prob.c:35: warning: suggest parentheses around comparison in operand of |
prob.c:35: warning: suggest parentheses around comparison in operand of &
prob.c:36: warning: suggest parentheses around comparison in operand of |
prob.c:49: warning: suggest parentheses around comparison in operand of &
setup.c: In function 'gettree':
setup.c:40: warning: suggest parentheses around comparison in operand of &
setup.c:131: warning: suggest parentheses around comparison in operand of &
setup.c:150: warning: suggest parentheses around comparison in operand of &
--------------------------------------------------------------------------
However, it runs fine:
--------------------------------------------------------------------------
$ ./ptree specimen.dat 1.0 2.0 3.0

        Data input  =    2
        Data input  =    1
        Data input  =    0
        Data input  =    2
        Data input  =    4
        Data input  =    0
        Data input  =    3
        Data input  =    6
        Data input  =    4
        Data input  =    0
        Data input  =    1
        Data input  =    5
        Data input  =    3
        Data input  =    0
        Data input  =    1
Probability of the unordered tree when theta = 1.000000 is 4.40749e-05

Probability of the unordered tree when theta = 2.000000 is 0.0001479

Probability of the unordered tree when theta = 3.000000 is 0.000167074
--------------------------------------------------------------------------

Try to compile on 64-bit machine and see if the out-of-stack problem doesn't
keep us from analyzing larger trees.  I think, currently for our ms simulations
of long locus ptree only worked fine for samples of size 4.  Try running on
samples of size 5 in 64-bit machines.  FIRST fix warning above.

RS
Tue Feb  3 22:44:08 NZDT 2009
