#!/bin/bash
SAMPLES=10
THETA=10
REPEATS=1
REPS=1
SIMSIZE=1001
cp seedms.$REPS seedms 
while [  $REPS -lt $SIMSIZE ]; do
echo The replicate is $REPS
#sleep 1
ms $SAMPLES $REPEATS -t $THETA -T > outms.$REPS 
#sleep 1
cat outms.$REPS | tail +4 | grep '(' > genms.$REPS
#sleep 1
cat outms.$REPS | grep -v '(' > hapms.$REPS 
#sleep 1
cat hapms.$REPS | msgt > hapgt.$REPS
#sleep 1
cat hapms.$REPS | msfreq > sfsms.$REPS
#sleep 1
let REPS=REPS+1 
cp seedms seedms.$REPS
done
rm seedms
date > finished
