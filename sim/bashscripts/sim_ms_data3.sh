#!/bin/bash
SAMPLES=10
THETA=10
REPEATS=1
REPS=1
SIMSIZE=1001
cp n_$SAMPLES.t_$THETA.seedms.$REPS seedms 
while [  $REPS -lt $SIMSIZE ]; do
echo The replicate is $REPS
#sleep 1
ms $SAMPLES $REPEATS -t $THETA -T > n_$SAMPLES.t_$THETA.outms.$REPS 
#sleep 1
cat n_$SAMPLES.t_$THETA.outms.$REPS | tail +4 | grep '(' > n_$SAMPLES.t_$THETA.genms.$REPS
#sleep 1
cat n_$SAMPLES.t_$THETA.outms.$REPS | grep -v '(' > n_$SAMPLES.t_$THETA.hapms.$REPS 
#sleep 1
cat n_$SAMPLES.t_$THETA.hapms.$REPS | msgt > n_$SAMPLES.t_$THETA.hapgt.$REPS
#sleep 1
cat n_$SAMPLES.t_$THETA.hapms.$REPS | msfreq > n_$SAMPLES.t_$THETA.sfsms.$REPS
#sleep 1
let REPS=REPS+1 
cp seedms n_$SAMPLES.t_$THETA.seedms.$REPS
done
rm seedms
date > finished
