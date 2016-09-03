#!/bin/bash
SAMPLES=90
THETA=10
REPS=1000
cp n_$SAMPLES.t_$THETA.r_$REPS.seedms seedms 
sleep 2
ms $SAMPLES $REPS -t $THETA -T > n_$SAMPLES.t_$THETA.r_$REPS.outms 
sleep 2
cat n_$SAMPLES.t_$THETA.r_$REPS.outms | tail +4 | grep '(' > n_$SAMPLES.t_$THETA.r_$REPS.genms
sleep 2
cat n_$SAMPLES.t_$THETA.r_$REPS.outms | grep -v '(' > n_$SAMPLES.t_$THETA.r_$REPS.hapms 
sleep 2
cat n_$SAMPLES.t_$THETA.r_$REPS.hapms | msgt > n_$SAMPLES.t_$THETA.r_$REPS.hapgt 
sleep 2
cat n_$SAMPLES.t_$THETA.r_$REPS.hapms | msfreq > n_$SAMPLES.t_$THETA.r_$REPS.sfsms
sleep 2
rm seedms
date > finished
