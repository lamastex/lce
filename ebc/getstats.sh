#!/bin/bash

#define the jobs where everything lies
#/Users/kt234/work/current_project/
DIR=/Users/raazesh/Machome/gtree90/IS_CLADE/MS_Simul_Data/n_10.t_10.r_0.g_0
SIMSIZE=1001
# REPS needs to be replaced by your $SGE_TASK_ID thingy...
REPS=1
#cat $DIR/hapms.$REPS | msstats | tail 
while [  $REPS -lt $SIMSIZE ]; do
#echo The replicate is $REPS
#check to make sure we don't over-write output
#if [ ! -e $DIR/srfgt/summaries ]
#then
cat $DIR/hapms.$REPS | msstats | tail +2 | colrm 1 96 | colrm 14
#fi
let REPS=REPS+1 
done



