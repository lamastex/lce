#!/bin/bash
SAMPLES=10
THETA=100
SITES=100000
REPEATS=10000
RHO=0
GROWTH=0
ms $SAMPLES $REPEATS -t 100 -G $GROWTH  -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
GROWTH=10
ms $SAMPLES $REPEATS -t 100 -G $GROWTH  -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
GROWTH=100
ms $SAMPLES $REPEATS -t 100 -G $GROWTH  -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
RHO=10
GROWTH=0
ms $SAMPLES $REPEATS -t 100 -G $GROWTH -r $RHO $SITES -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
GROWTH=10
ms $SAMPLES $REPEATS -t 100 -G $GROWTH -r $RHO $SITES -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
GROWTH=100
ms $SAMPLES $REPEATS -t 100 -G $GROWTH -r $RHO $SITES -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
RHO=100
GROWTH=0
ms $SAMPLES $REPEATS -t 100 -G $GROWTH -r $RHO $SITES -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
GROWTH=10
ms $SAMPLES $REPEATS -t 100 -G $GROWTH -r $RHO $SITES -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
GROWTH=100
ms $SAMPLES $REPEATS -t 100 -G $GROWTH -r $RHO $SITES -T | tee outms_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.g_$GROWTH.r_$RHO.m_$SITES.reps_$REPEATS
#ms $SAMPLES $REPEATS -t $THETA -G $GROWTH  | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.r_$RHO.g_$GROWTH
#ms $SAMPLES $REPEATS -t $THETA -r $RHO $SITES  | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.r_$RHO.g_$GROWTH
#ms $SAMPLES $REPEATS -t $THETA -G $GROWTH -r $RHO $SITES  | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_$THETA.r_$RHO.g_$GROWTH
#############
#ms_multitheta $SAMPLES $REPEATS -t 3 1 10 100 -G $GROWTH  | grep -v '(' | msfreq >> SFS_n_$SAMPLES.t_1_10_100.r_$RHO.g_$GROWTH

