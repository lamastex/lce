#!/bin/bash

TREES_SAMPLES=100
SIEVE_SIZE=100
THETA_POINTS=100
GROWTH_POINTS=100
SFS_SAMPLES=100
SFS_BURNIN=10
GRPARAMIN=0.0
GRPARAMAX=100.0

./ebc_MT_SPIETA -b 29 -d sfsms.1 -n $SIEVE_SIZE -N 10000000 -q 0 -J 0 -g $GRPARAMIN -G $GRPARAMAX -t 0.001 -T 100.0 -I $SFS_SAMPLES -R $TREES_SAMPLES -s 0 -m SPI30moves -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -z 101 -o out -E MTP

./ebc_MT_SPIETA.0 -b 29 -d sfsms.1 -n $SIEVE_SIZE -N 10000000 -q 0 -J 0 -g $GRPARAMIN -G $GRPARAMAX -t 0.001 -T 100.0 -I $SFS_SAMPLES -R $TREES_SAMPLES -s 0 -m SPI30moves -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -z 101 > out.0 

