#!/bin/bash
TREES_SAMPLES=500000
TREES_SAMPLES=500000
TREES_SAMPLES=100
TREES_SAMPLES_PER_THETA=1
THETA_POINTS=101
GROWTH_POINTS=101
SFS_SAMPLES=10000#0000
SFS_BURNIN=10000#10000#0
THPARAMIN=1.0
THPARAMAX=20.0
THPARAMIN=0.0001
THPARAMAX=100.0
GRPARAMIN=0.0
GRPARAMAX=100.0
#GRPARAMAX=10.0
THIN_OUT_SFS=1
ALPHA_PRIOR=10.0
NoDATASETS=100
SIEVE_SIZE=1000
SIEVE_SIZE=1
COAL_CONDL_SFS=0
SIEVE_MAX_TRIES=10000000
QUIET=55
QUIET=1
QUIET=0
#DATA=/Users/raazesh/Machome/gtree90/IS_CLADE/MS_Simul_Data/n_30.t_10.r_0.g_0/sfsms.712
#DATA=/Users/raazesh/Machome/gtree90/IS_CLADE/MS_Simul_Data/n_10.t_10.r_0.g_0/sfsms
#DATA=/Users/raazesh/Machome/gtree90/IS_CLADE/MS_Simul_Data/n_90.t_10.r_0.g_0/SFS
BINS=29
SAMPLEID=30
BINS=7
SAMPLEID=08
#DIR=/Users/raazesh/Machome/gtree90/IS_CLADE/MS_Simul_Data/n_$SAMPLEID.t_10.r_0.g_0
DIR=/home/raaz/Machome/gtree90/IS_CLADE/MS_Simul_Data/n_$SAMPLEID.t_10.r_0.g_0
#DATA=$DIR/sfsms
DATA=$DIR/SFS_n_$SAMPLEID
COUNTDATA=$DIR/CountSFS_n_$SAMPLEID
OUTFILE=./OUT_n_$SAMPLEID/PrSi
MOVESFILE=./SPI.moves/SPI_n_$SAMPLEID.moves

#mkdir $DIR/SFSsOfSPI_n_$SAMPLEID
#rm $OUTFILE.naninfs
SIMSIZE=2
REPS=1
SEED=3787795$REPS
while [  $REPS -lt $SIMSIZE ]; do
#	echo _____________________________
#	echo The replicate is $REPS
#	echo _____________________________
	
#DATA=sfsms.712
#echo -n $REPS 

#SIEVE_SIZE=`head -n $REPS $COUNTDATA | tail -n 1`
#echo count is $SIEVE_SIZE
#echo -n $SIEVE_SIZE


COAL_CONDL_SFS=0
#time ./ebc_MT_SPIETA -C 0 -H $ALPHA_PRIOR -b $BINS -d $DATA.$REPS -D $NoDATASETS -Y $REPS -n $SIEVE_SIZE -N $SIEVE_MAX_TRIES -q $QUIET -J 10 -g $GRPARAMIN -G $GRPARAMAX -t $THPARAMIN -T $THPARAMAX -I $SFS_SAMPLES -R $TREES_SAMPLES -s $TREES_SAMPLES_PER_THETA -m $MOVESFILE -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -K $THIN_OUT_SFS -h 100000 -i 500000 -o $OUTFILE.$REPS -E $OUTFILE.MTP_$TREES_SAMPLES.n_$SAMPLEID  -z $SEED -c $COAL_CONDL_SFS #| grep -e '(' | sed 's/:[0-9]*.[0-9]*//g' | sed 's/;//g' > C$COAL_CONDL_SFS.trees.$SAMPLEID
#time ./ebc_MT_SPIETA -C 0 -H $ALPHA_PRIOR -b $BINS -d $DATA.$REPS -D $NoDATASETS -Y $REPS -n $SIEVE_SIZE -N $SIEVE_MAX_TRIES -q $QUIET -J 10 -g $GRPARAMIN -G $GRPARAMAX -t $THPARAMIN -T $THPARAMAX -I $SFS_SAMPLES -R $TREES_SAMPLES -s $TREES_SAMPLES_PER_THETA -m $MOVESFILE -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -K $THIN_OUT_SFS -h 100000 -i 500000 -o $OUTFILE.$REPS -E $OUTFILE.MTP_$TREES_SAMPLES.n_$SAMPLEID  -z 22023$REPS -c $COAL_CONDL_SFS #| grep -v '#' | tee $DIR/SFSsOfSPI_n_$SAMPLEID/SFSsOfSPI_n_$SAMPLEID.$REPS

#cat C$COAL_CONDL_SFS.trees.$SAMPLEID | sed 's/[0-9]*//g' | sed 's/,//g' > C$COAL_CONDL_SFS.treesT.$SAMPLEID

./ebc_MT_SPIETA -C 0 -H $ALPHA_PRIOR -b $BINS -d $DATA -D $NoDATASETS -Y $REPS -n $SIEVE_SIZE -N $SIEVE_MAX_TRIES -q $QUIET -J 10 -g $GRPARAMIN -G $GRPARAMAX -t $THPARAMIN -T $THPARAMAX -I $SFS_SAMPLES -R $TREES_SAMPLES -s $TREES_SAMPLES_PER_THETA -m $MOVESFILE -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -K $THIN_OUT_SFS -h 100000 -i 500000 -o $OUTFILE -E $OUTFILE.MTP_$TREES_SAMPLES.n_$SAMPLEID  -z $SEED -c $COAL_CONDL_SFS  | tee ./OUT_n_$SAMPLEID/SimRes_SFS_$SIEVE_SIZE
#./ebc_MT_SPIETA -C 0 -H $ALPHA_PRIOR -b $BINS -d $DATA.$REPS -D $NoDATASETS -Y $REPS -n $SIEVE_SIZE -N $SIEVE_MAX_TRIES -q $QUIET -J 10 -g $GRPARAMIN -G $GRPARAMAX -t $THPARAMIN -T $THPARAMAX -I $SFS_SAMPLES -R $TREES_SAMPLES -s $TREES_SAMPLES_PER_THETA -m $MOVESFILE -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -K $THIN_OUT_SFS -h 100000 -i 500000 -o $OUTFILE.$REPS -E $OUTFILE.MTP_$TREES_SAMPLES.n_$SAMPLEID  -z $SEED -c $COAL_CONDL_SFS  | tee ./OUT_n_$SAMPLEID/SimRes_SFS_$SIEVE_SIZE
#| grep -e '(' | sed 's/:[0-9]*.[0-9]*//g' | sed 's/;//g' > C$COAL_CONDL_SFS.trees.$SAMPLEID

#cat C$COAL_CONDL_SFS.trees.$SAMPLEID | sed 's/[0-9]*//g' | sed 's/,//g' > C$COAL_CONDL_SFS.treesT.$SAMPLEID
#time ./ebc_MT_SPIETA -C 0 -H $ALPHA_PRIOR -b $BINS -d $DATA.$REPS -D $NoDATASETS -Y $REPS -n $SIEVE_SIZE -N $SIEVE_MAX_TRIES -q $QUIET -J 10 -g $GRPARAMIN -G $GRPARAMAX -t $THPARAMIN -T $THPARAMAX -I $SFS_SAMPLES -R $TREES_SAMPLES -s $TREES_SAMPLES_PER_THETA -m $MOVESFILE -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -K $THIN_OUT_SFS -h 100000 -i 500000 -o $OUTFILE.$REPS -E $OUTFILE.MTP_$TREES_SAMPLES.n_$SAMPLEID  -z 22023$REPS -c $COAL_CONDL_SFS #| grep -v '#' | tee $DIR/SFSsOfSPI_n_$SAMPLEID/SFSsOfSPI_n_$SAMPLEID.$REPS
#COAL_CONDL_SFS=0
#echo std coalescent
#./ebc_MT_SPIETA -C 2 -H $ALPHA_PRIOR -b $BINS -d $DATA.$REPS -D $NoDATASETS -Y $REPS -n $SIEVE_SIZE -N $SIEVE_MAX_TRIES -q $QUIET -J 10 -g $GRPARAMIN -G $GRPARAMAX -t $THPARAMIN -T $THPARAMAX -I $SFS_SAMPLES -R $TREES_SAMPLES -s $TREES_SAMPLES_PER_THETA -m $MOVESFILE -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -K $THIN_OUT_SFS -h 100000 -i 500000 -o $OUTFILE.$REPS -E $OUTFILE.MTP_$TREES_SAMPLES.n_$SAMPLEID  -z 22023$REPS -c $COAL_CONDL_SFS #| grep -v '#' | tee $DIR/SFSsOfSPI_n_$SAMPLEID/SFSsOfSPI_n_$SAMPLEID.$REPS
#echo -n $REPS >> $OUTFILE.naninfs
#echo -n "   " >> $OUTFILE.naninfs
#cat $OUTFILE.$REPS | grep 'nan' >> $OUTFILE.naninfs
#cat $OUTFILE.$REPS | grep 'inf' >> $OUTFILE.naninfs
#echo " _" >> $OUTFILE.naninfs
#cat $DIR/SFSsOfSPI_n_$SAMPLEID/SFSsOfSPI_n_$SAMPLEID.$REPS | wc -l | sed 's/ //g' >> $DIR/SFSsOfSPI_n_$SAMPLEID/CountSFSHeurSrch_n_$SAMPLEID

#./ebc_MT_SPIETA -H $ALPHA_PRIOR -b $BINS -d $DATA.$REPS -n $SIEVE_SIZE -N 1000000 -q 0 -J 10 -g $GRPARAMIN -G $GRPARAMAX -t 0.001 -T 100.0 -I $SFS_SAMPLES -R $TREES_SAMPLES -s 0 -m SPI10moves -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -z 10231 -K $THIN_OUT_SFS -h 100000 -i 500000 -o out -E MTP 
let REPS=REPS+1
done
#time ./ebc_MT_SPIETA -H $ALPHA_PRIOR -b 29 -d $DATA -n $SIEVE_SIZE -N 1000000 -q 1 -J 10 -g $GRPARAMIN -G $GRPARAMAX -t 0.001 -T 100.0 -I $SFS_SAMPLES -R $TREES_SAMPLES -s 0 -m SPIEXT30moves -M 3 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -z 10231 -K $THIN_OUT_SFS -h 100000 -i 500000 -o out -E MTP 
#time ./ebc_MT_SPIETA.0 -b 29 -d sfsms.1 -n $SIEVE_SIZE -N 10000000 -q 0 -J 0 -g $GRPARAMIN -G $GRPARAMAX -t 0.001 -T 100.0 -I $SFS_SAMPLES -R $TREES_SAMPLES -s 0 -m SPI30moves -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -z 101 -o out.0 -E MTP

#./ebc_MT_SPIETA.0 -b 29 -d sfsms.1 -n $SIEVE_SIZE -N 10000000 -q 0 -J 0 -g $GRPARAMIN -G $GRPARAMAX -t 0.001 -T 100.0 -I $SFS_SAMPLES -R $TREES_SAMPLES -s 0 -m SPI30moves -M 2 -P $THETA_POINTS -Q $GROWTH_POINTS -A 1 -W $SFS_BURNIN -z 101 > out.0 

