#!/usr/bin/env bash
#SBATCH --job-name=OpenMPAlgo_perf
#SBATCH  --mail-user=olivier.coulaud@inria.fr
#SBATCH -m a
#
#SBATCH -p court
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#
w
#
# Task algorithm Improvement
#
METH='final'
data_dir=/projets/scalfmm/data/openmpImprovements ;
project_dir=$HOME/Dev/src/ScalFMM/scalfmm ;

#
# PlaFRIM environment
#
source $HOME/Config/bashrc.bash
#
EXEC="Utils/Release/ChebyshevInterpolationCmpAlgo"
FILEPERF="RES-Cheb-openmp"
OPTION="-show-info -show-compile"
#
DATE=`date +%y-%m-%d-%H%M`
echo $DATE
#
#
#   INTEL
#
#module add compiler/gcc/5.1.0 compiler/intel/64/2016_beta
#cd $project_dir/BuildIntel2016
#COMP=INTEL
#export KMP_AFFINITY=scatter
#
#   GCC
#
module add compiler/gcc/5.1.0 intel/mkl/64/11.2/2015.3.187  
cd $project_dir/BuildGCC51
COMP=GCC
#
# For eztrace
#
module add mpi/openmpi/gcc/1.8.4 trace/eztrace/1.0
module li
#
#  Compilation 
#
make ChebyshevInterpolationCmpAlgo 
#
#
PER_SIZE=`cat /proc/cpuinfo |grep processor |wc -l`
PER_SIZE="4 8 12 16 20 24"
HOST=`hostname`
NBTEST=3  # 0..NBTEST
DISTRIB=(unitcube unitsphere prolate plummer)
DEPTH=(6           7           9            8)
ALGO=(basic balanced task  tasknew sectiontask sectiontasknew)
#ALGO=( task )
#
pwd
export OMP_PROC_BIND=true
NUM=`git rev-list HEAD --count`
echo   $DISTRIB
REP=${DATE}-${COMP}-${METH}
mkdir  ${REP}
for dist in `seq 0 $NBTEST`;
do
    echo "Distribution:  " $dist ${DISTRIB[$dist]}   "  depth octree "  ${DEPTH[$dist]}
    FILE=${data_dir}/${DISTRIB[$dist]}-1M-ref.bfma
    FILEGNU=${DISTRIB[$dist]}.txt
    mkdir   ${REP}/${DISTRIB[$dist]} 
    touch  $FILEGNU
    #
    #    Loop on algorithm 
    for a in `seq 0 5`;
    do 
echo  "# Core GlabalTIME P2PTIME M2LTIME  MALLOC   FREE MEMUse ENERGY PotentialError ForceError"> $FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out
	#  One computation
	for l in $PER_SIZE ;
	do
	    OUTPUT=${REP}/${DISTRIB[$dist]}/${FILEPERF}-${DISTRIB[$dist]}-${ALGO[$a]}-${l}.out
	    echo "ScalFMM git version: " $NUM > $OUTPUT
	    echo "Running  on " ${l} " threads with algo "   ${ALGO[$a]}
	    echo $EXEC ${OPTION} -f $FILE  -depth  ${DEPTH[$dist]}  -subdepth 4  -t $l  -cmp  -algo ${ALGO[$a]}>>  $OUTPUT
	    $EXEC ${OPTION} -f $FILE  -depth  ${DEPTH[$dist]}  -subdepth 4  -t $l  -cmp  -algo ${ALGO[$a]}>>  $OUTPUT
	    #
	    eztrace_stats /tmp/coulaud_eztrace_log_rank_1 > TRACE
	    TIME=`grep "@Algorithm"  $OUTPUT | awk '{print $4}'`
	    TIMEP2P=`grep "P2P and L2P"  $OUTPUT | grep "second" | awk '{print $4}'`
	    TIMEM2L=`grep "M2L"  $OUTPUT | grep "second" | awk '{print $2}'`
	    Energy=`grep "Energy"  $OUTPUT | awk '{print $3}'`
	    P1=`grep " Potential "  $OUTPUT | awk '{print $7}'`
	    P2=`grep " F "  $OUTPUT | awk '{print $7}'`
	    MALLOC=`grep Total  TRACE | grep malloc | awk '{print $2}'`
	    FREE=`grep Total  TRACE | grep malloc | awk '{print $7}'`
	    MEM=`grep Total  TRACE | grep malloc | awk '{print $15}'`
	    echo " " $l "   " $TIME "  " $TIMEP2P "  " $TIMEM2L "  " $MALLOC "  " $FREE "  " $MEM "  " $Energy "     " $P1  "  " $P2
	    echo " " $l "   " $TIME "  " $TIMEP2P "  " $TIMEM2L "  " $MALLOC "  " $FREE "  " $MEM "  "$Energy "     " $P1  "  " $P2  >>  $FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out
	    mv /tmp/coulaud_eztrace_log_rank_1   ${REP}/${DISTRIB[$dist]}/trace-${DISTRIB[$dist]}-${ALGO[$a]}-${l}

	    rm TRACE
	done
	echo  ${DISTRIB[$dist]}/$FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out  >> $FILEGNU
	mv    $FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out  ${REP}/${DISTRIB[$dist]}
    echo "----------------------------------------------------------"
    done
    echo "---------------------------- END ALGO $a ---------------------------------------"
   mv   $FILEGNU  ${REP}
done
echo " ----------------------   FINISH -------------------------------"
