#!/usr/bin/env bash
#SBATCH --job-name=OpenMPAlgo_perf
#SBATCH  --mail-user=olivier.coulaud@inria.fr
#SBATCH -m a
#
#SBATCH -p court
#SBATCH --time=60:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#

if test $# -eq 1
then
     echo "Good number of parameter"
else
    echo "Wrong number of parameter "
    echo  "The script needs 1 arguments coreponding to"
    echo  "  Arg1 the type of the interpolation method (0 Cheb, 1 Uniform)"
    exit
fi
data_dir=/projets/scalfmm/data/openmpImprovements ;
#project_dir=$HOME/Dev/src/ScalFMM/scalfmm ;

#
# PlaFRIM environment
#
source $HOME/Config/bashrc.bash
w
#
if test $1 -eq 0
then
  EXEC="Utils/Release/ChebyshevInterpolationCmpAlgo"
  Interp="Cheb"
else
  EXEC="Utils/Release/UniformInterpolationCmpAlgo"
  Interp="Unif"
fi
FILEPERF=RES-${Interp}-openmp

OPTION="-show-info -show-compile"
#
DATE=`date +%y-%m-%d-%H%M`
echo $DATE
#
#
#   INTEL 2015
#
# module add compiler/gcc/5.1.0 compiler/intel/64/2015.3.187
# cd $project_dir/BuildIntel2015
#
#   INTEL 2016
#
#module add compiler/gcc/5.1.0 compiler/intel/64/2016_beta
#cd $project_dir/BuildIntel2016
#
#   GCC  + MKL
#
module add compiler/gcc/7.2.0   intel/mkl/64/2018/update1
#

module li
#
#
PER_SIZE=`cat /proc/cpuinfo |grep processor |wc -l`
PER_SIZE="4 8 12 16 20 24"
HOST=`hostname`
NBTEST=3  # 0..NBTEST
DISTRIB=(unitcube unitsphere prolate plummer)
DEPTH=(6           7           9            8)    # 1M
#DEPTH=(7           8           10            9)   # 10M

if test ${#DISTRIB[*]}  -ne ${#DEPTH[*]}
then
    echo " DEPTH and DISTRIB arrays don't have se same size"
    exit
fi

ALGO=(task newtask sectiontask omp4)
#
pwd
export OMP_PROC_BIND=true
NUM=`git rev-list HEAD --count`
echo   $DISTRIB
mkdir  ${DATE}
for dist in ${!DISTRIB[*]}    ;
do
    echo "Distribution:  " $dist ${DISTRIB[$dist]}   "  depth octree "  ${DEPTH[$dist]}
    FILE=${data_dir}/${DISTRIB[$dist]}-1M-ref.bfma      #1M
#    FILE=${data_dir}/${DISTRIB[$dist]}-10M.bfma         #10M
    FILEGNU=${DISTRIB[$dist]}-${Interp}.txt
    mkdir   ${DATE}/${DISTRIB[$dist]} 
    touch  $FILEGNU
    #
    #    Loop on algorithm 
    for a in ${!ALGO[*]}    ;
    do 
echo  "# Core GlabalTIME P2PTIME M2LTIME  MALLOC   FREE MEMUse ENERGY PotentialError ForceError"> $FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out
	#  One computation
	for l in $PER_SIZE ;
	do
	    OUTPUT=${FILEPERF}-${DISTRIB[$dist]}-${ALGO[$a]}-${l}.out
	    echo "ScalFMM git version: " $NUM > $OUTPUT
	    echo "Running  on " ${l} " threads with algo "   ${ALGO[$a]}
	    $EXEC ${OPTION} -f $FILE  -depth  ${DEPTH[$dist]}  -subdepth 4  -t $l  -cmp  -algo ${ALGO[$a]}>>  $OUTPUT  #1M
#	    $EXEC ${OPTION} -f $FILE  -depth  ${DEPTH[$dist]}  -subdepth 4  -t $l  -algo ${ALGO[$a]}>>  $OUTPUT      #10M
	    #
	    TIME=`grep "@Algorithm"  $OUTPUT | awk '{print $4}'`
	    TIMEP2P=`grep "P2P"  $OUTPUT | grep "second" | awk '{print $4}'`
	    TIMEM2L=`grep "M2L"  $OUTPUT | grep "second" | awk '{print $2}'`
	    Energy=`grep "Energy"  $OUTPUT | awk '{print $3}'`
	    P1=`grep " Potential "  $OUTPUT | awk '{print $7}'`
	    P2=`grep " F "  $OUTPUT | awk '{print $7}'`
	    echo " " $l "   " $TIME "  " $TIMEP2P "  " $TIMEM2L "  " $Energy "     " $P1  "  " $P2
	    echo " " $l "   " $TIME "  " $TIMEP2P "  " $TIMEM2L "  " $Energy "     " $P1  "  " $P2  >>  $FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out
	    mv  $OUTPUT   ${DATE}/${DISTRIB[$dist]}
	done
	echo  ${DISTRIB[$dist]}/$FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out  >> $FILEGNU
	mv $FILEPERF-${DISTRIB[$dist]}-${ALGO[$a]}.out  ${DATE}/${DISTRIB[$dist]}
    echo "----------------------------------------------------------"
    done
    echo "---------------------------- END ALGO ---------------------------------------"
   mv   $FILEGNU  ${DATE}
done
echo " ----------------------   FINISH -------------------------------"
