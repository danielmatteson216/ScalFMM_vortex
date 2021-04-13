#!/bin/bash -v
#
#PBS -N MG_water
# 
# preciser le temps en heures, minutes, secondes
#PBS -l walltime=71:59:00
# nombre de noeuds et de coeurs
#PBS -l nodes=1:ppn=8
#PBS -l naccesspolicy=singlejob
#
# PlaFRIM environment
source $HOME/.bashrc
#
FMMPER_GEN_EXE="Tests/Release/DirectAlgorithm"

FILE="../Data/forceNacl_128_dlpolyPer.bin"

SCALFMM_WD=~/Dev/src/ScalFMM/scalfmmT/BuildGCC

usage ()
{
  echo  "--------------------------------------------------------------------"
  echo usage: $0 DEGRE
  echo
  echo DEGRE in the approximation of the expansion or the inrepolation GRID
  echo    used to store the output in FMM_PER_${DEGRE}.out
echo  "--------------------------------------------------------------------"
exit
}
#
# [ $1 ] || usage
# dd=$1
#
# RUN 1 influence de la taille de la boite sur la precision des calculs
#    Regarde size, Nx ny nz energie total, dipole et le temps de calcul
#
PER_SIZE="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15"
PER_SIZE="0 1 "
#
pwd
cd ${SCALFMM_WD}

FMMPER_EXE=${FMMPER_GEN_EXE}
rm -f  FMM_PER_DIRECT.out
echo  "# PerSize DOMAIN_SIZE  ENERGY"> FMM_PER_DIRECT.out

for l in $PER_SIZE
do
    OUTPUT=OUTPUT-DIRECT-${l}.out
    filenameOUT=Nacl-128-per=${l}.bin

    echo "Running  per = " ${l}
    echo  "${FMMPER_EXE} -bin -f $FILE  -h 4 -sh 2 -per $l -direct  >  $OUTPUT "
#
    ${FMMPER_EXE} -bin -fin $FILE  -h 4 -sh S2  -per ${l} -fin $FILE -bin  -fout $filenameOUT   >  $OUTPUT
#
    DOMAINSIZE=`grep "Simulated box:"  $OUTPUT | awk '{print $3}'`
    ENERGYDIRECT=`grep  "Energy DIRECT="  $OUTPUT |awk '{ print $3 }'`
    TIMEDIRECT=`grep "@DIRECT"  $OUTPUT |awk '{ print $5 }'`
#    TIMEFMM=`grep "@FMM"  $OUTPUT |awk '{ print $5 }'`
#    ERROR=`grep "Energy EWALD - Energy DIRECT|/directEnergy="  $OUTPUT |awk '{ print $7}'`
    echo " " $l   $DOMAINSIZE "   "  $TIMEDIRECT " "   $ENERGYDIRECT " "  >>  FMM_PER_DIRECT.out
done



