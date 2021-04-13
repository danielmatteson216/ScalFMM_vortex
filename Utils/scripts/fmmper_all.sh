#!/bin/bash
#
#PBS -N MG_water
# 
# preciser le temps en heures, minutes, secondes
#PBS -l walltime=01:59:00
# nombre de noeuds et de coeurs
#PBS -l nodes=1:ppn=12
#PBS -l naccesspolicy=singlejob
#
# PlaFRIM environment
w
source $HOME/.bashrc
#
FMMPER_GEN_EXE="Tests/Release/testSphericalEwalAlgorithm"
FILE="../Data/forceNacl_2000_dlpolyPer.txt"
FILE="../Data/forceNacl_128_dlpolyPer.txt"
FILE="../Data/forceNacl_2000_dlpolyPer.bin"
cd ~/Dev/src/ScalFMM/scalfmmT/BuildGCC
#
# RUN 1 influence de la taille de la boite sur la precision des calculs
#    Regarde size, Nx ny nz energie total, dipole et le temps de calcul
#
DEGRE="3 5 7 9 11 13"
PER_SIZE="-1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18"
#
pwd
for dd in $DEGRE
do
 FMMPER_EXE=${FMMPER_GEN_EXE}_${dd}
rm -f FMM_PER_${DEGRE}.out
echo  "# PerSize DOMAIN_SIZE  ENERGY"> FMM_PER_${dd}.out

for l in $PER_SIZE
do
    OUTPUT=OUTPUT-${dd}-${l}.out
    echo "Running  per = " ${l}
    echo  "$FMMPER_EXE -f $FILE  -h 4 -sh 2 -per $l  >  $OUTPUT "
#
    $FMMPER_EXE -f $FILE  -h 4 -sh 2 -per $l  >  $OUTPUT
#
    DOMAINSIZE=`grep "Simulated box:"  $OUTPUT | awk '{print $3}'`
    ENERGY=`grep "Energy FMM="  $OUTPUT |awk '{ print $4 }'`
    ENERGYE=`grep " Energy EWALD="  $OUTPUT |awk '{ print $4 }'`
    ERROR=`grep "|Energy EWALD -Energy FMM|/Energy EWALD= "  $OUTPUT |awk '{ print $7}'`
    echo " " $l   $DOMAINSIZE "   " $ENERGY " "   $ENERGYE " " $ERROR " " >>  FMM_PER_${dd}.out
done
done


