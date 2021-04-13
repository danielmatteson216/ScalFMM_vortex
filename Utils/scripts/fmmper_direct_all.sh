#!/bin/bash
#
#PBS -N 
# 
# preciser le temps en heures, minutes, secondes
#PBS -l walltime=71:59:00
# nombre de noeuds et de coeurs
#PBS -l nodes=1:ppn=8
#PBS -l naccesspolicy=singlejob
#
# PlaFRIM environment
w
SCALFMM_WD=~/Dev/src/ScalFMM/scalfmmT/BuildGCC
#source $HOME/.bashrc
#
FMMPER_GEN_EXE="Tests/Release/DirectAlgorithm"

cd ${SCALFMM_WD}
#
DEGRE="2 4 6 8 10 12 14"
#
pwd
for dd in $DEGRE
do
  bash ../Utils/fmmper_direct.sh ${dd}
done


