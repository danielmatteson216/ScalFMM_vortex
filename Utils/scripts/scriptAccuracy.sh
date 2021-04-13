#!/bin/bash
#
EXEC=Utils/Release/testAccuracyChebFMM
source $HOME/.bashrc
FILE=/Users/coulaud/Dev/src/ScalFMM/scalfmmT/Data/noDist/CEA/2000REFParticules.bfma
#
rm output*
for o in {2..13}
do
${EXEC} --datehost --flags -show-params -order $o -f $FILE -depth 5 -subdepth 3 -t 1 > output-accuracy-$o
done



