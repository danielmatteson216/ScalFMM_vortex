#!/bin/bash
#
#if ($#
filename=$1
gnuscript=~/mount/Plafrim1/Dev/src/ScalFMM/scalfmm/Utils/scripts/cmpAlgoMemoryPlot.gnu
FILEBASIC=`cat $filename  | grep basic`
FILETASK=`cat $filename  | grep "\-task"`
FILEBALANCED=`cat $filename  | grep balanced `
FILESECTIONTASK=`cat $filename  | grep section `

outputeps=`basename $filename  txt`"eps"
TITLE=`basename $filename  .txt`
echo "gnuplot"  $outputeps
gnuplot -e "FILEEPS='$outputeps'; TITLE='$TITLE'; FILEBASIC='$FILEBASIC' ; FILETASK='$FILETASK';  FILEBALANCED='$FILEBALANCED';  FILESECTIONTASK='$FILESECTIONTASK' "   $gnuscript

