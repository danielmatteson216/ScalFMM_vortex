#!/bin/bash
#
#if ($#
filename=$1
gnuscript=~/mount/Plafrim1/Dev/src/ScalFMM/scalfmm/Utils/scripts/histTime1.gnu
FILEBASIC=`cat $filename  | grep basic`
FILETASK=`cat $filename  | grep "\-task"`
FILEBALANCED=`cat $filename  | grep balanced `
FILESECTIONTASK=`cat $filename  | grep section `

//
allfile=`cat $filename`
echo $allfile

for file in $allfile ;
do
outputeps=`basename $file  .out`"Time.eps"
TITLE=`basename $filename  .txt`

echo "gnuplot"  $outputeps
gnuplot -e "FILEEPS='$outputeps'; TITLE='$TITLE'; FILE='$file' "   $gnuscript

done
