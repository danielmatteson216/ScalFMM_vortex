#!/bin/bash

# @SCALFMM_PRIVATE

if [ "$#" -ne 1 ]; then
    echo "Use current directory as source"
    here="."
else
    echo "Use first parameter as source"
    here=$1
fi    

echo "Makes the files included in $here and below public"

for ext in "cpp" "hpp" "h" "c" ; do
    allfiles=`find $here -name \*.$ext -print`
    for thefile in $allfiles ; do
        echo "Update : $thefile"
        grep -v "@SCALFMM_PRIVATE" $thefile > $thefile.$$.tmp
        mv $thefile.$$.tmp $thefile
    done
done    
