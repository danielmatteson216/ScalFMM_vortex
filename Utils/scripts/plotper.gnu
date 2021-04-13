
#set terminal postscript enhanced color 
#set output FMMEnergy.ps'
set key left
set encoding iso_8859_1 

set title "FMM Energy"
set xlabel "per size (per)"
set ylabel "Energy (kcal/mol)"
 plot "FMM_PER.out" u ($1,$3) w l t "FMM Energy"

