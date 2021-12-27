set grid
set ylabel "RHS"
set xlabel "x grid "
plot "RHS.dat" w lp title "RHS" ,  "SOL.dat" w lp title "SOL" , "EX_SOL.dat" w lp title "EX_SOL"
pause -1
