 set yrange [-3e8:3e8]
 set xlabel "Position x (en m)"
 set ylabel "Déviateur tenseur s (en N.m-2)"
do for [i=1:3] {
     filename = sprintf('data_sol/sol_approx_%03d.dat', i)
     plot filename using 1:2 with lines title "Solution"
     pause 0.1
 }
 pause -1
