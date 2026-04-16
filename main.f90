program Lagrange

   use precision
   use public
   use temps
   use cond_init
   use ecri

   implicit none
   real(kind=PR), allocatable, dimension(:,:) :: phi
   real(PR), allocatable, dimension(:,:) :: s_solution
   real(PR), allocatable, dimension(:) :: tab_temps

   call condition_initiales( phi )

   call boucle_temps( phi(1:N,0:Nx+1) )

   call ecriture( phi(1:N,1:Nx) )

    

   open(unit = 20, file = "anim.gnu", status = "replace", action = "write")

   write(20,*) "set yrange [-3e8:3e8]" 
   write(20,*) 'set xlabel "Position x (en m)"'
   write(20,*) 'set ylabel "Déviateur tenseur s (en N.m-2)"'        
   write(20,'(A,I0,A)') "do for [i=1:", cpt,"] {"
   write(20,*) "    filename = sprintf('data_sol/sol_approx_%03d.dat', i)"
   write(20,*) '    plot filename using 1:2 with lines title "Solution"'
   write(20,*) "    pause 0.1"
   write(20,*) "}"
   write(20,*) "pause -1"

   close(20)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   open(unit = 20, file = "gif.gnu", status = "replace", action = "write")

   write(20,*) "set terminal gif animate delay 10 loop 0"
   write(20,*) "set output 'animation.gif'"

   write(20,*) "set yrange [-3e8:3e8]"  
   write(20,*) 'set xlabel "Position x (en m)"'
   write(20,*) 'set ylabel "Déviateur tenseur s (en N.m-2)"'    
   write(20,'(A,I0,A)') "do for [i=1:", cpt,"] {"
   write(20,*) "    filename = sprintf('data_sol/sol_approx_%03d.dat', i)"
   write(20,*) '    plot filename using 1:2 with lines title "Solution"'
   write(20,*) "}"

   write(20,*) "unset output"

   close(20)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   open(unit = 20, file = "trace_rho.gnu", status = "replace", action = "write")

   write(20,*) 'set xlabel "Position x (en m)"'
   write(20,*) 'set ylabel "Densité rho (en kg.m-3)"'  
   write(20,*) 'plot "rho.dat" using 1:2 with lines title "Densité"'
   write(20,*) "pause -1"

   close(20)


end program Lagrange
