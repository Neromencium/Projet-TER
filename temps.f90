module temps

   use precision
   use public
   use cond_lim
   use eos_dt
   use ecri

   implicit none

   private

   public :: boucle_temps

   
contains

   subroutine boucle_temps( phi )

      real(kind=PR),dimension(1:N,0:Nx+1),intent(inout):: phi
      real(kind=PR):: dt, temps, radial_return
      real(kind=PR),dimension(1:Nx+1) :: u_int       ! Tableau des vitesses nodales
      real(kind=PR),dimension(1:Nx+1) :: stress_int  ! Tableau des pressions aux interfaces
      real(kind=PR),dimension(1:Nx) :: imped, stress, deviateur_stress ! Tableau des impédances, contraintes et deviateur des contraintes
  
      integer:: kk, i
      real(PR) :: Delta_x

      character(len=100) :: filename

      temps   = 0._PR
      kk      = 0
      dt      = 1e-12
      deviateur_stress(:) = 0.0 
      cpt = 0


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !do kk = 1,1000
      do while ( temps .lt. T_final )

         ! Calcul des impedances et contraintes
         do i = 1, Nx
          imped(i)  = phi(1,i) * sound_speed(phi(1:N,i))
          stress(i) = - pressure(phi(1:N,i)) + deviateur_stress(i)
         end do

         ! Calcul condition limite
         call condition_limites(phi(1:N,1:Nx), imped(1:Nx), stress(1:Nx), u_int(1:Nx+1), stress_int(1:Nx+1), kk)

         ! Calcul des flux numériques
         do i = 2, Nx
          u_int(i)      = ( (imped(i)*phi(2,i) + imped(i-1)*phi(2,i-1)) + (stress(i)-stress(i-1)) ) / (imped(i)+imped(i-1))
          stress_int(i) = -(imped(i)*stress(i-1) + imped(i-1)*stress(i) + imped(i)*imped(i-1)*(phi(2,i) - phi(2,i-1))) / (imped(i)+imped(i-1))
         end do

         ! Calcul pas de temps
         call calcul_dt(dt, phi(1:N,1:Nx) )

         ! Mise à jour des positions des noeuds
         do i = 1, Nx+1
          tab_nodes(i)%x = tab_nodes(i)%x + dt * u_int(i)
         end do

         ! Calcul des nouvelles densités
         do i = 1, Nx
          phi(1,i) = tab_cells(i)%mass / (tab_nodes(i+1)%x - tab_nodes(i)%x)
         end do
 
         ! Calcul des viteses et des energies totales
         do i = 1, Nx
           phi(2,i) = phi(2,i) - (dt / tab_cells(i)%mass) * (stress_int(i+1) - stress_int(i)) 
           phi(3,i) = phi(3,i) - (dt / tab_cells(i)%mass) * (stress_int(i+1) * u_int(i+1) - stress_int(i) * u_int(i)) 
         end do

         ! Prise en compte de la plasticité (ne pas modifier) (Critere de Von Mises pris en compte par methode de retour radial)
         do i = 1, Nx
          radial_return = abs(deviateur_stress(i)) / (sqrt(2.0/3.0)*yield_strength)
          if (radial_return > 1.0) then 
           deviateur_stress(i) = deviateur_stress(i) / radial_return
          end if
         end do 

        
         !!!!!!!!!!!!!!!!!!!
        !Partie pour calculer la solution

         do i = 2,Nx
             Delta_x = tab_nodes(i+1)%x - tab_nodes(i-1)%x 
             deviateur_stress(i) = deviateur_stress(i) + 4._PR/3 * 27.6e+9 * dt/Delta_x * (u_int(i+1)-u_int(i-1))             
         end do
                
         !!!!!!!!!!!!!!!!!!!
         !Partie pour écrire les fichiers

         if (kk == 100*cpt) then
             write(filename,'("data_sol/sol_approx_", I3.3, ".dat")') cpt + 1
             open(unit = 20, file = filename, status = "replace", action = "write")
             do i = 2, SIZE(deviateur_stress)
                 write(20,*) tab_cells(i)%x, deviateur_stress(i)          
             end do
             close(20)
             
             cpt = cpt + 1

         end if        
           

         !!!!!!!!!!!!!!!!!!!

         temps = temps + dt
         kk = kk + 1
         print*,"Physical time:",temps,'Time step:',dt,"Iteration:",kk         

      end do
      

      !Ecriture au temps final
      write(filename,'("data_sol/sol_approx_", I3.3, ".dat")') cpt + 1
      open(unit = 20, file = filename, status = "replace", action = "write")
      do i = 2, SIZE(deviateur_stress)
          write(20,*) tab_cells(i)%x, deviateur_stress(i)             
      end do
      close(20)


   end subroutine boucle_temps

end module temps
