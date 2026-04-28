module cond_lim

   use precision
   use public
   use type_def
   use eos_dt

   implicit none

   private

   public :: condition_limites

contains

   subroutine condition_limites( phi, imped, stress, u_int, stress_int, temps)
      implicit none
      real(kind=PR),dimension(1:N,1:Nx),intent(in):: phi
      real(kind=PR),dimension(1:Nx),intent(in):: imped, stress ! Tableau des impedances et des contraintes
      real(kind=PR),dimension(1:Nx+1),intent(inout) :: u_int ! Tableau des vitesses nodales
      real(kind=PR),dimension(1:Nx+1),intent(inout) :: stress_int ! Tableau des pressions aux interfaces
      real(PR), intent(in) :: temps
      
      integer :: compteur, j
      real(PR) :: reel_t, reel_p, reel_t1, reel_p1

      compteur = 0

      if (test_case == 'Sod') then
      ! Bord gauche 
      ! Condition de vitesse nulle imposee (ici vitesse nulle) 
      u_int(1) = 0.0
      stress_int(1) = - stress(1) - imped(1) * ( u_int(1) - phi(2,1) )
      ! Bord droit       
      ! Condition de vitesse nulle imposee (ici vitesse nulle) 
      u_int(Nx+1) = 0.0
      stress_int(Nx+1) = -stress(Nx) + imped(Nx) * ( u_int(Nx+1) - phi(2,Nx) )

      else if (test_case == 'Impact') then
      ! Bord gauche 
      ! Condition de pression imposee 
      stress_int(1) = 1e-6
      u_int(1) = phi(2,1) + ( stress_int(1) + stress(1) ) / imped(1)
      ! Bord droit       
      ! Condition de pression imposee 
      stress_int(Nx+1) = 1e-6
      u_int(Nx+1) = phi(2,Nx) + ( stress_int(Nx+1) + stress(Nx) ) / imped(Nx)

      else if (test_case == 'Laser') then
      

      open(unit = 10, file = "pression_laser.dat", action = "read")        
      do j = 1, compteur 
         read(10,*) reel_t, reel_p
      end do
      if (j <= 1502) then
         read(10,*) reel_t1, reel_p1
      end if
      close(10)


      if (temps > reel_t1) then
         stress_int(1) = reel_p1
         compteur = compteur + 1
      else
         stress_int(1) = reel_p
      end if
      

      u_int(1) = phi(2,1) + ( stress_int(1) + stress(1) ) / imped(1)

      stress_int(Nx+1) = 0
      u_int(Nx+1) = phi(2,Nx) + ( stress_int(Nx+1) + stress(Nx) ) / imped(Nx)
       
      end if
   end subroutine condition_limites

end module cond_lim
