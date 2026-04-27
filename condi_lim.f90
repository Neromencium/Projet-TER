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
      real(PR), intent(in) :: temps !Indice d'incrémentation'

      real(PR) :: reel

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
      
      !open(unit = 10, file = "", action = "read")        
      !read(10,*) reel, phi(2,1)
      !close(10)

      !stress_int(1) = 1e-6
      !u_int(1) = phi(2,1) + ( stress_int(1) + stress(1) ) / imped(1)

      !stress_int(Nx+1) = 1e-6
      !u_int(Nx+1) = phi(2,Nx) + ( stress_int(Nx+1) + stress(Nx) ) / imped(Nx)
       
      end if
   end subroutine condition_limites

end module cond_lim
