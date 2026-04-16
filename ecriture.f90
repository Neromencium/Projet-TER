module ecri

   use precision
   use public
   use type_def
   use eos_dt

   implicit none

   private

   public :: ecriture

contains

   subroutine ecriture( phi_final )
      implicit none
      real(kind=PR), dimension(1:N,1:Nx),intent(in):: phi_final
      integer:: i

      ! Deduction des positions des centres des mailles
      do i = 1, Nx
         ! Calculer x(i) 
         tab_cells(i)%x = 0.5 * (tab_nodes(i)%x + tab_nodes(i+1)%x)
      end do

      open (20, file='rho.dat', action='write', position='rewind', form='formatted', access='sequential')
      do i = 1, Nx
         write(20,*)tab_cells(i)%x,phi_final(1,i)
      end do
      close(20)

      open (20, file='speed.dat', action='write', position='rewind', form='formatted', access='sequential')
      do i = 1, Nx
         write(20,*)tab_cells(i)%x,phi_final(2,i)
      end do
      close(20)

      open (20, file='pressure.dat', action='write', position='rewind', form='formatted', access='sequential')
      do i = 1, Nx
         write(20,*)tab_cells(i)%x,pressure(phi_final(1:N,i))
      end do
      close(20)

   end subroutine ecriture


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ecri
