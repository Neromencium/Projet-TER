module cond_init

   use precision
   use public
   use type_def
   use public

   implicit none

   private

   public :: condition_initiales

contains

   subroutine condition_initiales(phi_initial)

      implicit none
      real(kind=PR),allocatable,dimension(:,:),intent(inout) :: phi_initial
      real(kind=PR),allocatable,dimension(:) :: rho, speed_x, internal_energy, pressure
      real(kind=PR) :: pressure_L, pressure_R, dx
      real(kind=PR) :: dist
      integer:: i

      ! Choix du cas test
      !test_case = 'Sod'
      !test_case = 'Impact'
      test_case = 'Laser'

      if (test_case == 'Sod') then
         Lxmin = 0._PR
         Lxmax = 1._PR

         N = 3 ! Number of equations (we start from 1)
         Nx = 400

         Lx = Lxmax - Lxmin
         dx = Lx / Nx

          gamma = 5._PR / 3. ! monoatomic
         !gamma = 7._PR / 5._PR ! diatomic

         CFL     = 0.25_PR
         T_final = 0.2_PR

         allocate( tab_cells(1:Nx), tab_nodes(1:Nx+1), phi_initial(1:N,0:Nx+1) )
         allocate( rho(1:Nx), speed_x(1:Nx), internal_energy(1:Nx), pressure(1:Nx) )

         phi_initial(1:N,1:Nx) = 0._PR

         pressure_L = 1._PR
         do i = 1,Nx/2
               rho(i) = 1._PR
               speed_x(i) = 0._PR
               pressure(i) = pressure_L
         end do

         pressure_R = 0.1_PR
         do i = Nx/2+1, Nx
               rho(i) = 0.125_PR
               speed_x(i) = 0._PR
               pressure(i) = pressure_R
         end do

      else if (test_case == 'Impact') then

         Lxmin = 0._PR
         Lxmax = 50e-3_PR

         N = 3 ! Number of equations (we start from 1)
         Nx = 2500

         Lx = Lxmax - Lxmin
         dx = Lx / Nx

         Gamma0 = 2._PR
         c0     = 5328._PR
         s1     = 1.338_PR 
         rho0   = 2785._PR
         shear_modulus = 27.6e+9
         yield_strength = 300.0e+6

         CFL     = 0.45_PR
         T_final = 5e-6_PR

         allocate( tab_cells(1:Nx), tab_nodes(1:Nx+1), phi_initial(1:N,0:Nx+1) )
         allocate( rho(1:Nx), speed_x(1:Nx), internal_energy(1:Nx), pressure(1:Nx) )

         phi_initial(1:N,1:Nx) = 0._PR

         pressure_L = 1e-6_PR
         do i = 1,Nx/10
            rho(i) = 2785._PR
            speed_x(i) = 800._PR
            pressure(i) = pressure_L
         end do

         pressure_R = 1e-6_PR
         do i = Nx/10+1, Nx
            rho(i) = 2785_PR
            speed_x(i) = 0._PR
            pressure(i) = pressure_R
         end do

      else if (test_case == 'Laser') then

      ! print*,'TO DO :)'

         Lxmin = 0._PR
         Lxmax = 0.001030_PR  
        
         N = 3 ! Number of equations (we start from 1)
         Nx = 400

         Lx = Lxmax - Lxmin
         dx = Lx / Nx
    

         allocate( tab_cells(1:Nx), tab_nodes(1:Nx+1), phi_initial(1:N,0:Nx+1) )
         allocate( rho(1:Nx), speed_x(1:Nx), internal_energy(1:Nx), pressure(1:Nx) )
      

         !rho_init = 3900.0D0
         !pre_init = 1e+5
         speed_x = 0.0D0
         Gamma0 = 0.1D0
         C0 = 7410
         s1 = 0.0D0 
         rho0 = 3900
         shear_modulus = 1.47e+11
         yield_strength = 7e+9

         CFL     = 0.45_PR
         T_final = 5.0D-7

         phi_initial(1:N,1:Nx) = 0._PR

         do i = 1, Nx
            rho(i) = rho0
            speed_x(i) = 0._PR
            pressure(i) = 0._PR
         end do

      end if

    
      !!!!!!!!!!!!!!!!!!!
        
  
       ! i = 1
       tab_nodes(1)%x = Lxmin  
       ! other i      
       do i = 2, Nx+1
         tab_nodes(i)%x = tab_nodes(i-1)%x + dx
      end do

      ! Deduction des positions des centres des mailles
      do i = 1, Nx
         ! Calculer x(i) 
         tab_cells(i)%x = 0.5 * (tab_nodes(i)%x + tab_nodes(i+1)%x)
      end do

      ! Initialisation commune
      do i = 1,Nx
         ! pressure = (gamma - 1._PR) * rho(i) * internal_energy(i)
         ! internal_energy = pressure / ((gamma-1._PR)*rho(i))
         phi_initial(1,i) = rho(i)
         phi_initial(2,i) = speed_x(i)
         if (test_case == 'Sod') then
          phi_initial(3,i) = pressure(i) / ((gamma-1._PR)*rho(i)) + 0.5 * speed_x(i)**2
         else if (test_case == 'Impact') then
          phi_initial(3,i) = pressure(i) / (Gamma0*rho(i)) + 0.5 * speed_x(i)**2
         end if 
      end do

      ! Calcul des masses
      do i = 1,Nx
         ! Calculer dx(i) 
         tab_cells(i)%mass = rho(i) * (tab_nodes(i+1)%x - tab_nodes(i)%x)
      end do      

      open (20, file='rho_initial.dat', action='write', position='rewind', form='formatted', access='sequential')
      do i = 1, Nx
         write(20,*)tab_cells(i)%x,rho(i)
      end do
      close(20)

   end subroutine condition_initiales

end module cond_init

