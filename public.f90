module public

   use precision
   use type_def

   implicit none

   private

   real(kind=PR), parameter, public:: pi = 2._PR * acos(0._PR)
   integer, public :: N, Nx
   real(kind=PR),public :: Lxmin, Lxmax, Lx, dx
   type(Cell),dimension(:), allocatable, public :: tab_cells
   type(Node),dimension(:), allocatable, public :: tab_nodes
   real(kind=PR),public :: CFL, T_final
   real(kind=PR),public :: gamma, Gamma0, c0, s1, rho0, shear_modulus, yield_strength
   character(50), public :: test_case
   integer, public :: cpt

end module public
