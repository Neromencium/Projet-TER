module type_def
   use precision
   implicit none
   !type declaration
   type Cell
      real(kind=PR) :: x, mass
   end type Cell
   type Node
      real(kind=PR) :: x
   end type Node

end module type_def
