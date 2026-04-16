module eos_dt

   use precision
   use public
   use type_def

   implicit none

   private

   public :: calcul_dt, pressure, sound_speed

contains

   subroutine calcul_dt(dt, phi)
      real(kind=PR),intent(inout):: dt
      real(kind=PR),dimension(1:N,1:Nx),intent(in):: phi
      real(kind=PR) :: dt_min, caracteristical_length
      integer :: i

      dt_min = 1e+0_PR
      do i = 1,Nx
         caracteristical_length = abs(tab_nodes(i+1)%x - tab_nodes(i)%x)
         dt_min = min( dt_min, CFL * caracteristical_length / sound_speed(phi(1:N,i)) )
      end do
      dt_min = min(dt_min, 1.05 * dt)
      dt = dt_min

   end subroutine calcul_dt

   function pressure(phi)
      implicit none
      real(kind=PR), dimension(1:N), intent(in) :: phi
      real(kind=PR) :: pressure, rho, internal_energy, mu, denom, num
      rho = phi(1)
      internal_energy = phi(N) - 0.5 * phi(2)**2
      if (test_case == 'Sod') then ! Model gaz parfait pour Sod
       pressure = (gamma - 1._PR) * rho * internal_energy
      else if (test_case == 'Impact') then ! Equation d'état Mie-gruneisen pour cas test impact
       mu = rho / rho0 - 1.0D0
       if (mu > 0) then ! compression
        num   =  rho0 * c0 * c0 * mu * (1.0 + (1.0-0.5*Gamma0) * mu)
        denom = (1.0 - (s1-1.0) * mu)**2 
        pressure = (num / denom) + rho0 * Gamma0 * internal_energy
       else ! traction
        pressure = rho0 * c0 * c0 * mu + rho0 * Gamma0 * internal_energy
       end if
      end if 
   end function pressure

   function sound_speed(phi)
      implicit none
      real(kind=PR), dimension(1:N), intent(in) :: phi
      real(kind=PR) :: sound_speed, rho, mu, denom, num
      rho = phi(1)
      if (test_case == 'Sod') then
       sound_speed = sqrt( gamma * pressure(phi(1:N)) / rho )
      else if (test_case == 'Impact') then
       mu    = rho / rho0 - 1.0D0
       num   = (mu+1) + (s1-Gamma0) * mu
       denom = (mu+1-s1*mu)**3
       sound_speed = c0 * c0 * (num/denom) + Gamma0 * pressure(phi(1:N)) / (rho0 * (mu+1.)**2)
       sound_speed = sqrt(sound_speed)
      end if
   end function sound_speed

end module eos_dt
