module langevin
use math
implicit none

double precision, parameter :: kB = 1.38064852d-23 ! J/K
double precision, parameter :: muB = 9.274009994d-24 ! J/T
double precision, parameter :: muB_by_kB = 0.6717140430498562 ! K/T

contains
  double precision function langevinSingleValue(B, p, Np)
    !Langevin: Spherical Homogenously Magnetized Particle with large moment
    double precision, intent(in) :: B
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: Ms, xi, x
    Ms = p(1)
    xi = p(2)

    if (abs(B) > 1e-9) then
      x = xi*B
      langevinSingleValue = Ms * (1d0/tanh(x) - 1d0/x)
    else
      langevinSingleValue = 0d0
    end if
  end function langevinSingleValue

  double precision function langevinMuSingleValue(B, p, Np)
    double precision, intent(in) :: B
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np
    double precision :: Ms, mu, T, xi

    Ms = p(1)
    mu = p(2)
    T = p(3)

    xi = mu/T*muB_by_kB
    langevinMuSingleValue = langevinSingleValue(B, (/Ms, xi/), 2)
  end function langevinMuSingleValue

  subroutine magnetization(B, Ms, mu, T, sigMu, NB, M)
    double precision, dimension(NB), intent(in) :: B
    double precision, intent(in) :: Ms, mu, T, sigMu
    integer, intent(in) :: NB
    double precision, dimension(NB), intent(out) :: M

    integer, parameter :: Np = 3
    double precision, dimension(Np) :: p

    double precision :: mu_min, mu_max
    integer :: iB

    call get_cutoff_lognormal(mu, sigMu, mu_min, mu_max)

    p = (/Ms, mu, T/)
    !$omp parallel
    !$omp do
    do iB=1, NB
      call integrate_size_distribution(B(iB), p, Np, &
            2, mu_min, mu_max, sigMu, &
            langevinMuSingleValue, lognormal,&
            M(iB))
    end do
    !$omp end do
    !$omp end parallel
  end subroutine magnetization
end module langevin

