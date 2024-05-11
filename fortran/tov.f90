module eos_init
    use constants
    use cubic_spline_interpolation

    implicit none

    private
    public :: eos_interp

    type(cubic_spline) :: eos_interp

    !TODO: initialize the module, open data file, convert units and interpolate data
    
    contains

    subroutine eos(P, e)
        real(8), intent(in) :: P
        real(8), intent(out) :: e

        call evaluate(eos_interp, P, e)
    end subroutine eos
end

module tov
    use constants
    use runge_kutta_4
    use eos_init

    implicit none

contains

function pressure_diffeq(r, x) result(dP)
    real(8) :: r
    real(8), dimension(2) :: x
    real(8) :: dP

    real(8) :: P, M, e

    P = x(1)
    M = x(2)
    call eos(P, e)

    dP = -(e*M/r**2)*(1 + P/e)*(1 + 4*pi*r**3*P/M)*(1 - 2*M/r)**(-1)
end function pressure_diffeq

function mass_diffeq(r, x) result(dM)
    real(8) :: r
    real(8), dimension(2) :: x
    real(8) :: dM

    real(8) :: P, M, e

    P = x(1)
    M = x(2)
    call eos(P, e)

    dM = 4*pi*r**2*e
end

function tov_diffeq(n, r, x) result(dx)
    integer :: n
    real(8) :: r
    real(8), dimension(n) :: x
    real(8), dimension(n) :: dx

    dx(1) = pressure_diffeq(r, x)
    dx(2) = mass_diffeq(r, x)
end function tov_diffeq

subroutine solve_tov(p0, h, sol)
    real(8) :: p0
    real(8) :: h
    type(diffeq_solution) :: sol

    real(8) :: t0
    real(8), dimension(2) :: x0

    t0 = 1.0d-8
    x0(1) = p0
    x0(2) = 1.0d-24

    call sove_diffeq(tov_diffeq, 2, t0, x0, h, sol)
end

end module tov