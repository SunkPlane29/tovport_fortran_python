module eos_init
    use constants
    use cubic_spline_interpolation

    implicit none

    type eos_data
        real(8), dimension(:, :), allocatable :: data
        integer :: n
    contains
        final :: eos_data_destructor
    end type eos_data

    private :: eos_file
    public :: eos_interp

    type(cubic_spline) :: eos_interp
    type(eos_data) :: eos_file

    !TODO: initialize the module, open data file, convert units and interpolate data

    contains

    subroutine read_eos(eos_filename)
        character(len=100), intent(in) :: eos_filename
        integer :: i, n, iostat, unit
        real :: P, e
        character(len=100) :: line

        open(newunit=unit, file=eos_filename, status='old', action='read')

        iostat = 0
        n = -1
        do while (iostat == 0)
            n = n + 1
            read(unit, *, iostat=iostat) P, e
        end do

        eos_file%n = n

        allocate(eos_file%data(2, n))
        rewind(unit)

        read(unit, *) eos_file%data

        eos_file%data = transpose(eos_file%data)
        eos_file%data(:, 1) = eos_file%data(:, 1)*MEVFM3_TO_PRESSURE_UNIT
        eos_file%data(:, 2) = eos_file%data(:, 2)*MEVFM3_TO_PRESSURE_UNIT

        close(unit)

        call interpolate(eos_file%data(:, 1), eos_file%data(:, 2), eos_interp)

    end subroutine read_eos

    pure subroutine eos(P, e)
        real(8), intent(in) :: P
        real(8), intent(out) :: e

        call evaluate(eos_interp, P, e)
    end subroutine eos

    subroutine eos_data_destructor(self)
        type(eos_data), intent(inout) :: self
        if (allocated(self%data)) deallocate(self%data)
    end subroutine eos_data_destructor
end

module tov
    use constants
    use runge_kutta_4
    use eos_init

    implicit none

contains

pure function pressure_diffeq(r, x) result(dP)
    real(8), intent(in) :: r
    real(8), dimension(2), intent(in) :: x
    real(8) :: dP

    real(8) :: P, M, e

    P = x(1)
    M = x(2)
    call eos(P, e)

    dP = -(e*M/r**2)*(1 + P/e)*(1 + 4*pi*r**3*P/M)*(1 - 2*M/r)**(-1)
end function pressure_diffeq

pure function mass_diffeq(r, x) result(dM)
    real(8), intent(in) :: r
    real(8), dimension(2), intent(in) :: x
    real(8) :: dM

    real(8) :: P, M, e

    P = x(1)
    M = x(2)
    call eos(P, e)

    dM = 4*pi*r**2*e
end

pure function tov_diffeq(n, r, x) result(dx)
    integer, intent(in) :: n
    real(8), intent(in) :: r
    real(8), dimension(n), intent(in) :: x
    real(8), dimension(n) :: dx

    dx(1) = pressure_diffeq(r, x)
    dx(2) = mass_diffeq(r, x)
end function tov_diffeq

pure function tov_terminate(n, i, r, x) result(b)
    integer, intent(in) :: n, i
    real(8), intent(in) :: r
    real(8), dimension(n), intent(in) :: x
    logical :: b

    if (x(1) <= 0.0d0) then
        b = .true.
    else
        b = .false.
    end if
end function tov_terminate

pure subroutine solve_tov(p0, sol)
    real(8), intent(in) :: p0
    type(diffeq_solution), intent(out) :: sol

    real(8) :: t0, h
    real(8), dimension(2) :: x0
    integer :: maxiter

    t0 = 1.0d-8
    x0(1) = p0
    x0(2) = 1.0d-24
    h = 1.0d0*SI_TO_LENGTH_UNIT
    maxiter = 100000

    call solve_diffeq(tov_diffeq, t0, x0, h, tov_terminate, maxiter, sol)

    sol%sol(:, 1) = sol%sol(:, 1)*LENGTH_UNIT_TO_SI*1.0d-3
    sol%sol(:, 2) = sol%sol(:, 2)*PRESSURE_UNIT_TO_MEVFM3
end subroutine solve_tov

pure function calculate_mr(p0) result(mr)
    real(8), intent(in) :: p0
    real(8), dimension(2) :: mr

    type(diffeq_solution) :: sol

    call solve_tov(p0, sol)

    mr(1) = sol%sol(sol%n, 3)
    mr(2) = sol%sol(sol%n, 1)
end function calculate_mr

subroutine solve_mrdiagram(n, p0, mrdiagram)
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: p0
    real(8), dimension(n, 3), intent(out) :: mrdiagram
    real(8), dimension(2) :: mr

    integer :: i

    mrdiagram = 0.0d0

    do concurrent (i = 1:n)
        mrdiagram(i, :) = [p0(i), calculate_mr(p0(i))]
    end do
end subroutine solve_mrdiagram

end module tov
