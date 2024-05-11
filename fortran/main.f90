program main
    use cubic_spline_interpolation
    use runge_kutta_4
    implicit none

    real(8) :: t0, h
    real(8), dimension(2) :: x0
    integer :: maxiter
    real(8), dimension(:, :), allocatable :: sol

    integer :: i

    t0 = 0.0d0
    x0 = [20.0d0, 20.0d0]
    h = 0.01d0
    maxiter = 100000
    call solve(f, t0, x0, h, terminate, maxiter, sol)

    open(unit=10, file='data.dat')

    print *, size(sol(:, 1))
    
    do i = 1, size(sol(:, 1))
        write(10, *) sol(i, 1), sol(i, 2), sol(i, 3)
    end do

    close(10)

end program main

function f(n, t, x) result(dx)
    implicit none

    integer :: n
    real(8) :: t
    real(8), dimension(n) :: x
    real(8), dimension(n) :: dx

    real(8) :: alpha, beta

    alpha = 0.01
    beta = 0.02

    dx(1) = (1 - alpha*x(2))*x(1)
    dx(2) = (-1 + beta*x(1))*x(2)
end

function terminate(n, i, tn, xn) result(b)
    implicit none

    integer :: n, i
    real(8) :: tn
    real(8), dimension(n) :: xn
    logical :: b

    if (tn > 15.0d0) then
        b = .true.
    else
        b = .false.
    end if
end