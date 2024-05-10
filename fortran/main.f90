program main
    use cubic_spline_interpolation
    use runge_kutta_4
    implicit none

    interface
        function f(t, x) result(dx)
            real(8), intent(in) :: t
            real(8), dimension(2), intent(in) :: x
            real(8), dimension(size(x)) :: dx
        end function f

        function terminate(i, tn, xn) result(b)
            integer, intent(in) :: i
            real(8), intent(in) :: tn
            real(8), dimension(2), intent(in) :: xn
            logical :: b
        end function terminate
    end interface

    real(8) :: t0, h
    real(8), dimension(2) :: x0
    integer :: maxiter
    real(8), dimension(:, :), allocatable :: sol

    t0 = 0.0d0
    x0 = [20.0d0, 20.0d0]
    h = 0.01d0
    maxiter = 100000
    call solve(f, t0, x0, h, maxiter, terminate, sol)

end program main

function f(t, x) result(dx)
    implicit none

    real(8), intent(in) :: t
    real(8), intent(in), dimension(2) :: x
    real(8), dimension(2) :: dx

    real(8) :: alpha, beta

    alpha = 0.01
    beta = 0.02

    dx(1) = (1 - alpha*x(2))*x(1)
    dx(2) = (-1 + beta*x(1))*x(2)
end

function terminate(i, tn, xn) result(b)
    implicit none

    integer, intent(in) :: i
    real(8), intent(in) :: tn
    real(8), dimension(2), intent(in) :: xn
    logical :: b

    b = tn > 15 
end