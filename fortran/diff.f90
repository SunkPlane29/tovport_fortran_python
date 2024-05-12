module runge_kutta_4
    implicit none

    type diffeq_solution
        real(8), dimension(:, :), allocatable :: sol
        integer :: n

    contains
        final :: diffeq_solution_destructor
    end type diffeq_solution
    
    interface
        function f(n, t, x) result(dx)
            integer, intent(in) :: n
            real(8), intent(in) :: t
            real(8), dimension(n), intent(in) :: x
            real(8), dimension(n) :: dx
        end function f

        function terminate(n, i, tn, xn) result(b)
            integer, intent(in) :: n, i
            real(8), intent(in) :: tn
            real(8), dimension(n), intent(in) :: xn
            logical :: b
        end function terminate
    end interface
contains

    pure subroutine step(f, n, t, x, h, next_x)
        interface
            pure function f(n, t, x) result(dx)
                integer, intent(in) :: n
                real(8), intent(in) :: t
                real(8), dimension(n), intent(in) :: x
                real(8), dimension(n) :: dx
            end function f
        end interface

        integer, intent(in) :: n
        real(8), intent(in) :: t, h
        real(8), dimension(n), intent(in) :: x
        real(8), dimension(n), intent(out) :: next_x

        real(8), dimension(n) :: kn1, kn2, kn3, kn4

        kn1 = f(n, t, x)
        kn2 = f(n, t + h/2, x + (h/2) * kn1)
        kn3 = f(n, t + h/2, x + (h/2) * kn2)
        kn4 = f(n, t + h, x + h * kn3)

        next_x = x + (h/6) * (kn1 + 2*kn2 + 2*kn3 + kn4)

    end subroutine step

    pure subroutine solve_diffeq(f, t0, x0, h, terminate, maxiter, sol)
        interface
            pure function f(n, t, x) result(dx)
                integer, intent(in) :: n
                real(8), intent(in) :: t
                real(8), dimension(n), intent(in) :: x
                real(8), dimension(n) :: dx
            end function f

            pure function terminate(n, i, tn, xn) result(b)
                integer, intent(in) :: n
                integer, intent(in) :: i
                real(8), intent(in) :: tn
                real(8), dimension(n), intent(in) :: xn
                logical :: b
            end function terminate
        end interface

        real(8), intent(in) :: t0, h
        real(8), dimension(:), intent(in) :: x0
        integer, intent(in) :: maxiter

        real(8) :: tn
        real(8), dimension(size(x0)) :: xn, next_x
        real(8), dimension(maxiter, size(x0)+1) :: tempsol
        type(diffeq_solution), intent(out) :: sol
        
        integer :: i, n

        n = size(x0)

        tn = t0
        xn = x0

        tempsol(1, :) = [tn, xn]

        do i = 2, maxiter
            call step(f, n, tn, xn, h, next_x)
            tn = tn + h
            if (terminate(n, i, tn, next_x)) exit
            xn = next_x

            tempsol(i, :) = [tn, xn]
        end do

        allocate(sol%sol(i, n+1))
        sol%sol = tempsol(1:i-1, :)
        sol%n = i-1

    end subroutine solve_diffeq

    subroutine diffeq_solution_destructor(self)
        type(diffeq_solution), intent(inout) :: self
        if (allocated(self%sol)) deallocate(self%sol)
    end subroutine diffeq_solution_destructor
end