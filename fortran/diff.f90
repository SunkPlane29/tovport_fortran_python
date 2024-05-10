module runge_kutta_4
    implicit none
    
contains

    subroutine step(f, t, x, h, next_x)
        interface
            function f(t, x) result(dx)
                real(8), intent(in) :: t
                real(8), dimension(:), intent(in) :: x
                real(8), dimension(size(x)) :: dx
            end function f
        end interface

        real(8), intent(in) :: t, h
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(size(x)) :: next_x

        real(8), dimension(size(x)) :: kn1, kn2, kn3, kn4

        kn1 = f(t, x)
        kn2 = f(t + h/2, x + (h/2) * kn1)
        kn3 = f(t + h/2, x + (h/2) * kn2)
        kn4 = f(t + h, x + h * kn3)

        next_x = x + (h/6) * (kn1 + 2*kn2 + 2*kn3 + kn4)

    end subroutine step

    subroutine solve(f, t0, x0, h, terminate, maxiter, sol)
        interface
            function f(t, x) result(dx)
                real(8), intent(in) :: t
                real(8), dimension(:), intent(in) :: x
                real(8), dimension(size(x)) :: dx
            end function f

            function terminate(i, tn, xn) result(b)
                integer, intent(in) :: i
                real(8), intent(in) :: tn
                real(8), dimension(:), intent(in) :: xn
                logical :: b
            end function terminate
        end interface

        real(8), intent(in) :: t0, h
        real(8), dimension(:), intent(in) :: x0
        integer, intent(in) :: maxiter

        real(8) :: tn
        real(8), dimension(size(x0)) :: xn, next_x
        real(8), dimension(maxiter, size(x0)+1) :: tempsol
        real(8), dimension(:, :), allocatable, intent(out) :: sol
        
        integer :: i, n

        n = size(x0)

        tn = t0
        xn = x0

        sol(1, :) = [tn, xn]

        do i = 2, maxiter
            call step(f, tn, xn, h, next_x)
            tn = tn + h
            if (terminate(i, tn, next_x)) exit
            xn = next_x

            sol(i, :) = [tn, xn]
        end do

        allocate(sol(i, n+1))
        sol(1:i, :) = tempsol(1:i, :)

    end subroutine solve
end