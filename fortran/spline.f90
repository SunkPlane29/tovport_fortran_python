module cubic_spline_interpolation
    use :: f95_lapack, only: la_gesv

    implicit none (type, external)

    external :: sgesv

    type :: cubic_spline
        integer :: n
        real(8), dimension(:), allocatable :: x, y, k, c, d
    end type

contains

    subroutine interpolate(x, y, n, cs)
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: x, y
        type(cubic_spline), intent(out) :: cs

        real(8), dimension(:), allocatable :: c, d
        
        real(8), dimension(n, n) :: M
        real(8), dimension(n) :: b
        
        integer :: i, j
        
        allocate(c(n-1), d(n-1))
        allocate(cs%x(n), cs%y(n), cs%k(n), cs%c(n-1), cs%d(n-1))
            
        M = 0.0
        b = 0.0

        M(1, 1) = 2.0/(x(2) - x(1))
        M(1, 2) = 1.0/(x(2) - x(1))
        b(1) = 3.0*(y(2) - y(1))/(x(2) - x(1))**2

        do i = 2, n-1
            do j = i-1, i+1
                if (j == i-1) then
                    M(i, j) = 1.0/(x(i) - x(i-1))
                else if (j == i) then
                    M(i, j) = 2.0*(1/(x(i) - x(i-1)) + 1/(x(i+1) - x(i)))
                else if (j == i+1) then
                    M(i, j) = 1.0/(x(i+1) - x(i))
                end if
            end do
            b(i) = 3.0*(y(i+1) - y(i))/(x(i+1) - x(i))**2 + 3.0*(y(i) - y(i-1))/(x(i) - x(i-1))**2
        end do

        M(n, n-1) = 1.0/(x(n) - x(n-1))
        M(n, n) = 2.0/(x(n) - x(n-1))
        b(n) = 3.0*(y(n) - y(n-1))/(x(n) - x(n-1))**2

        call la_gesv(M, b)

        c = 0.0
        d = 0.0

        do i = 1, n-1
            c(i) = b(i)*(x(i+1) - x(i)) - (y(i+1) - y(i))
            d(i) = -b(i+1)*(x(i+1) - x(i)) + (y(i+1) - y(i))
        end do

        cs%n = n
        cs%x = x
        cs%y = y
        cs%k = b
        cs%c = c
        cs%d = d

    end subroutine interpolate

    subroutine evaluate(cs, x, y)
        type(cubic_spline), intent(in) :: cs
        real(8), intent(in) :: x
        real(8), intent(out) :: y
        real(8) :: t

        integer :: i

        do i = 1, cs%n-1
            if (x >= cs%x(i) .and. x <= cs%x(i+1)) then
                t = (x - cs%x(i))/(cs%x(i+1) - cs%x(i))
                y = (1.0 - t)*cs%y(i) + t*cs%y(i+1) + t*(1.0 - t)*(cs%c(i)*(1.0 - t) + cs%d(i)*t)
                return
            end if
        end do

        y = 0.0
    end subroutine evaluate

end module cubic_spline_interpolation