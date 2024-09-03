module cubic_spline_interpolation
    use :: f95_lapack, only: la_gesv

    implicit none (type, external)

    type :: cubic_spline
        integer :: n
        real(8), dimension(:), allocatable :: x, y, k, c, d

        contains
            final :: csdestructor
    end type

contains

    subroutine interpolate(x, y, cs)
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(size(x)), intent(in) :: y
        type(cubic_spline), intent(out) :: cs

        real(8), dimension(:), allocatable :: c, d

        real(8), dimension(size(x), size(x)) :: M
        real(8), dimension(size(x)) :: b

        integer :: i, j, n

        n = size(x)

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

    pure subroutine binary_search(v, x, i)
        real(8), dimension(:), intent(in) :: v
        real(8), intent(in) :: x
        integer, intent(out) :: i

        integer :: n
        integer :: l, r, m

        n = size(v)
        l = 1
        r = n

        do while (l <= r)
            m = (l + r)/2
            if (v(m) == x) then
                i = m
                return
            else if (v(m) < x) then
                l = m + 1
            else
                r = m - 1
            end if
        end do

        i = l-1

    end subroutine binary_search

    !TODO: performance issue here, should do binary search
    pure subroutine evaluate(cs, x, y)
        type(cubic_spline), intent(in) :: cs
        real(8), intent(in) :: x
        real(8), intent(out) :: y
        real(8) :: t

        integer :: i

        call binary_search(cs%x, x, i)
        if (i == 0) then
            y = 0.0d0
            return
        end if

        t = (x - cs%x(i))/(cs%x(i+1) - cs%x(i))
        y = (1.0 - t)*cs%y(i) + t*cs%y(i+1) + t*(1.0 - t)*(cs%c(i)*(1.0 - t) + cs%d(i)*t)
    end subroutine evaluate

    subroutine csdestructor(self)
        type(cubic_spline), intent(inout) :: self

        if (allocated(self%x)) deallocate(self%x)
        if (allocated(self%y)) deallocate(self%y)
        if (allocated(self%k)) deallocate(self%k)
        if (allocated(self%c)) deallocate(self%c)
        if (allocated(self%d)) deallocate(self%d)
    end subroutine csdestructor

end module cubic_spline_interpolation
