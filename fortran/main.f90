program main
    use cubic_spline_interpolation
    implicit none

    real(8), dimension(3) :: x, y
    integer :: n
    type(cubic_spline) :: cs
    real(8) :: y_interp

    integer :: i

    x = [-1.0, 0.0, 3.0]
    y = [0.5, 0.0, 3.0]
    n = size(x)

    call interpolate(x, y, n, cs)
    call evaluate(cs, 2.0d0, y_interp)

    print '("y_interp = ", f0.4)', y_interp
    
end program main