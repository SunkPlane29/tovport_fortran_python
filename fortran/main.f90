program main
    use constants
    use eos_init
    use tov
    implicit none

    integer, parameter :: n = 200
    real(8), dimension(n, 3) :: mrdiagram
    real(8), dimension(n) :: p0

    real(8) :: delp0, p0i, p0f
    integer :: i, unit

    real(8), dimension(10) :: testv

    p0i = 5.0d0
    p0f = 800.0d0
    delp0 = (log(p0f) - log(p0i)) / n
    do i = 1, n
        p0(i) = log(p0i) + delp0 * (i - 1)
    end do
    p0 = exp(p0)*MEVFM3_TO_PRESSURE_UNIT

    call read_eos("../eos_1.csv")

    call solve_mrdiagram(n, p0, mrdiagram)

    open(newunit=unit, file="out/mrdiagram.dat")
    do i = 1, n
        write(unit, *) mrdiagram(i, 1), mrdiagram(i, 2), mrdiagram(i, 3)
    end do
    close(unit)

end program main
