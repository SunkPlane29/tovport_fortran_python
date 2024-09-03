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

    p0i = 1.0d0 * MEVFM3_TO_PRESSURE_UNIT
    p0f = 600.0d0 * MEVFM3_TO_PRESSURE_UNIT
    delp0 = (p0f - p0i) / n
    do i = 1, n
        p0(i) = p0i + delp0 * (i - 1)
    end do

    call read_eos("../eos.csv")

    call solve_mrdiagram(n, p0, mrdiagram)

    open(newunit=unit, file="out/mrdiagram.dat")
    do i = 1, n
        write(unit, *) mrdiagram(i, 1), mrdiagram(i, 2), mrdiagram(i, 3)
    end do
    close(unit)

end program main
