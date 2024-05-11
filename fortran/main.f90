program main
    use constants
    use eos_init
    use tov
    implicit none

    type(diffeq_solution) :: tov_sol
    integer :: i, unit

    call read_eos("../eos.csv")
    call solve_tov(1.0d0*MEVFM3_TO_PRESSURE_UNIT, tov_sol)

    open(newunit=unit, file="out/tov.dat")
    do i = 1, tov_sol%n
        write(unit, *) tov_sol%sol(i, :)
    end do
    close(unit)
    
end program main