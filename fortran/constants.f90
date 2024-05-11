module constants
    real(8), parameter :: si_c = 2.99792458d8
    real(8), parameter :: si_G = 6.67430d-11
    real(8), parameter :: si_MSOLAR = 1.98855d30
    real(8), parameter :: si_hbar = 6.582119569d-16

    real(8), parameter :: pi = atan(1.d0)*4.d0

    real(8), parameter :: MASS_UNIT_TO_SI = si_MSOLAR
    real(8), parameter :: LENGTH_UNIT_TO_SI = si_G*si_MSOLAR/(si_c**2)
    real(8), parameter :: TIME_UNIT_TO_SI = si_G*si_MSOLAR/(si_c**3)
    real(8), parameter :: SI_TO_MASS_UNIT = 1.d0/MASS_UNIT_TO_SI
    real(8), parameter :: SI_TO_LENGTH_UNIT = 1.d0/LENGTH_UNIT_TO_SI
    real(8), parameter :: SI_TO_TIME_UNIT = 1.d0/TIME_UNIT_TO_SI

    real(8), parameter :: PRESSURE_UNIT_TO_SI = MASS_UNIT_TO_SI*(LENGTH_UNIT_TO_SI**(-1))*(TIME_UNIT_TO_SI**(-2))
    real(8), parameter :: SI_TO_PRESSURE_UNIT = 1.d0/PRESSURE_UNIT_TO_SI

    real(8), parameter :: hbarc = (si_c*10.0d0**15)*(si_hbar*10.0d0**(-6))

    real(8), parameter :: FM4_TO_MEV4 = hbarc**4
    real(8), parameter :: MEV4_TO_FM4 = 1.d0/FM4_TO_MEV4

    real(8), parameter :: FM4_TO_MEVFM3 = hbarc**3
    real(8), parameter :: MEVFM3_TO_FM4 = 1.d0/FM4_TO_MEVFM3

    real(8), parameter :: MEV4_TO_MEVFM3 = hbarc**(-3)
    real(8), parameter :: MEVFM3_TO_MEV4 = 1.d0/MEV4_TO_MEVFM3

    real(8), parameter :: EV_TO_JOULE = 1.602176634d-19
    real(8), parameter :: JOULE_TO_EV = 1.d0/EV_TO_JOULE

    real(8), parameter :: MEV4_TO_JOULEM3 = ((si_hbar*10.0d0**(-6))*si_c)**(-3) * (EV_TO_JOULE*10.0d0**6)
    real(8), parameter :: JOULEM3_TO_MEV4 = 1.d0/MEV4_TO_JOULEM3

    real(8), parameter :: MEVFM3_TO_PRESSURE_UNIT = MEVFM3_TO_MEV4 * MEV4_TO_JOULEM3 * SI_TO_PRESSURE_UNIT
    real(8), parameter :: PRESSURE_UNIT_TO_MEVFM3 = 1.d0/MEVFM3_TO_PRESSURE_UNIT
end module constants