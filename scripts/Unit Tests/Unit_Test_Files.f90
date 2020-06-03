program Write_Unit_Test_Files
    implicit none
    integer            :: ud_file, IFAIL, i
    double precision   :: x, y, z
    integer            :: step, species

    double precision, parameter :: length_scale = 1.0d-9 ! Length scale (1 nanometer)
    double precision, parameter :: time_scale = 1.0d-12 ! Time scale (1 ps)

    integer, parameter :: species_elec     = 1 ! Electron
    integer, parameter :: species_hole     = 2 ! Hole

    integer, parameter :: N = 100

    double precision :: len = 100.d0*length_scale
    double precision :: d = 1000.0d0*length_scale

    !------------------------------------------------------------------------------------------------------------

    open(newunit=ud_file, iostat=IFAIL, file='Unit_Test_1.bin', &
    status='REPLACE', action='WRITE', access='STREAM')

    x = 1.0d-9
    y = -2.34d-9
    z = 1.0d-9
    step = 1
    species = species_elec
    write(unit=ud_file) x, y, z, step, species

    x = -4.4d-9
    y = -1.44d-9
    z = 1.0d-9
    step = 10
    species = species_elec
    write(unit=ud_file) x, y, z, step, species

    close(unit=ud_file, status='keep')

    !------------------------------------------------------------------------------------------------------------

    call random_init(.false., .false.)

    open(newunit=ud_file, iostat=IFAIL, file='Unit_Test_rand.bin', &
    status='REPLACE', action='WRITE', access='STREAM')

    do i = 1, N
        call random_number(x)
        call random_number(y)
        call random_number(z)

        x = (x * len - len/2)
        y = (y * len - len/2)
        z = z*d

        step = i
        species = species_elec

        write(unit=ud_file) x, y, z, step, species
    end do

    close(unit=ud_file, status='keep')

    !------------------------------------------------------------------------------------------------------------
end program Write_Unit_Test_Files