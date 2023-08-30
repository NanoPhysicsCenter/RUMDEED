!
! Kristinn Torfason
! Vacuum electronics molecular dynamics simulations
! Units test program
program Vacuum_Unit_Test
! Vacuum
  use iso_fortran_env
#if defined(_OPENMP)
  use omp_lib
#endif
  use mod_global
  use mod_verlet
  use mod_photo_emission
  use mod_field_emission
  use mod_field_emission_v2
  use mod_emission_tip
  use mod_field_emission_2D
  use mod_field_thermo_emission
  use mod_pair
  use mod_unit_tests
#if defined(__INTEL_COMPILER)
  use IFPORT ! Needed for getpid()
#endif
  implicit none
#if defined(__PGI)
  interface
    integer function getpid()
    end function getpid
  end interface
#endif

  print *, 'Unit tests'
end program Vacuum_Unit_Test
