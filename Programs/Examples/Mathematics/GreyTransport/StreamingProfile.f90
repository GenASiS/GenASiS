program StreamingProfile

  use Basics
  use StreamingProfile_Form

  implicit none

  type ( StreamingProfileForm ), allocatable :: &
    SP

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'StreamingProfile', DimensionalityOption = '1D' )

  allocate ( SP )
  call SP % Initialize ( PROGRAM_HEADER % Name )
  call SP % Evolve ( )
  deallocate ( SP)

  deallocate ( PROGRAM_HEADER )

end program StreamingProfile
