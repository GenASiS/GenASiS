program Fluid_ASC__Form_Test

  use Basics
  use Mathematics
  use Spaces
  use Fluid_ASC__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Fluid_ASC__Form_Test'

  type ( Atlas_SC_Form ), allocatable :: &
    A
  type ( Geometry_ASC_Form ), allocatable :: &
    GA
  type ( Fluid_ASC_Form ), allocatable :: &
    FA_D, FA_I

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( A )
  call A % Initialize ( 'Atlas_SC', PROGRAM_HEADER % Communicator )
  call A % CreateChart ( )

  allocate ( GA )
  call GA % Initialize ( A, 'GALILEAN' )
  call A % SetGeometry ( GA )

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  allocate ( FA_D, FA_I )
  call FA_D % Initialize ( A, 'DUST', NameShortOption = 'Fluid_D' )
  call FA_I % Initialize ( A, 'IDEAL', NameShortOption = 'Fluid_I' )
  deallocate ( FA_I, FA_D )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  deallocate ( GA )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program Fluid_ASC__Form_Test
