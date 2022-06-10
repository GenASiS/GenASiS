program Geometry_CSL__Form_Test

  use Basics
  use Mathematics
  use Geometry_CSL__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Geometry_CSL__Form_Test'

  type ( Atlas_SC_Form ), allocatable :: &
    A_G, A_N, A_N_S
  type ( Geometry_CSL_Form ), allocatable :: &
    GC_G, GC_N, GC_N_S

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( A_G, A_N, A_N_S )
  call A_G % Initialize ( 'Atlas_SC_G', PROGRAM_HEADER % Communicator )
  call A_N % Initialize ( 'Atlas_SC_N', PROGRAM_HEADER % Communicator )
  call A_N_S % Initialize ( 'Atlas_SC_N_S', PROGRAM_HEADER % Communicator )
  call A_G % CreateChart ( )
  call A_N % CreateChart ( )
  call A_N_S % CreateChart ( )

  select type ( C_G => A_G % Chart )
  class is ( Chart_SL_Template )
  select type ( C_N => A_N % Chart )
  class is ( Chart_SL_Template )
  select type ( C_N_S => A_N_S % Chart )
  class is ( Chart_SL_Template )

  associate ( nValues => C_G % nProperCells + C_G % nGhostCells )

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  allocate ( GC_G, GC_N, GC_N_S )
  call GC_G % Initialize ( C_G, 'Geometry_G', 'GALILEAN', nValues )
  call GC_N % Initialize ( C_N, 'Geometry_N', 'NEWTONIAN', nValues )
  call GC_N_S % Initialize &
         ( C_N_S, 'Geometry_N_S', 'NEWTONIAN_STRESS', nValues )
  deallocate ( GC_N_S, GC_N, GC_G )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  end associate !-- nValues
  end select !-- C_N_S
  end select !-- C_N
  end select !-- C_G
  deallocate ( A_N_S, A_N, A_G )
  deallocate ( PROGRAM_HEADER )

end program Geometry_CSL__Form_Test
