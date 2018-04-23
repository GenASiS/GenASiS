program Geometry_CSL__Form_Test

  use Basics
  use Mathematics
  use Geometry_CSL__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Geometry_CSL__Form_Test'

  type ( Atlas_SC_Form ), allocatable :: &
    A_G, A_N
  type ( Geometry_CSL_Form ), allocatable :: &
    GC_G, GC_N

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( A_G, A_N )
  call A_G % Initialize ( 'Atlas_SC_G', PROGRAM_HEADER % Communicator )
  call A_N % Initialize ( 'Atlas_SC_N', PROGRAM_HEADER % Communicator )
  call A_G % CreateChart ( )
  call A_N % CreateChart ( )

  select type ( C_G => A_G % Chart )
  class is ( Chart_SL_Template )
  select type ( C_N => A_N % Chart )
  class is ( Chart_SL_Template )

  associate ( nValues => C_G % nProperCells + C_G % nGhostCells )

  allocate ( GC_G, GC_N )
  call GC_G % Initialize ( C_G, 'Geometry_CSL_G', 'GALILEAN', nValues )
  call GC_N % Initialize ( C_N, 'Geometry_CSL_N', 'NEWTONIAN', nValues )

  deallocate ( GC_N, GC_G )
  end associate !-- nValues
  end select !-- C_N
  end select !-- C_G
  deallocate ( A_N, A_G )
  deallocate ( PROGRAM_HEADER )

end program Geometry_CSL__Form_Test
