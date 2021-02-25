program AtlasHeader_Form_Test

  use Basics
  use AtlasBasics

  type ( AtlasHeaderForm ), allocatable :: &
    AH_Base, &
    AH_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'AtlasHeader_Form_Test', DimensionalityOption = '3D_1D' )

  allocate ( AH_Base )
  associate ( A => AH_Base )
  call A % Initialize &
         ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           iDimensionalityOption = 1 )
  call A % Show ( )
  end associate !-- A

  allocate ( AH_Fiber )
  associate ( A => AH_Fiber )
  call A % Initialize ( 'Fiber', iDimensionalityOption = 2 )
  call A % Show ( )
  end associate !-- A

  deallocate ( AH_Fiber )
  deallocate ( AH_Base )
  deallocate ( PROGRAM_HEADER )

end program AtlasHeader_Form_Test
