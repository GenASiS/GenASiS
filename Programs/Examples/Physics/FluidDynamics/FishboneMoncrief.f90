program FishboneMoncrief

  use GenASiS
  use FishboneMoncrief_Form

  implicit none

  type ( FishboneMoncriefForm ), allocatable :: &
    FM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FishboneMoncrief', DimensionalityOption = '2D' )

  allocate ( FM )
  call FM % Initialize ( PROGRAM_HEADER % Name )
  call FM % Evolve ( )
  deallocate ( FM )

  deallocate ( PROGRAM_HEADER )

end program FishboneMoncrief
