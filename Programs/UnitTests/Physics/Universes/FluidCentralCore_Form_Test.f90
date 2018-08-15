program FluidCentralCore_Form_Test

  use Basics
  use FluidCentralCore_Form

  implicit none

  type ( FluidCentralCoreForm ), allocatable :: &
    FCC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FluidCentralCore_Form_Test', DimensionalityOption = '2D' )

  allocate ( FCC )
  call FCC % Initialize &
         ( PROGRAM_HEADER % Name, &
           FluidType = 'DUST', &
           GeometryType = 'NEWTONIAN', &
           DimensionlessOption = .true. )

  call FCC % OpenManifoldStreams ( VerboseStreamOption = .true. )
  call FCC % Write ( )

  call FCC % CoarsenSingularities ( )

  deallocate ( FCC )

  deallocate ( PROGRAM_HEADER )

end program FluidCentralCore_Form_Test
