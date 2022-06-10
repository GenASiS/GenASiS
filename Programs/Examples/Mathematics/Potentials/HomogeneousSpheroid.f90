program HomogeneousSpheroid

  use Basics
  use Mathematics
  use HomogeneousSpheroid_Form

  implicit none

  type ( HomogeneousSpheroidForm ), allocatable :: &
    HS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'HomogeneousSpheroid', DimensionalityOption = '2D' )

  if ( trim ( PROGRAM_HEADER % Dimensionality )  ==  '1D' ) then
    call Show ( 'This program must be run in 2D or 3D', CONSOLE % ERROR )
    call Show ( 'HomogeneousSpheroid', 'program', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )
  end if

  allocate ( HS )
  call HS % Initialize ( )
  call HS % Compute ( )
  deallocate ( HS )

  deallocate ( PROGRAM_HEADER )

end program HomogeneousSpheroid
