program ChartHeader_Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'ChartHeader_Form_Test', DimensionalityOption = '2D' )

  deallocate ( PROGRAM_HEADER )

end program ChartHeader_Form_Test
