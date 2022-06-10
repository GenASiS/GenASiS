program Chart_H__Form_Test

  !-- Chart_Header__Form_Test

  use Basics
  use StructuredGrids

  implicit none

  type ( Chart_H_Form ), allocatable :: &
    C_Base, &
    C_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Chart_H__Form_Test', DimensionalityOption = '2D_1D' )

  allocate ( C_Base )
  call C_Base % Initialize_H &
         ( NameOption = 'C_Base', &
           iDimensionalityOption = 1 )

  allocate ( C_Fiber )
  call C_Fiber % Initialize_H &
         ( CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL', &
           NameOption = 'C_Fiber', &
           iDimensionalityOption = 2 )

  call C_Base % Show ( )
  call C_Fiber % Show ( )

  deallocate ( C_Fiber )
  deallocate ( C_Base )
  deallocate ( PROGRAM_HEADER )

end program Chart_H__Form_Test
