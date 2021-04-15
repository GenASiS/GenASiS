program FieldSet_CH__Form_Test

  !-- FieldSet_ChartHeader__Form_Test

  use Basics
  use Manifolds
  use FieldSets

  implicit none

  integer ( KDI ) :: &
    nFields
  type ( Integer_1D_Form ), dimension ( 1 ) :: &
    VectorIndices
  type ( MeasuredValueForm ), dimension ( 5 ) :: &
    FieldUnit
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( FieldSet_CH_Form ), allocatable :: &
    FSC, &
    FSC_234, &
    FSC_5

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_CH__Form_Test', DimensionalityOption = '2D' )

  allocate ( C )
  call C % Initialize_H ( PeriodicOption = [ .true., .true., .true. ] )

  nFields  =  5

  FieldUnit ( 1 )      =  UNIT % MASS_DENSITY_MKS
  FieldUnit ( 2 : 4 )  =  UNIT % SPEED_MKS
  FieldUnit ( 5 )      =  UNIT % JOULE

  call VectorIndices ( 1 ) % Initialize ( [ 2, 3, 4 ] )

  allocate ( FSC )
  allocate ( FSC_234 )
  allocate ( FSC_5 )
  call FSC % Initialize_H &
         ( C, &
           UnitOption = FieldUnit, &
           VectorIndicesOption = VectorIndices, &
           nFieldsOption = nFields )
  call FSC_234 % Initialize_H &
         ( FSC, NameOption = 'Fields_234', iaSelectedOption = [ 2, 3, 4 ] )
  call FSC_5 % Initialize_H &
         ( FSC, NameOption = 'Fields_5', iaSelectedOption = [ 5 ] )

  call   C     % Show ( )
  call FSC     % Show ( )
  call FSC_234 % Show ( )
  call FSC_5   % Show ( )

  deallocate ( FSC_5 )
  deallocate ( FSC_234 )
  deallocate ( FSC )
  deallocate ( C )
  deallocate ( PROGRAM_HEADER )

end program FieldSet_CH__Form_Test
