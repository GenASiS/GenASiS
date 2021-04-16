program Stream_CH__Form_Test

  !-- Stream_ChartHeader__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none

  integer ( KDI ) :: &
    nFields
  type ( Integer_1D_Form ), dimension ( 1 ) :: &
    VectorIndices
  type ( MeasuredValueForm ), dimension ( 5 ) :: &
    FieldUnit
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( FieldSet_CH_Form ), allocatable :: &
    FSC, &
    FSC_234, &
    FSC_5
  type ( Stream_CH_Form ), allocatable :: &
    SC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_CH__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

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
  call FSC_234 % Clone &
         ( FSC, NameOption = 'Fields_234', iaSelectedOption = [ 2, 3, 4 ] )
  call FSC_5 % Clone &
         ( FSC, NameOption = 'Fields_5', iaSelectedOption = [ 5 ] )

  allocate ( SC )
  call SC % Initialize_H ( C, GIS )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call SC % AddFieldSet ( FSC )
  call SC % AddFieldSet ( FSC_234 )
  call SC % AddFieldSet ( FSC_5 )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call   C     % Show ( )
  call FSC     % Show ( )
  call FSC_234 % Show ( )
  call FSC_5   % Show ( )
  call  SC     % Show ( )

  deallocate ( SC )
  deallocate ( FSC_5 )
  deallocate ( FSC_234 )
  deallocate ( FSC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Stream_CH__Form_Test
