program Geometry_F_CH__Form_Test

  !-- Geometry_Flat_ChartHeader__Form_Test

  use Basics
  use Manifolds
  use Streams
  use Geometries

  implicit none

  ! type ( MeasuredValueForm ), dimension ( 5 ) :: &
  !   FieldUnit
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( Stream_CH_Form ), allocatable :: &
    SC
  type ( Geometry_F_CH_Form ), allocatable :: &
    GC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Geometry_F_CH__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize_H ( PeriodicOption = [ .true., .true., .true. ] )

  ! FieldUnit ( 1 )      =  UNIT % MASS_DENSITY_MKS
  ! FieldUnit ( 2 : 4 )  =  UNIT % SPEED_MKS
  ! FieldUnit ( 5 )      =  UNIT % JOULE

  allocate ( GC )
  call GC % Initialize_H &
         ( C ) !UnitOption = FieldUnit )

  allocate ( SC )
  call SC % Initialize_H ( C, GIS )

  ! call CONSOLE % SetVerbosity ( 'INFO_2' )
  ! call SC % AddFieldSet ( FSC )
  ! call SC % AddFieldSet ( FSC_234 )
  ! call SC % AddFieldSet ( FSC_5 )
  ! call CONSOLE % SetVerbosity ( 'INFO_1' )

  call   C % Show ( )
  ! call FSC     % Show ( )
  ! call FSC_234 % Show ( )
  ! call FSC_5   % Show ( )
  call  SC % Show ( )

  deallocate ( GC )
  deallocate ( SC )
  ! deallocate ( FSC_5 )
  ! deallocate ( FSC_234 )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Geometry_F_CH__Form_Test
