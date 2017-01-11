module ChartHeader_SL__Form_Test__Form

  use Basics
  use AtlasBasics
  use ChartHeader_SL__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'ChartHeader_SL__Form_Test'

  type, public :: ChartHeader_SL__Form_Test_Form
    type ( AtlasHeaderForm ), allocatable :: &
      AtlasHeaderBase, &
      AtlasHeaderFiber
    type ( ChartHeader_SL_Form ), allocatable :: &
      ChartHeaderBase, &
      ChartHeaderFiber
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type ChartHeader_SL__Form_Test_Form

contains

  
  subroutine Initialize ( CHFT )

    class ( ChartHeader_SL__Form_Test_Form ), intent ( inout ), target :: &
      CHFT

    real ( KDR ), dimension ( 3 ) :: &
      Scale
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    logical ( KDL ), dimension ( 3 ) :: &
      IsPeriodic
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    !-- Base Atlas 

    allocate &
      ( CHFT % AtlasHeaderBase, &
        CHFT % ChartHeaderBase )
    associate &
      ( A => CHFT % AtlasHeaderBase, &
        C => CHFT % ChartHeaderBase )

    call A % Initialize &
           ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
             iDimensionalityOption = 1 )

    IsPeriodic = .true.

    call C % Initialize ( A, IsPeriodic, iChart = 1 )

    call A % Show ( )
    call C % Show ( )
    end associate !-- AH, etc.

    !-- Fiber Atlas

    allocate &
      ( CHFT % AtlasHeaderFiber, &
        CHFT % ChartHeaderFiber )
    associate &
      ( A => CHFT % AtlasHeaderFiber, &
        C => CHFT % ChartHeaderFiber )

    call A % Initialize ( 'Fiber', iDimensionalityOption = 2 )
    call A % SetBoundaryConditionsFace ( [ 'REFLECTING', 'REFLECTING' ], 1 )

    Scale = 0.0_KDR
    Scale ( 1 ) = 5.0_KDR * UNIT % MEV

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % MEV

    IsPeriodic = .false.

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    call C % Initialize &
           ( A, IsPeriodic, iChart = 1, SpacingOption = Spacing, &
             CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = Scale )

    call A % Show ( )
    call C % Show ( )
    end associate !-- AH, etc.

  end subroutine Initialize


  subroutine Finalize ( CHFT )

    type ( ChartHeader_SL__Form_Test_Form ), intent ( inout ) :: &
      CHFT

    deallocate ( CHFT % ChartHeaderFiber )
    deallocate ( CHFT % AtlasHeaderFiber )

    deallocate ( CHFT % ChartHeaderBase )
    deallocate ( CHFT % AtlasHeaderBase )

  end subroutine Finalize


end module ChartHeader_SL__Form_Test__Form



program ChartHeader_SL__Form_Test

  use Basics
  use ChartHeader_SL__Form_Test__Form
  
  implicit none

  type ( ChartHeader_SL__Form_Test_Form ), allocatable :: &
    CHFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, DimensionalityOption = '3D_1D' )
    
  allocate ( CHFT )
  call CHFT % Initialize ( )
  deallocate ( CHFT )

  deallocate ( PROGRAM_HEADER )

end program ChartHeader_SL__Form_Test
