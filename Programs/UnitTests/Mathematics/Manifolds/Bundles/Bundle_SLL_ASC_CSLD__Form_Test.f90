module Bundle_SLL_ASC_CSLD__Form_Test__Form

  use Basics
  use Atlases
  use Bundle_SLL_ASC_CSLD__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Bundle_SLL_ASC_CSLD__Form_Test'

  type, public :: Bundle_SLL_ASC_CSLD_Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( Atlas_SC_Form ), allocatable :: &
      Base
    type ( Bundle_SLL_ASC_CSLD_Form ), allocatable :: &
      Bundle
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Bundle_SLL_ASC_CSLD_Form_Test_Form

contains

  
  subroutine Initialize ( BFT, Name )

    class ( Bundle_SLL_ASC_CSLD_Form_Test_Form ), intent ( inout ), &
      target :: &
        BFT
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ), dimension ( 3 ) :: &
      Scale
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    !-- Base 

    allocate &
      ( BFT % Base, &
        BFT % GridImageStream )
    associate &
      ( Base => BFT % Base, &
        GIS  => BFT % GridImageStream )

    call Base % Initialize &
           ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
             iDimensionalityOption = 1 )
    call Base % CreateChart ( nCellsOption = [ 16, 16, 16 ] )  
    call Base % SetGeometry ( )

    !-- Bundle

    allocate ( BFT % Bundle )
    associate ( Bundle => BFT % Bundle )

    call Bundle % Initialize ( Base, 'Fiber' )
    call Bundle % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    Scale = 0.0_KDR
    Scale ( 1 ) = 5.0_KDR * UNIT % MEGA_ELECTRON_VOLT

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % MEGA_ELECTRON_VOLT

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    call Bundle % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateLabelOption = [ 'Energy  ', 'VarTheta', 'VarPhi  ' ], &
             CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = Scale, nCellsOption = [ 8, 8, 8 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !--- Open stream

    call GIS % Initialize &
           ( PROGRAM_HEADER % Name, CommunicatorOption = Base % Communicator )
    call Base % OpenStream ( GIS, '1', iStream = 1 )
    call Bundle % OpenStream ( GIS, iStream = 1 )

    !-- Write

    call GIS % Open ( GIS % ACCESS_CREATE )
    call Bundle % MarkFibersWritten ( )
    call Base % Write ( iStream = 1 )
    call GIS % Close ( )

    !-- Base's GIS must be closed before call to Bundle % Write ( ).
    call Bundle % Write ( iStream = 1 )

    end associate !-- B, etc.
    end associate !-- Base, etc.

  end subroutine Initialize


  subroutine Finalize ( BFT )

    type ( Bundle_SLL_ASC_CSLD_Form_Test_Form ), intent ( inout ) :: &
      BFT

    deallocate ( BFT % Bundle )
    deallocate ( BFT % Base )
    deallocate ( BFT % GridImageStream )

  end subroutine Finalize


end module Bundle_SLL_ASC_CSLD__Form_Test__Form



program Bundle_SLL_ASC_CSLD__Form_Test

  use Basics
  use Bundle_SLL_ASC_CSLD__Form_Test__Form
  
  implicit none

  type ( Bundle_SLL_ASC_CSLD_Form_Test_Form ), allocatable :: &
    BFT
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, DimensionalityOption = '2D_1D' )

  allocate ( BFT )
  call BFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( BFT )

  deallocate ( PROGRAM_HEADER )

end program Bundle_SLL_ASC_CSLD__Form_Test
