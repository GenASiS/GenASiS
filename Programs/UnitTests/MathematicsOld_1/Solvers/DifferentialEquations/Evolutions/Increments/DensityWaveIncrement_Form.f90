module DensityWaveIncrement_Form

  use Basics
  use Manifolds
  use Fields
  use ProtoFields
  use StorageDivergence_Form
  use IncrementDivergence_FV__Form

  implicit none
  private

  type, public, extends ( DensityWaveForm ) :: DensityWaveIncrementForm
    integer ( KDI ) :: &
      nRampCycles = 1
    type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
      BoundaryFluence_CSL
    type ( StorageForm ), allocatable :: &
      Increment_S
    type ( StorageDivergenceForm ), allocatable :: &
      Storage
    type ( IncrementDivergence_FV_Form ), allocatable :: &
      Increment
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, public, pass :: &
      PrepareBoundaryFluence
    procedure, public, pass :: &
      ComputeTimeStep
    procedure, private, pass :: &
      ComputeTimeStepLocal
  end type DensityWaveIncrementForm

    private :: &
      ComputeTimeStepKernel_CSL

contains


  subroutine Initialize ( DW, Name )

    class ( DensityWaveIncrementForm ), intent ( inout ) :: &
      DW
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iF  !-- iField
    real ( KDR ) :: &
      TimeStep
    character ( LDL ), dimension ( : ), allocatable :: &
      ConservedVariable
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( ProtoCurrentForm ), pointer :: &
      PC

    !-- Initialize parent

    call DW % DensityWaveForm % Initialize ( Name )

    !-- Initialize members and values needed for Increment

    associate ( PCA => DW % ProtoCurrent )
    PC => PCA % ProtoCurrent ( )
    end associate !-- PCA

    select type ( C => DW % Atlas % Chart )
    class is ( Chart_SLD_Form )
    G => C % Geometry ( )

    call DW % PrepareBoundaryFluence ( C, PC )

    allocate ( ConservedVariable ( PC % N_CONSERVED ) )
    ConservedVariable &
      = [ ( PC % Variable ( PC % iaConserved ( iF ) ), &
            iF = 1, PC % N_CONSERVED ) ] 

    allocate ( DW % Increment_S )
    associate ( S_Increment => DW % Increment_S )
    call S_Increment % Initialize &
           ( [ PC % nValues, PC % N_CONSERVED ], &
             VariableOption = ConservedVariable, &
             ClearOption = .true., NameOption = 'Increment' )
    
    if ( PC % AllocatedDevice ) &
      call S_Increment % AllocateDevice ( )

    call C % AddFieldImage ( S_Increment, iStream = 1 )

    !-- TimeStep

    call DW % ComputeTimeStep ( TimeStep )
    call Show ( TimeStep, 'TimeStep', DW % IGNORABILITY )

    !-- Increment

    allocate ( DW % Storage )
    allocate ( DW % Increment )
    associate &
      ( I => DW % Increment, &
        S => DW % Storage, &
        Weight_RK => 0.5_KDR )

    call S % Allocate &
           ( AllocateDevice = .true., nCurrent = PC % nVariables, &
             nConserved = PC % N_CONSERVED, &
             nReconstructed = PC % N_PRIMITIVE, &
             nSolverSpeeds = PC % N_SOLVER_SPEEDS, &
             nGeometry = G % nVariables, nValues = PC % nValues )

    call I % Initialize ( DW % ProtoCurrent % Chart )
    call I % SetStorage ( DW % Storage )
    call I % SetBoundaryFluence ( DW % BoundaryFluence_CSL )

    call I % Compute ( S_Increment, TimeStep, Weight_RK )

    associate ( PCS => PC % Sources )
    PCS % Value ( :, 1 : PCS % N_FIELDS_C )  &
      =  Weight_RK  *  S_Increment % Value  /  TimeStep
    end associate !-- PCS

    end associate !-- S_Increment
    end associate !-- I, etc.
    end select !-- C
    nullify ( PC, G )

  end subroutine Initialize


  subroutine Finalize ( DW )

    type ( DensityWaveIncrementForm ), intent ( inout ) :: &
      DW

    if ( allocated ( DW % Increment ) ) &
      deallocate ( DW % Increment )
    if ( allocated ( DW % Storage ) ) &
      deallocate ( DW % Storage )
    if ( allocated ( DW % Increment_S ) ) &
      deallocate ( DW % Increment_S )
    if ( allocated ( DW % BoundaryFluence_CSL ) ) &
      deallocate ( DW % BoundaryFluence_CSL )

  end subroutine Finalize


  subroutine PrepareBoundaryFluence ( DW, Grid, Current )

    class ( DensityWaveIncrementForm ), intent ( inout ) :: &
      DW
    class ( * ), intent ( in ) :: &
      Grid
    class ( CurrentTemplate ), intent ( in ) :: &
      Current

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension
      iF      !-- iField
    integer ( KDI ), dimension ( 3 ) :: &
      nSurface

    select type ( Grid )
    class is ( Chart_SLD_Form )

      if ( allocated ( DW % BoundaryFluence_CSL ) ) &
        deallocate ( DW % BoundaryFluence_CSL )

      associate &
        ( C => Grid % Atlas % Connectivity, &
          nDimensions => Grid % nDimensions )
      allocate &
        ( DW % BoundaryFluence_CSL ( Current % N_CONSERVED, C % nFaces ) )
      do iD = 1, nDimensions
        jD = mod ( iD, 3 ) + 1
        kD = mod ( jD, 3 ) + 1
        nSurface ( iD ) = 1
        nSurface ( jD ) = Grid % nCellsBrick ( jD ) 
        nSurface ( kD ) = Grid % nCellsBrick ( kD )
        do iF = 1, Current % N_CONSERVED
          call DW % BoundaryFluence_CSL ( iF, C % iaInner ( iD ) ) &
                 % Initialize ( nSurface, ClearOption = .true. )
          call DW % BoundaryFluence_CSL ( iF, C % iaOuter ( iD ) ) &
                 % Initialize ( nSurface, ClearOption = .true. )
        end do !-- iF
      end do !-- iD
      end associate !-- C, etc.

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'DensityWaveIncrement_Form', 'module', CONSOLE % ERROR )
      call Show ( 'PrepareBoundaryFluence', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

  end subroutine PrepareBoundaryFluence


  subroutine ComputeTimeStep ( DW, TimeStep )

    class ( DensityWaveIncrementForm ), intent ( in ) :: &
      DW
    real ( KDR ), intent ( out ) :: &
      TimeStep

    real ( KDR ) :: &
      RampFactor
    type ( CollectiveOperation_R_Form ) :: &
      CO

    TimeStep = huge ( 0.0_KDR )

    call DW % ComputeTimeStepLocal ( TimeStep )

    call CO % Initialize &
           ( DW % Atlas % Communicator, nOutgoing = [ 1 ], nIncoming = [ 1 ] )
    CO % Outgoing % Value ( 1 ) = TimeStep

    call CO % Reduce ( REDUCTION % MIN )

    RampFactor &
      = min ( real ( DW % iCycle + 1, KDR ) / DW % nRampCycles, 1.0_KDR )
    TimeStep &
      = RampFactor * CO % Incoming % Value ( 1 )

  end subroutine ComputeTimeStep


  subroutine ComputeTimeStepLocal ( DW, TimeStep )

    class ( DensityWaveIncrementForm ), intent ( in ) :: &
      DW
    real ( KDR ), intent ( inout ) :: &
      TimeStep

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( ProtoCurrentForm ), pointer :: &
      PC

    select type ( CSL => DW % Atlas % Chart )
    class is ( Chart_SL_Template )

    associate ( PCA => DW % ProtoCurrent )

    G  => CSL % Geometry ( )
    PC => PCA % ProtoCurrent ( )
    
    call ComputeTimeStepKernel_CSL &
           ( CSL % IsProperCell, &
             PC % Value ( :, PC % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             PC % Value ( :, PC % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             PC % Value ( :, PC % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             PC % Value ( :, PC % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             PC % Value ( :, PC % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             PC % Value ( :, PC % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
             G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), & 
             G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
             G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
             G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), & 
             G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
             CSL % nDimensions, TimeStep )

    end associate !-- PCA
    end select !-- CSL

    nullify ( PC, G )

  end subroutine ComputeTimeStepLocal


  subroutine ComputeTimeStepKernel_CSL &
               ( IsProperCell, FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                 dX_L_1, dX_L_2, dX_L_3, dX_R_1, dX_R_2, dX_R_3, &
                 nDimensions, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      dX_L_1, dX_L_2, dX_L_3, &
      dX_R_1, dX_R_2, dX_R_3
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    real ( KDR ), intent ( inout ) :: &
      TimeStep

    real ( KDR ) :: &
      TimeStepInverse

    select case ( nDimensions )
    case ( 1 )
      TimeStepInverse &
        = maxval ( max ( FEP_1, -FEM_1 ) / ( dX_L_1 + dX_R_1 ), &
                   mask = IsProperCell )
    case ( 2 )
      TimeStepInverse &
        = maxval (   max ( FEP_1, -FEM_1 ) / ( dX_L_1 + dX_R_1 ) &
                   + max ( FEP_2, -FEM_2 ) / ( dX_L_2 + dX_R_2 ), &
                   mask = IsProperCell )
    case ( 3 )
      TimeStepInverse &
        = maxval (   max ( FEP_1, -FEM_1 ) / ( dX_L_1 + dX_R_1 ) &
                   + max ( FEP_2, -FEM_2 ) / ( dX_L_2 + dX_R_2 ) &
                   + max ( FEP_3, -FEM_3 ) / ( dX_L_3 + dX_R_3 ), &
                   mask = IsProperCell )
    end select !-- nDimensions

    TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStepKernel_CSL


end module DensityWaveIncrement_Form
