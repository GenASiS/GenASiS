module DensityWaveIncrement_Form

  use Basics
  use Manifolds
  use Fields
  use ProtoFields
  use IncrementDivergence_FV__Form

  implicit none
  private

  type, public, extends ( DensityWaveForm ) :: DensityWaveIncrementForm
    integer ( KDI ) :: &
      nRampCycles = 1
    type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
      BoundaryFluence_CSL
    type ( VariableGroupForm ), allocatable :: &
      Increment_VG
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

    call DW % PrepareBoundaryFluence ( C, PC )

    allocate ( ConservedVariable ( PC % N_CONSERVED ) )
    ConservedVariable &
      = [ ( PC % Variable ( PC % iaConserved ( iF ) ), &
            iF = 1, PC % N_CONSERVED ) ] 

    allocate ( DW % Increment_VG )
    associate ( VG_Increment => DW % Increment_VG )
    call VG_Increment % Initialize &
           ( [ PC % nValues, PC % N_CONSERVED ], &
             VariableOption = ConservedVariable, NameOption = 'Increment' )

    call C % AddFieldImage ( VG_Increment, iStream = 1 )

    !-- TimeStep

    call DW % ComputeTimeStep ( TimeStep )
    call Show ( TimeStep, 'TimeStep', DW % IGNORABILITY )

    !-- Increment

    allocate ( DW % Increment )
    associate ( I => DW % Increment )

    call I % Initialize ( Name )
    call I % Set ( DW % BoundaryFluence_CSL )
    call I % Set ( Weight_RK = 0.5_KDR )

    call I % Compute ( VG_Increment, C, PC, TimeStep )

    end associate !-- VG_Increment
    end associate !-- I
    end select !-- C
    nullify ( PC )

  end subroutine Initialize


  subroutine Finalize ( DW )

    type ( DensityWaveIncrementForm ), intent ( inout ) :: &
      DW

    if ( allocated ( DW % Increment ) ) &
      deallocate ( DW % Increment )
    if ( allocated ( DW % Increment_VG ) ) &
      deallocate ( DW % Increment_VG )
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
             G % Value ( :, G % WIDTH ( 1 ) ), &
             G % Value ( :, G % WIDTH ( 2 ) ), & 
             G % Value ( :, G % WIDTH ( 3 ) ), &
             CSL % nDimensions, TimeStep )

    end associate !-- PCA
    end select !-- CSL

    nullify ( PC, G )

  end subroutine ComputeTimeStepLocal


  subroutine ComputeTimeStepKernel_CSL &
               ( IsProperCell, FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                 dX_1, dX_2, dX_3, nDimensions, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      dX_1, dX_2, dX_3
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    real ( KDR ), intent ( inout ) :: &
      TimeStep

    real ( KDR ) :: &
      TimeStepInverse

    select case ( nDimensions )
    case ( 1 )
      TimeStepInverse &
        = maxval ( max ( FEP_1, -FEM_1 ) / dX_1, mask = IsProperCell )
    case ( 2 )
      TimeStepInverse &
        = maxval (   max ( FEP_1, -FEM_1 ) / dX_1 &
                   + max ( FEP_2, -FEM_2 ) / dX_2, &
                   mask = IsProperCell )
    case ( 3 )
      TimeStepInverse &
        = maxval (   max ( FEP_1, -FEM_1 ) / dX_1 &
                   + max ( FEP_2, -FEM_2 ) / dX_2 &
                   + max ( FEP_3, -FEM_3 ) / dX_3, &
                   mask = IsProperCell )
    end select !-- nDimensions

    TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStepKernel_CSL


end module DensityWaveIncrement_Form
