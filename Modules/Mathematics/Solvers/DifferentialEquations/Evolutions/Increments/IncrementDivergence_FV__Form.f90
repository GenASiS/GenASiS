!-- IncrementDivergence_FV computes a first-order (in time) update due the
!   divergence in a conservation law.

module IncrementDivergence_FV__Form

  !-- IncrementDivergence_FiniteVolume_Form
  
  use Basics
  use Manifolds
  use Operations
  use Fields
  use StorageDivergence_Form

  implicit none
  private

  type, public :: IncrementDivergence_FV_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
!       iTimerReconstruction_G, &
!       iTimerBoundary, &
!       iTimerReconstruction_CSL, &
!       iTimerFromPrimitive, &
! !       iTimerGradient, &
! !       iTimerReconstructionKernel, &
      iStream
    real ( KDR ) :: &
      Weight_RK = - huge ( 0.0_KDR ) !-- RungeKutta weight
    real ( KDR ), dimension ( : ), pointer :: &
      UseLimiter => null ( )
    type ( Real_1D_Form ), dimension ( : ), pointer :: &
      dLogVolumeJacobian_dX => null ( )
    type ( Real_3D_Form ), dimension ( :, : ), pointer :: &
      BoundaryFluence_CSL => null ( )
    logical ( KDL ) :: &
      UseIncrementStream
    character ( LDL ) :: &
      ReconstructionType
    character ( LDF ) :: &
      Name = ''
    type ( StorageForm ), dimension ( : ), allocatable :: &
      Output
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    class ( ChartTemplate ), pointer :: &
      Chart => null ( )
    class ( FieldChartTemplate ), pointer :: &
      CurrentChart => null ( )
    class ( GeometryFlatForm ), pointer :: &
      Geometry => null ( )
    class ( CurrentTemplate ), pointer :: &
      Current => null ( )
    class ( StorageDivergenceForm ), pointer :: &
      Storage => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SetStorage
    procedure, public, pass :: &
      ClearStorage
    procedure, private, pass :: &
      SetBoundaryFluence_CSL
    generic, public :: &
      SetBoundaryFluence => SetBoundaryFluence_CSL
    procedure, public, pass :: &
      ClearBoundaryFluence
    procedure, private, pass :: &
      SetMetricDerivativesFlat
    generic, public :: &
      SetMetricDerivatives => SetMetricDerivativesFlat
    procedure, public, pass :: &
      SetUseLimiter
    procedure, public, pass :: &
      ClearMetricDerivatives
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type IncrementDivergence_FV_Form

    private :: &
      SetCurrent, &
      PrepareOutput, &
      WriteOutput, &
      ComputeReconstruction, &
      ComputeFluxes

      private :: &
        ComputeReconstructionLinear_CSL, &
        ComputeReconstructionParabolic_CSL, &
        ComputeIncrement_CSL

        private :: &
          ComputeReconstructionLinear_CSL_Kernel, &
          ComputeReconstructionParabolic_CSL_Kernel, &
          ComputeLogDerivative_CSL_Kernel, &
          ComputeIncrement_CSL_Kernel, &
          RecordBoundaryFluence_CSL

          private :: &
            RecordBoundaryFluence_CSL_Kernel

  integer ( KDI ), private :: &
    iCURRENT    = 1, &
    iCURRENT_IL = 2, &
    iCURRENT_IR = 3, &
    iFLUX_IL    = 4, &
    iFLUX_IR    = 5, &
    iFLUX_I     = 6, &
    iINCREMENT  = 7, &
    nOUTPUT = 7
    
  interface
  
    module subroutine ComputeReconstructionLinear_CSL_Kernel &
                 ( V, dVdX, dX_L, dX_R, iD, oV, V_IL, V_IR, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        V, &
        dVdX, &
        dX_L, dX_R
      integer ( KDI ), intent ( in ) :: &
        iD, &
        oV   
      real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
        V_IL, V_IR
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeReconstructionLinear_CSL_Kernel

    module subroutine ComputeReconstructionParabolic_CSL_Kernel &
                 ( V, X, dX_L, dX_R, XA, X2A, iD, oV, V_IL, V_IR, &
                   UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        V, &
        X, &
        dX_L, dX_R, &
        XA, X2A
      integer ( KDI ), intent ( in ) :: &
        iD, &
        oV   
      real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
        V_IL, V_IR
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeReconstructionParabolic_CSL_Kernel

    module subroutine ComputeLogDerivative_CSL_Kernel &
                 ( A_I, V, iD, oV, dLVdX, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        A_I, &
        V
      integer ( KDI ), intent ( in ) :: &
        iD, &
        oV   
      real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
        dLVdX
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeLogDerivative_CSL_Kernel
    
    module subroutine ComputeIncrement_CSL_Kernel &
                 ( dU, F_I, A_I, V, dT, iD, oV, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
        dU
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        F_I, &
        A_I, &
        V
      real ( KDR ), intent ( in ) :: &
        dT
      integer ( KDI ), intent ( in ) :: &
        iD, &
        oV
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeIncrement_CSL_Kernel
    
    module subroutine RecordBoundaryFluence_CSL_Kernel &
                 ( BF, F, Factor, nB, oB, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
        BF
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        F
      real ( KDR ), intent ( in ) :: &
        Factor
      integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
        nB, &
        oB
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine RecordBoundaryFluence_CSL_Kernel
      
  end interface

contains


  subroutine Initialize ( I, CurrentChart )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    class ( FieldChartTemplate ), intent ( in ), target :: &
      CurrentChart

    character ( LDF ) :: &
      OutputDirectory

    I % IGNORABILITY = CurrentChart % IGNORABILITY
    I % Name = 'IncrementDivergence_' // trim ( CurrentChart % Name )

    call Show ( 'Initializing an IncrementDivergence_FV', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    I % CurrentChart => CurrentChart
    I % Chart        => CurrentChart % Chart

    ! call PROGRAM_HEADER % AddTimer &
    !        ( '__Reconstruction_G', I % iTimerReconstruction_G )
    ! call PROGRAM_HEADER % AddTimer &
    !        ( '__Boundary', I % iTimerBoundary )
    ! call PROGRAM_HEADER % AddTimer &
    !        ( '__Reconstruction_CSL', I % iTimerReconstruction_CSL )
    ! call PROGRAM_HEADER % AddTimer &
    !        ( '__FromPrimitive', I % iTimerFromPrimitive )

    ! ! call PROGRAM_HEADER % AddTimer &
    ! !        ( 'ComputeReconstructionLinear_CSL', I % iTimerReconstruction_CSL )
    ! call PROGRAM_HEADER % AddTimer &
    !        ( 'Gradient', I % iTimerGradient )
    ! call PROGRAM_HEADER % AddTimer &
    !        ( 'ReconstructionKernel', I % iTimerReconstructionKernel )

    I % ReconstructionType = 'LINEAR'
    call PROGRAM_HEADER % GetParameter &
           ( I % ReconstructionType, 'ReconstructionType' )

    I % UseIncrementStream = .false.
    call PROGRAM_HEADER % GetParameter &
           ( I % UseIncrementStream, 'UseIncrementStream' )

    if ( I % UseIncrementStream ) then

      I % iStream = ATLAS % MAX_STREAMS

      OutputDirectory = '../Output/'
      call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )

      allocate ( I % GridImageStream )
      associate ( GIS => I % GridImageStream )
      call GIS % Initialize &
         ( 'Increment', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           WorkingDirectoryOption = OutputDirectory )
      end associate !-- GIS

    end if

  end subroutine Initialize


  subroutine SetStorage ( I, StorageDivergence )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    class ( StorageDivergenceForm ), intent ( in ), target :: &
      StorageDivergence

    I % Storage => StorageDivergence

  end subroutine SetStorage


  subroutine ClearStorage ( I )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I

    nullify ( I % Storage )

  end subroutine ClearStorage


  subroutine SetBoundaryFluence_CSL ( I, BoundaryFluence_CSL )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ), target :: &
        BoundaryFluence_CSL

    I % BoundaryFluence_CSL => BoundaryFluence_CSL
    
  end subroutine SetBoundaryFluence_CSL


  subroutine ClearBoundaryFluence ( I )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I

    nullify ( I % BoundaryFluence_CSL )

  end subroutine ClearBoundaryFluence


  subroutine SetMetricDerivativesFlat ( I, dLogVolumeJacobian_dX )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( Real_1D_Form ), dimension ( : ), intent ( in ), target :: &
      dLogVolumeJacobian_dX

    I % dLogVolumeJacobian_dX => dLogVolumeJacobian_dX

  end subroutine SetMetricDerivativesFlat


  subroutine SetUseLimiter ( I, UseLimiter )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      UseLimiter

    I % UseLimiter => UseLimiter

  end subroutine SetUseLimiter


  subroutine ClearMetricDerivatives ( I )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I

    nullify ( I % dLogVolumeJacobian_dX ) 

  end subroutine ClearMetricDerivatives


  subroutine Compute ( I, Increment, TimeStep, Weight_RK )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( StorageForm ), intent ( inout ) :: &
      Increment  !-- Assume Increment is already cleared!
    real ( KDR ), intent ( in ) :: &
      TimeStep
    real ( KDR ), intent ( in ) :: &
      Weight_RK

    integer ( KDI ) :: &
      iD  !-- iDimension
    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_DTH

    Timer => PROGRAM_HEADER % TimerPointer ( I % Storage % iTimerDivergence )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call Show ( 'Computing an IncrementDivergence', I % IGNORABILITY + 3 )
    call Show ( I % Name, 'Name', I % IGNORABILITY + 3 )

    if ( I % UseIncrementStream ) &
      call Show ( '>>> Entering Increment % Compute' )

    call SetCurrent ( I )

    call PrepareOutput ( I )

    !-- Assume Increment is already cleared!

    I % Weight_RK = Weight_RK

    do iD = 1, I % Chart % nDimensions

      if ( I % UseIncrementStream ) &
        call Show ( iD, '>>> iDimension' )

      call ComputeReconstruction ( I, iD )
      call ComputeFluxes ( I, iD )
      
      select type ( Chart => I % Chart )
      class is ( Chart_SL_Template )
        call ComputeIncrement_CSL ( I, Increment, Chart, TimeStep, iD )
      end select !-- Grid
      
      if ( I % UseIncrementStream ) then
        Timer_DTH => PROGRAM_HEADER % TimerPointer &
                       ( I % Storage % iTimerDataToHost )
        call Timer_DTH % Start ( )
        call I % Current % UpdateHost ( )
        call Increment % UpdateHost ( )
        call Timer_DTH % Stop ( )
        call Copy ( I % Current % Value, &
                    I % Output ( iCURRENT ) % Value )
        call Copy ( Increment % Value, &
                    I % Output ( iINCREMENT ) % Value )
        call WriteOutput ( I, I % Chart )
      end if

    end do !-- iD

    if ( I % UseIncrementStream ) &
      call Show ( '>>> Leaving Increment % Compute' )

    if ( associated ( Timer ) ) call Timer % Stop ( )
    
  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % GridImageStream ) ) &
      deallocate ( I % GridImageStream )
    if ( allocated ( I % Output ) ) &
      deallocate ( I % Output )

    if ( I % Name == '' ) return

    call Show ( 'Finalizing an IncrementDivergence_FV', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )
    
  end subroutine Finalize


  subroutine SetCurrent ( I )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I

    select type ( Chart => I % Chart )
    class is ( Chart_SL_Template )
      I % Geometry => Chart % Geometry ( )
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCurrent', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    select type ( CurrentChart => I % CurrentChart )
    class is ( Field_CSL_Template )

      select type ( Current => CurrentChart % Field )
      class is ( CurrentTemplate )
        I % Current => Current
      class default
        call Show ( 'Field is not a Current', CONSOLE % ERROR )
        call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetCurrent', 'subroutine', CONSOLE % ERROR ) 
        call PROGRAM_HEADER % Abort ( )
      end select !-- Current

    class default
      call Show ( 'CurrentChart type not found', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCurrent', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    if ( .not. associated ( I % Storage ) ) then
      call Show ( 'Storage is not set', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCurrent', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine SetCurrent


  subroutine PrepareOutput ( I )

    type ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iF, &  !-- iField
      iFG    !-- iFieldStorage
    character ( LDL ), dimension ( : ), allocatable :: &
      PrimitiveVariable, &
      ConservedVariable

    if ( .not. I % UseIncrementStream ) &
      return
    if ( allocated ( I % Output ) ) &
      return

    associate ( C => I % Current )

    allocate ( I % Output ( nOUTPUT ) )

    allocate ( ConservedVariable ( C % N_CONSERVED ) )
    ConservedVariable &
      = [ ( C % Variable ( C % iaConserved ( iF ) ), &
            iF = 1, C % N_CONSERVED ) ]

    !-- Field groups

    call I % Output ( iCURRENT ) % Initialize &
           ( [ C % nValues, C % nVariables ], &
             VariableOption = C % Variable, NameOption = 'Current' )
    call I % Output ( iCURRENT_IL ) % Initialize &
           ( [ C % nValues, C % nVariables ], &
             VariableOption = C % Variable, NameOption = 'Current_IL' )
    call I % Output ( iCURRENT_IR ) % Initialize &
           ( [ C % nValues, C % nVariables ], &
             VariableOption = C % Variable, NameOption = 'Current_IR' )
    call I % Output ( iFLUX_IL ) % Initialize &
           ( [ C % nValues, C % N_CONSERVED ], &
             VariableOption = ConservedVariable, NameOption = 'Flux_IL' )
    call I % Output ( iFLUX_IR ) % Initialize &
           ( [ C % nValues, C % N_CONSERVED ], &
             VariableOption = ConservedVariable, NameOption = 'Flux_IR' )
    call I % Output ( iFLUX_I ) % Initialize &
           ( [ C % nValues, C % N_CONSERVED ], &
             VariableOption = ConservedVariable, NameOption = 'Flux_I' )
    call I % Output ( iINCREMENT ) % Initialize &
           ( [ C % nValues, C % N_CONSERVED ], &
             VariableOption = ConservedVariable, NameOption = 'Update' )

    end associate !-- C

    !-- Open stream

    select type ( Chart => I % Chart )
    class is ( Chart_SL_Template )
      if ( .not. allocated ( Chart % Stream ( I % iStream ) % Element ) ) &
      then 
        call Chart % OpenStream &
               ( I % GridImageStream, 'Increment', I % iStream )
        do iFG = 1, nOUTPUT
          call Chart % AddFieldImage &
                 ( I % Output ( iFG ), I % iStream )
        end do !-- iFG
      end if !-- allocated stream
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'PrepareOutput', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )        
    end select !-- Grid

  end subroutine PrepareOutput


  subroutine WriteOutput ( I, Grid )

    type ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    class ( * ), intent ( inout ) :: &
      Grid

    if ( .not. I % UseIncrementStream ) &
      return

    associate ( GIS => I % GridImageStream )
    call GIS % Open ( GIS % ACCESS_CREATE )

    select type ( Grid )
    class is ( Chart_SL_Template )
      call Grid % Write ( I % iStream )
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'PrepareOutput', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )        
    end select !-- Grid

    call GIS % Close ( )
    end associate !-- GIS

  end subroutine WriteOutput


  subroutine ComputeReconstruction ( I, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iDD_22, iDD_33, &
      iUU_22, iUU_33
   type ( StorageForm ) :: &
     Reconstructed
    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_DTH

    Timer => PROGRAM_HEADER % TimerPointer &
               ( I % Storage % iTimerReconstruction )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call Show ( 'Computing Reconstruction', I % IGNORABILITY + 4 )

    associate &
      ( C    => I % Current, &
        A    => I % Chart % Atlas, &
        G    => I % Geometry, &
        C_IL => I % Storage % Current_IL, &
        C_IR => I % Storage % Current_IR, &
        G_I  => I % Storage % Geometry_I )

    select type ( CSL => I % Chart )
    class is ( ChartHeader_SL_Form )
    call G % ComputeReconstruction &
           ( G_I, CSL, I % Chart % nDimensions, iDimension )
    end select !-- CSL

    call Reconstructed % Initialize &
           ( C, iaSelectedOption = C % iaReconstructed )

    associate &
      ( iaI => A % Connectivity % iaInner ( iDimension ), &
        iaO => A % Connectivity % iaOuter ( iDimension ) )

    select type ( A )
    class is ( Atlas_SC_Template )
      call A % ApplyBoundaryConditions ( Reconstructed, iDimension, iaI )
      call A % ApplyBoundaryConditions ( Reconstructed, iDimension, iaO )
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeReconstructionLinear', 'subroutine', CONSOLE % ERROR )
    end select !-- A

    select type ( Chart => I % Chart )
    class is ( Chart_SL_Template )
      select case ( trim ( I % ReconstructionType ) )
      case ( 'LINEAR' )
        call ComputeReconstructionLinear_CSL &
               ( I, Reconstructed, Chart, iDimension )
      case ( 'PARABOLIC' )
!call Show ( '>>> 1' )
        call ComputeReconstructionParabolic_CSL &
               ( I, Reconstructed, Chart, iDimension )
!call Show ( '>>> 2' )
!call Show ( '>>> Stopping after Reconstruction', CONSOLE % WARNING )
!call PROGRAM_HEADER % Abort ( )
      case default
        call Show ( 'ReconstructionType not recognized', CONSOLE % ERROR )
        call Show ( I % ReconstructionType, 'ReconstructionType', &
                    CONSOLE % ERROR )
        call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeReconstruction', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- ReconstructionType
    end select !-- Grid

    if ( trim ( C % ReconstructedType ) == 'PRIMITIVE' ) then
!call Show ( '>>> 3' )
      call C % ComputeFromPrimitive ( C_IL, G, G_I )
!call Show ( '>>> 4' )
      call C % ComputeFromPrimitive ( C_IR, G, G_I )
!call Show ( '>>> 5' )
    else if ( trim ( C % ReconstructedType ) == 'CONSERVED' ) then
      call C % ComputeFromConserved ( C_IL, G, G_I )
      call C % ComputeFromConserved ( C_IR, G, G_I )
    end if

    end associate !-- iaI, iaO

    if ( I % UseIncrementStream ) then
      Timer_DTH => PROGRAM_HEADER % TimerPointer &
                     ( I % Storage % iTimerDataToHost )
      call Timer_DTH % Start ( )
      call C_IL % UpdateHost ( )
      call C_IR % UpdateHost ( )
      call Timer_DTH % Stop ( )
      call Copy ( C_IL % Value, I % Output ( iCURRENT_IL ) % Value )
      call Copy ( C_IR % Value, I % Output ( iCURRENT_IR ) % Value )
    end if

    end associate !-- C, etc.

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeReconstruction


  subroutine ComputeFluxes ( I, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      iDimension

    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_DTH

    Timer => PROGRAM_HEADER % TimerPointer ( I % Storage % iTimerFluxes )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call Show ( 'Computing Fluxes', I % IGNORABILITY + 4 )

    associate &
      ( C    => I % Current, &
        G    => I % Geometry, &
        C_IL => I % Storage % Current_IL, &
        C_IR => I % Storage % Current_IR, &
        F_I  => I % Storage % Flux_I, &
        F_IL => I % Storage % Flux_IL, &
        F_IR => I % Storage % Flux_IR, &
        DF_I => I % Storage % DiffusionFactor_I, &
        SS_I => I % Storage % SolverSpeeds_I, &
        G_I  => I % Storage % Geometry_I )

    call C % ComputeFluxes &
           ( I, F_I, F_IL, F_IR, SS_I, DF_I, I % Chart, G, C_IL, C_IR, G_I, &
             iDimension )

    if ( I % UseIncrementStream ) then
      Timer_DTH => PROGRAM_HEADER % TimerPointer &
                     ( I % Storage % iTimerDataToHost )
      call Timer_DTH % Start ( )
      call F_IL % UpdateHost ( )
      call F_IR % UpdateHost ( )
      call F_I  % UpdateHost ( )
      call Timer_DTH % Stop ( )
      call Copy ( F_IL % Value, I % Output ( iFLUX_IL ) % Value )
      call Copy ( F_IR % Value, I % Output ( iFLUX_IR ) % Value )
      call Copy ( F_I % Value, I % Output ( iFLUX_I ) % Value )
    end if

    end associate !-- C, etc.

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeFluxes


  subroutine ComputeReconstructionLinear_CSL &
               ( I, Reconstructed, CSL, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    class ( StorageForm ), intent ( in ) :: &
      Reconstructed
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iF  !-- iField
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dX_L, dX_R, &
      V, &
      dVdX, &
      V_IL, &
      V_IR, &
      A_I, &
      dLVdX

!    associate &
!      ( Timer => PROGRAM_HEADER % Timer ( I % iTimerReconstruction_CSL ) )
!    call Timer % Start ( )

    associate &
      ( C    => I % Current, &
        G    => I % Geometry, &
        Grad => I % Storage % GradientReconstructed, &
        C_IL => I % Storage % Current_IL, &
        C_IR => I % Storage % Current_IR, &
        G_I  => I % Storage % Geometry_I )

!    associate ( Timer_G => PROGRAM_HEADER % Timer ( I % iTimerGradient ) )
!    call Timer_G % Start ( )
    if ( C % UseLimiter ) then
      if ( associated ( I % UseLimiter ) ) then
        call Grad % Compute &
               ( CSL, Reconstructed, iDimension, &
                 UseLimiterOption = I % UseLimiter, &
                 LimiterParameterOption = C % LimiterParameter )
      else
        call Grad % Compute &
               ( CSL, Reconstructed, iDimension, &
                 LimiterParameterOption = C % LimiterParameter )
      end if
    else
      call Grad % Compute ( CSL, Reconstructed, iDimension )
    end if
!    call Timer_G % Stop
!    end associate !-- Timer_G  

    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_LEFT_U ( iDimension ) ), dX_L )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_RIGHT_U ( iDimension ) ), dX_R )

    !-- Reconstruct Current

!    associate ( Timer_RK => PROGRAM_HEADER % Timer &
!                              ( I % iTimerReconstructionKernel ) )
!    call Timer_RK % Start ( )

    associate ( iaR => C % iaReconstructed )
    do iF = 1, C % N_RECONSTRUCTED
      call CSL % SetVariablePointer &
             ( C % Value ( :, iaR ( iF ) ), V )
      call CSL % SetVariablePointer &
             ( Grad % Output % Value ( :, iF ), dVdX )
      call CSL % SetVariablePointer &
             ( C_IL % Value ( :, iaR ( iF ) ), V_IL )
      call CSL % SetVariablePointer &
             ( C_IR % Value ( :, iaR ( iF ) ), V_IR )
      call ComputeReconstructionLinear_CSL_Kernel &
             ( V, dVdX, dX_L, dX_R, iDimension, &
               CSL % nGhostLayers ( iDimension ), V_IL, V_IR, &
               UseDeviceOption = C % AllocatedDevice )
    end do !-- iF
    end associate !-- iaR

!    call Timer_RK % Stop
!    end associate !-- Timer_RK

    !-- VolumeJacobian derivative

    if ( associated ( I % dLogVolumeJacobian_dX ) ) then
      if ( size ( I % dLogVolumeJacobian_dX ) >= iDimension ) then
        call CSL % SetVariablePointer &
               ( G % Value ( :, G % VOLUME ), V )
        call CSL % SetVariablePointer &
               ( G % Value ( :, G % AREA_INNER_D ( iDimension ) ), A_I )
        call CSL % SetVariablePointer &
               ( I % dLogVolumeJacobian_dX ( iDimension ) % Value, dLVdX )
        call ComputeLogDerivative_CSL_Kernel &
               ( A_I, V, iDimension, CSL % nGhostLayers ( iDimension ), &
                 dLVdX, UseDeviceOption = G % AllocatedDevice )
      end if
    end if

    nullify ( dX_L, dX_R, V, dVdX, V_IL, V_IR, A_I, dLVdX )

    end associate !-- C, etc.

!    call Timer % Stop
!    end associate !-- Timer
  
  end subroutine ComputeReconstructionLinear_CSL


  subroutine ComputeReconstructionParabolic_CSL &
               ( I, Reconstructed, CSL, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    class ( StorageForm ), intent ( in ) :: &
      Reconstructed
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iF  !-- iField
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      X, &
      dX_L, dX_R, &
      XA, X2A, &
      V, &
      V_IL, &
      V_IR, &
      A_I, &
      dLVdX

    associate &
      ( C    => I % Current, &
        G    => I % Geometry, &
        C_IL => I % Storage % Current_IL, &
        C_IR => I % Storage % Current_IR, &
        G_I  => I % Storage % Geometry_I )

    call CSL % SetVariablePointer &
           ( G % Value ( :, G % CENTER_U ( iDimension ) ), X )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_LEFT_U ( iDimension ) ), dX_L )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_RIGHT_U ( iDimension ) ), dX_R )

    !-- Reconstruct Current

    associate ( iaR => C % iaReconstructed )
!call Show ( '>>> 1.1' )
!call Show ( C % N_RECONSTRUCTED, '>>> N_RECONSTRUCTED' )
    do iF = 1, C % N_RECONSTRUCTED
      call CSL % SetVariablePointer &
             ( C % Value ( :, iaR ( iF ) ), V )
      call CSL % SetVariablePointer &
             ( C_IL % Value ( :, iaR ( iF ) ), V_IL )
      call CSL % SetVariablePointer &
             ( C_IR % Value ( :, iaR ( iF ) ), V_IR )
!call Show ( iF, '>>> iF' )
      call ComputeReconstructionParabolic_CSL_Kernel &
             ( V, X, dX_L, dX_R, X, X ** 2, iDimension, &
               CSL % nGhostLayers ( iDimension ), V_IL, V_IR, &
               UseDeviceOption = C % AllocatedDevice )
!call Show ( '>>> 1.2' )
    end do !-- iF
    end associate !-- iaR

    !-- VolumeJacobian derivative

    if ( associated ( I % dLogVolumeJacobian_dX ) ) then
      if ( size ( I % dLogVolumeJacobian_dX ) >= iDimension ) then
        call CSL % SetVariablePointer &
               ( G % Value ( :, G % VOLUME ), V )
        call CSL % SetVariablePointer &
               ( G % Value ( :, G % AREA_INNER_D ( iDimension ) ), A_I )
        call CSL % SetVariablePointer &
               ( I % dLogVolumeJacobian_dX ( iDimension ) % Value, dLVdX )
        call ComputeLogDerivative_CSL_Kernel &
               ( A_I, V, iDimension, CSL % nGhostLayers ( iDimension ), &
                 dLVdX, UseDeviceOption = G % AllocatedDevice )
      end if
    end if

    nullify ( X, dX_L, dX_R, XA, X2A, V, V_IL, V_IR, A_I, dLVdX )

    end associate !-- C, etc.

  end subroutine ComputeReconstructionParabolic_CSL


  subroutine ComputeIncrement_CSL &
               ( I, Increment, CSL, TimeStep, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( StorageForm ), intent ( inout ) :: &
      Increment
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iF  !-- iField
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dU, &
      F_I, &
      A_I, &
      V
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( I % Storage % iTimerIncrement )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call Show ( 'Computing Increment', I % IGNORABILITY + 4 )

    associate &
      ( C => I % Current, &
        G => I % Geometry )

    call CSL % SetVariablePointer &
           ( G % Value ( :, G % VOLUME ), V )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % AREA_INNER_D ( iDimension ) ), A_I )

    if ( .not. associated ( I % BoundaryFluence_CSL ) ) then
      call Show ( 'Increment % Set has not been called', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeIncrement_CSL', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    do iF = 1, C % N_CONSERVED
      call CSL % SetVariablePointer &
             ( Increment % Value ( :, iF ), dU )
      call CSL % SetVariablePointer &
             ( I % Storage % Flux_I % Value ( :, iF ), F_I )
      call ComputeIncrement_CSL_Kernel &
             ( dU, F_I, A_I, V, TimeStep, iDimension, &
               CSL % nGhostLayers ( iDimension ), &
               UseDeviceOption = Increment % AllocatedDevice )
      if ( associated ( I % BoundaryFluence_CSL ) ) &
        call RecordBoundaryFluence_CSL &
               ( I % BoundaryFluence_CSL, CSL, F_I, I % Weight_RK, &
                 TimeStep, iDimension, iF )
    end do !-- iF

    end associate !-- C, etc.
    nullify ( dU, F_I, A_I, V )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeIncrement_CSL


  subroutine RecordBoundaryFluence_CSL &
               ( BF, CSL, F_I, Weight_RK, dT, iD, iC )

    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      BF
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_I
    real ( KDR ), intent ( in ) :: &
      Weight_RK, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iD, &  !-- iDimension
      iC     !-- iConserved

    integer ( KDI ) :: &
      jD, kD, &   !-- jDimension, kDimension
      nCells
    integer ( KDI ), dimension ( 3 ) :: &
      oB, & !-- oBoundary
      nB    !-- nBoundary
    logical ( KDL ) :: &
      RecordInner, &
      RecordOuter

    select type ( CSL )
    class is ( Chart_SLL_Form )
      RecordInner = .true.
      RecordOuter = .true.
      nCells = CSL % nCells ( iD )
    class is ( Chart_SLD_Form )
      RecordInner = ( CSL % iaBrick ( iD ) == 1 )
      RecordOuter = ( CSL % iaBrick ( iD ) == CSL % nBricks ( iD ) )
      nCells = CSL % nCellsBrick ( iD )
    end select !-- CSL

    jD = mod ( iD, 3 ) + 1
    kD = mod ( jD, 3 ) + 1
    
    nB ( iD ) = 1
    nB ( jD ) = CSL % nCellsBrick ( jD )
    nB ( kD ) = CSL % nCellsBrick ( kD )

    if ( RecordInner ) then
      associate ( iCI => CSL % Atlas % Connectivity % iaInner ( iD ) )
      associate ( BF_Inner => BF ( iC, iCI ) % Value )
      oB = CSL % nGhostLayers
      call RecordBoundaryFluence_CSL_Kernel &
             ( BF_Inner, F_I, Weight_RK * dT, nB, oB, &
               UseDeviceOption = BF ( iC, iCI ) % AllocatedDevice )
      end associate !-- BF_Inner
      end associate !-- iCI
    end if !-- iaBrick ( iD ) == 1

    if ( RecordOuter ) then
      associate ( iCO => CSL % Atlas % Connectivity % iaOuter ( iD ) )
      associate ( BF_Outer => BF ( iC, iCO ) % Value )
      oB        = CSL % nGhostLayers
      oB ( iD ) =  oB ( iD ) + nCells
      call RecordBoundaryFluence_CSL_Kernel &
             ( BF_Outer, F_I, Weight_RK * dT, nB, oB, &
               UseDeviceOption = BF ( iC, iCO ) % AllocatedDevice )
      end associate !-- BF_Outer
      end associate !-- iCO
    end if !-- iaBrick ( iD ) == nBricks ( iD )

  end subroutine RecordBoundaryFluence_CSL


end module IncrementDivergence_FV__Form
