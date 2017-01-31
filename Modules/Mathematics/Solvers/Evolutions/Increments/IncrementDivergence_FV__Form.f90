!-- IncrementDivergence_FV computes a first-order (in time) update due the
!   divergence in a conservation law.

module IncrementDivergence_FV__Form

  !-- IncrementDivergence_FiniteVolume_Form

  use Basics
  use Manifolds
  use Operations
  use Fields

  implicit none
  private

  type, public :: IncrementDivergence_FV_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerReconstruction, &
!      iTimerReconstruction_CSL, &
      iTimerFromPrimitive, &
!      iTimerBoundary, &
      iTimerGradient, &
      iTimerReconstructionKernel, &
      iTimerFluxes, &
      iTimerIncrement, &
      iStream
    integer ( KDI ) :: &
      ALPHA_PLUS  = 1, &
      ALPHA_MINUS = 2, &
      N_MODIFIED_SPEEDS = 2
    real ( KDR ) :: &
      LimiterParameter, &
      Weight_RK = - huge ( 0.0_KDR ) !-- RungeKutta weight
    logical ( KDL ) :: &
      UseLimiterParameter = .true., &
      UseIncrementStream = .false.
    character ( LDF ) :: &
      Name = ''
    type ( Real_1D_Form ), dimension ( : ), pointer :: &
      dLogVolumeJacobian_dX => null ( )
    type ( Real_3D_Form ), dimension ( :, : ), pointer :: &
      BoundaryFluence_CSL => null ( )
    type ( VariableGroupForm ), allocatable :: &
      ModifiedSpeeds_I, &  !-- ModifiedSpeed_Inner
      DiffusionFactor_I
    type ( VariableGroupForm ), dimension ( : ), allocatable :: &
      Output
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      SetLimiterParameter
    procedure, private, pass :: &
      SetBoundaryFluence_CSL
    procedure, private, pass :: &
      Set_dLog_VJ
    procedure, private, pass :: &
      Set_Weight_RK
    generic, public :: &
      Set => SetLimiterParameter, SetBoundaryFluence_CSL, Set_dLog_VJ, &
             Set_Weight_RK
    procedure, private, pass :: &
      Clear_CSL
    generic, public :: &
      Clear => Clear_CSL
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type IncrementDivergence_FV_Form

    private :: &
      PrepareOutput, &
      WriteOutput, &
      ComputeReconstruction, &
      ComputeFluxes

      private :: &
        ComputeReconstruction_CSL, &
        ComputeIncrement_CSL

        private :: &
          ComputeReconstruction_CSL_Kernel, &
          ComputeLogDerivative_CSL_Kernel, &
          ComputeFluxesKernel, &
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

contains


  subroutine Initialize ( I, NameSuffix )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    character ( * ), intent ( in ) :: &
      NameSuffix

    character ( LDF ) :: &
      OutputDirectory

    I % IGNORABILITY = CONSOLE % INFO_4
    I % Name = 'IncrementDivergence_FV_' // trim ( NameSuffix )

    call Show ( 'Initializing an IncrementDivergence_FV', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    I % LimiterParameter = 1.4_KDR
    call PROGRAM_HEADER % GetParameter &
           ( I % LimiterParameter, 'LimiterParameter' )

    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeReconstruction', I % iTimerReconstruction )
    ! call PROGRAM_HEADER % AddTimer &
    !        ( 'ComputeReconstruction_CSL', I % iTimerReconstruction_CSL )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeFromPrimitive', I % iTimerFromPrimitive )
    ! call PROGRAM_HEADER % AddTimer &
    !        ( 'ApplyBoundary', I % iTimerBoundary )
    call PROGRAM_HEADER % AddTimer &
           ( 'Gradient', I % iTimerGradient )
    call PROGRAM_HEADER % AddTimer &
           ( 'ReconstructionKernel', I % iTimerReconstructionKernel )

    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeFluxes', I % iTimerFluxes )
    call PROGRAM_HEADER % AddTimer &
           ( 'ApplyFluxes', I % iTimerIncrement )

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


  subroutine SetLimiterParameter ( I, UseLimiterParameter )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ) :: &
      UseLimiterParameter

    I % UseLimiterParameter = UseLimiterParameter

  end subroutine SetLimiterParameter


  subroutine SetBoundaryFluence_CSL ( I, BoundaryFluence_CSL )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ), target :: &
      BoundaryFluence_CSL

    I % BoundaryFluence_CSL => BoundaryFluence_CSL

  end subroutine SetBoundaryFluence_CSL


  subroutine Set_dLog_VJ ( I, dLogVolumeJacobian_dX )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( Real_1D_Form ), dimension ( : ), intent ( in ), target :: &
      dLogVolumeJacobian_dX

    I % dLogVolumeJacobian_dX => dLogVolumeJacobian_dX

  end subroutine Set_dLog_VJ


  subroutine Set_Weight_RK ( I, Weight_RK )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      Weight_RK

    I % Weight_RK = Weight_RK

  end subroutine Set_Weight_RK


  subroutine Clear_CSL ( I )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I

    I % BoundaryFluence_CSL   => null ( )
    I % dLogVolumeJacobian_dX => null ( )

    I % Weight_RK =  - huge ( 1.0_KDR )

    I % UseLimiterParameter = .true.

  end subroutine Clear_CSL


  subroutine Compute ( I, Increment, Grid, C, TimeStep )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      Increment  !-- Assume Increment is already cleared!
    class ( * ), intent ( inout ), target :: &
      Grid
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      nDimensions
    type ( VariableGroupForm ), allocatable :: &
      G_I, &         !-- Geometry_Inner
      C_IL, C_IR, &  !-- Current_InnerLeft, Current_InnerRight
      F_I            !-- Flux_Inner
    class ( GeometryFlatForm ), pointer :: &
      G

    if ( I % UseIncrementStream ) &
      call Show ( '>>> Entering Increment % Compute' )

    call PrepareOutput ( I, Grid, C )

    select type ( Grid )
    class is ( Chart_SLD_Form )
      nDimensions = Grid % nDimensions
      G => Grid % Geometry ( )
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    allocate ( G_I )
    call G_I % Initialize ( [ G % nValues, G % nVariables ] )

    allocate ( C_IL, C_IR )
    call C_IL % Initialize &
           ( [ C % nValues, C % nVariables ], &
             VectorIndicesOption = C % VectorIndices, ClearOption = .true. )
    call C_IR % Initialize &
           ( [ C % nValues, C % nVariables ], &
             VectorIndicesOption = C % VectorIndices, ClearOption = .true. )

    allocate ( F_I )
    call F_I % Initialize ( [ C % nValues, C % N_CONSERVED ] )

    !-- Assume Increment is already cleared!

    do iD = 1, nDimensions

      if ( I % UseIncrementStream ) &
        call Show ( iD, '>>> iDimension' )

      select type ( Grid )
      class is ( Chart_SLD_Form )
        call ComputeReconstruction &
               ( I, C_IL, C_IR, G_I, C, Grid, iD )
        call ComputeFluxes &
               ( I, F_I, C, Grid, G, C_IL, C_IR, G_I, iD )
        call ComputeIncrement_CSL &
               ( I, Increment, C, F_I, G_I, Grid, TimeStep, iD )
      end select !-- Grid

      if ( I % UseIncrementStream ) then
        call Copy ( C % Value, &
                    I % Output ( iCURRENT ) % Value )
        call Copy ( Increment % Value, &
                    I % Output ( iINCREMENT ) % Value )
        call WriteOutput ( I, Grid )
      end if

    end do !-- iD

    nullify ( G )

    if ( I % UseIncrementStream ) &
      call Show ( '>>> Leaving Increment % Compute' )

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


  subroutine PrepareOutput ( I, Grid, C )

    type ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    class ( * ), intent ( inout ) :: &
      Grid
    class ( CurrentTemplate ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      iFG    !-- iFieldGroup
    character ( LDL ), dimension ( : ), allocatable :: &
      PrimitiveVariable, &
      ConservedVariable

    if ( .not. I % UseIncrementStream ) &
      return
    if ( allocated ( I % Output ) ) &
      return

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

    !-- Open stream

    select type ( Grid )
    class is ( Chart_SL_Template )
      if ( .not. allocated ( Grid % Stream ( I % iStream ) % Element ) ) &
      then 
        call Grid % OpenStream &
               ( I % GridImageStream, 'Increment', I % iStream )
        do iFG = 1, nOUTPUT
          call Grid % AddFieldImage &
                 ( I % Output ( iFG ), I % iStream )
        end do !-- iFG
      end if !-- allocated stream
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
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


  subroutine ComputeReconstruction &
               ( I, C_IL, C_IR, G_I, C, Grid, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      C_IL, C_IR, &  !-- Current
      G_I           !-- VolumeJacobian
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( * ), intent ( in ), target :: &
      Grid
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iDD_22, iDD_33, &
      iUU_22, iUU_33, &
      nDimensions
    type ( VariableGroupForm ), allocatable :: &
      P, &
      P_IL, &  !-- Primitive_InnerLeft
      P_IR     !-- Primitive_InnerRight
    class ( AtlasHeaderForm ), pointer :: &
      A
    class ( GeometryFlatForm ), pointer :: &
      G

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( I % iTimerReconstruction ) )
    call Timer % Start ( )

    select type ( Grid )
    class is ( Chart_SL_Template )
      nDimensions = Grid % nDimensions
      A => Grid % Atlas
      G => Grid % Geometry ( )
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeReconstruction', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    call G % ComputeReconstruction ( G_I, nDimensions, iDimension )

    allocate ( P )
    allocate ( P_IL )
    allocate ( P_IR )
    call P % Initialize ( C, iaSelectedOption = C % iaPrimitive )
    call P_IL % Initialize ( C_IL, iaSelectedOption = C % iaPrimitive )
    call P_IR % Initialize ( C_IR, iaSelectedOption = C % iaPrimitive )

    associate &
      ( iaI => A % Connectivity % iaInner ( iDimension ), &
        iaO => A % Connectivity % iaOuter ( iDimension ) )

!    associate &
!      ( Timer_B => PROGRAM_HEADER % Timer ( I % iTimerBoundary ) )
    
!    call Timer_B % Start ( )
    select type ( A )
    class is ( Atlas_SC_Template )
      call A % ApplyBoundaryConditions ( P, iDimension, iaI )
      call A % ApplyBoundaryConditions ( P, iDimension, iaO )
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Update_CF_E__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeReconstruction', 'subroutine', CONSOLE % ERROR )
    end select !-- A
!    call Timer_B % Stop

    select type ( Grid )
    class is ( Chart_SL_Template )
      call ComputeReconstruction_CSL &
             ( I, C_IL, C_IR, C, P, G_I, Grid, iDimension )
    end select !-- Grid

    associate &
      ( Timer_FP => PROGRAM_HEADER % Timer ( I % iTimerFromPrimitive ) )
    call Timer_FP % Start ( )
    call C % ComputeFromPrimitive ( C_IL % Value, G, G_I % Value )
    call C % ComputeFromPrimitive ( C_IR % Value, G, G_I % Value )
    call Timer_FP % Stop ( )
    end associate !-- Timer_FP

    end associate !-- iaI, iaO

    nullify ( A, G )

    if ( I % UseIncrementStream ) then
      call Copy ( C_IL % Value, I % Output ( iCURRENT_IL ) % Value )
      call Copy ( C_IR % Value, I % Output ( iCURRENT_IR ) % Value )
    end if

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeReconstruction


  subroutine ComputeFluxes &
               ( I, Flux_I, C, Grid, G, C_IL, C_IR, G_I, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      Flux_I
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( * ), intent ( in ) :: &
      Grid
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( VariableGroupForm ), intent ( in ) :: &
      C_IL, C_IR, &  !-- ConservedFields
      G_I            !-- Geometry
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iF  !-- iField
    type ( VariableGroupForm ), allocatable :: &
      Flux_IL, &
      Flux_IR

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( I % iTimerFluxes ) )
    call Timer % Start ( )

    allocate ( I % ModifiedSpeeds_I )
    call I % ModifiedSpeeds_I % Initialize &
           ( [ C % nValues, I % N_MODIFIED_SPEEDS ] )!, ClearOption = .true. )

    allocate ( I % DiffusionFactor_I )
    call I % DiffusionFactor_I % Initialize &
           ( [ C % nValues, C % N_CONSERVED ] )!, ClearOption = .true. )

    allocate ( Flux_IL, Flux_IR )
    call Flux_IL % Initialize &
           ( [ C % nValues, C % N_CONSERVED ] )!, ClearOption = .true. )
    call Flux_IR % Initialize &
           ( [ C % nValues, C % N_CONSERVED ] )!, ClearOption = .true. )

    call C % ComputeRiemannSolverInput &
           ( I, I % DiffusionFactor_I, &
             I % ModifiedSpeeds_I % Value ( :, I % ALPHA_PLUS ), &
             I % ModifiedSpeeds_I % Value ( :, I % ALPHA_MINUS ), &
             Grid, C_IL, C_IR, iDimension )

    call C % ComputeRawFluxes &
           ( Flux_IL % Value, G, C_IL % Value, G_I % Value, iDimension )
    call C % ComputeRawFluxes &
           ( Flux_IR % Value, G, C_IR % Value, G_I % Value, iDimension )

!    call Clear ( Flux_I % Value )

    associate ( iaC => C % iaConserved )    
    do iF = 1, C % N_CONSERVED
      call ComputeFluxesKernel &
             ( Flux_IL % Value ( :, iF ), &
               Flux_IR % Value ( :, iF ), &
               C_IL % Value ( :, iaC ( iF ) ), &
               C_IR % Value ( :, iaC ( iF ) ), &
               I % ModifiedSpeeds_I % Value ( :, I % ALPHA_PLUS ), &
               I % ModifiedSpeeds_I % Value ( :, I % ALPHA_MINUS ), &
               I % DiffusionFactor_I % Value ( :, iF ), &
               Flux_I % Value ( :, iF ) )
    end do !-- iF
    end associate !-- iaC

    if ( I % UseIncrementStream ) then
      call Copy ( Flux_IL % Value, I % Output ( iFLUX_IL ) % Value )
      call Copy ( Flux_IR % Value, I % Output ( iFLUX_IR ) % Value )
      call Copy ( Flux_I % Value, I % Output ( iFLUX_I ) % Value )
    end if

    deallocate ( I % DiffusionFactor_I )
    deallocate ( I % ModifiedSpeeds_I )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeFluxes


  subroutine ComputeReconstruction_CSL &
               ( I, C_IL, C_IR, C, P, G_I, CSL, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( in ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      C_IL, C_IR  !-- Current
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    type ( VariableGroupForm ), intent ( in ) :: &
      P, &
      G_I
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iF  !-- iField
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dX, &
      V, &
      dVdX, &
      V_IL, &
      V_IR, &
      V_I, &
      dLVdX
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( GradientForm ), allocatable :: &
      Gradient_C

!    associate &
!      ( Timer => PROGRAM_HEADER % Timer ( I % iTimerReconstruction_CSL ) )
!    call Timer % Start ( )

    associate ( Timer_G => PROGRAM_HEADER % Timer ( I % iTimerGradient ) )
    call Timer_G % Start ( )
    allocate ( Gradient_C )
    call Gradient_C % Initialize ( P, 'Primitive' )
    if ( I % UseLimiterParameter ) then
      call Gradient_C % Compute &
             ( CSL, iDimension, LimiterParameterOption = I % LimiterParameter )
    else
      call Gradient_C % Compute ( CSL, iDimension )
    end if
    call Timer_G % Stop
    end associate !-- Timer_G  

    G => CSL % Geometry ( )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH ( iDimension ) ), dX )

    !-- Reconstruct Current

    associate ( Timer_RK => PROGRAM_HEADER % Timer &
                              ( I % iTimerReconstructionKernel ) )
    call Timer_RK % Start ( )

    associate ( iaP => C % iaPrimitive )
    do iF = 1, C % N_PRIMITIVE
      call CSL % SetVariablePointer &
             ( C % Value ( :, iaP ( iF ) ), V )
      call CSL % SetVariablePointer &
             ( Gradient_C % Output % Value ( :, iF ), dVdX )
      call CSL % SetVariablePointer &
             ( C_IL % Value ( :, iaP ( iF ) ), V_IL )
      call CSL % SetVariablePointer &
             ( C_IR % Value ( :, iaP ( iF ) ), V_IR )
      call ComputeReconstruction_CSL_Kernel &
             ( V, dVdX, dX, V_IL, V_IR, iDimension, &
               CSL % nGhostLayers ( iDimension ) )
    end do !-- iF
    end associate !-- iaP

    call Timer_RK % Stop
    end associate !-- Timer_RK

    deallocate ( Gradient_C )

    !-- VolumeJacobian derivative

    if ( associated ( I % dLogVolumeJacobian_dX ) ) then
      if ( size ( I % dLogVolumeJacobian_dX ) >= iDimension ) then
        !-- dX already set
        call CSL % SetVariablePointer &
               ( G % Value ( :, G % VOLUME_JACOBIAN ), V )
        call CSL % SetVariablePointer &
               ( G_I % Value ( :, G % VOLUME_JACOBIAN ), V_I )
        call CSL % SetVariablePointer &
               ( I % dLogVolumeJacobian_dX ( iDimension ) % Value, dLVdX )
        call ComputeLogDerivative_CSL_Kernel &
               ( V_I, V, dX, iDimension, CSL % nGhostLayers ( iDimension ), &
                 dLVdX )
      end if
    end if

    nullify ( G, V, dVdX, dX, V_IL, V_IR, V_I, dLVdX )

!    call Timer % Stop
!    end associate !-- Timer
  
  end subroutine ComputeReconstruction_CSL


  subroutine ComputeIncrement_CSL &
               ( I, Increment, C, Flux_I, G_I, CSL, TimeStep, iDimension )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    type ( VariableGroupForm ), intent ( in ) :: &
      Flux_I, &
      G_I
    class ( Chart_SLD_Form ), intent ( in ) :: &
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
      VJ_I, &
      VJ, &
      dX
    class ( GeometryFlatForm ), pointer :: &
      G

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( I % iTimerIncrement ) )
    call Timer % Start ( )

    G => CSL % Geometry ( )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH ( iDimension ) ), dX )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % VOLUME_JACOBIAN ), VJ )
    call CSL % SetVariablePointer &
           ( G_I % Value ( :, G % VOLUME_JACOBIAN ), VJ_I )

    if ( .not. associated ( I % BoundaryFluence_CSL ) ) then
      call Show ( 'Increment % Set has not been called', CONSOLE % ERROR )
      call Show ( 'IncrementDivergence_FV__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeIncrement_CSL', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    do iF = 1, C % N_CONSERVED
      call CSL % SetVariablePointer ( Increment % Value ( :, iF ), dU )
      call CSL % SetVariablePointer ( Flux_I % Value ( :, iF ), F_I )
      call ComputeIncrement_CSL_Kernel &
             ( dU, F_I, VJ_I, VJ, dX, TimeStep, iDimension, &
               CSL % nGhostLayers ( iDimension ) )
      call RecordBoundaryFluence_CSL &
             ( I % BoundaryFluence_CSL, CSL, F_I, VJ_I, I % Weight_RK, &
               TimeStep, iDimension, iF )
    end do !-- iF

    nullify ( G )
    nullify ( dU, F_I, VJ_I, VJ, dX )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeIncrement_CSL


  subroutine ComputeReconstruction_CSL_Kernel &
               ( V, dVdX, dX, V_IL, V_IR, iD, oV )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      V, &
      dVdX, &
      dX
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      V_IL, &
      V_IR
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV   
    
    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV

!    V_IL  =  cshift ( V  +  0.5_KDR * dX * dVdX, shift = -1, dim = iD )

!    V_IR  =  V  -  0.5_KDR * dX * dVdX

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = -1
      
    !$OMP parallel do private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          V_IL ( iV, jV, kV )  &
            =  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
               +  0.5_KDR  *  dX ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                           *  dVdX ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
      
    !$OMP parallel do private ( iV, jV, kV )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          V_IR ( iV, jV, kV )  &
            =  V ( iV, jV, kV )  &
               -  0.5_KDR  *  dX ( iV, jV, kV )  *  dVdX ( iV, jV, kV )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
    
  end subroutine ComputeReconstruction_CSL_Kernel


  subroutine ComputeLogDerivative_CSL_Kernel ( V_I, V, dX, iD, oV, dLVdX )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      V_I, &
      V, &
      dX
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV   
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dLVdX

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
      
    iaS = 0
    iaS ( iD ) = +1
      
    !$OMP parallel do private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          dLVdX ( iV, jV, kV ) &
            =  (    V_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                 -  V_I ( iV, jV, kV )  ) &
               / ( V ( iV, jV, kV ) * dX ( iV, jV, kV ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
    
  end subroutine ComputeLogDerivative_CSL_Kernel


  subroutine ComputeFluxesKernel &
               ( F_IL, F_IR, U_IL, U_IR, AP_I, AM_I, DF_I, F_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      F_IL, F_IR, &
      U_IL, U_IR, &
      AP_I, &
      AM_I, &
      DF_I
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      F_I

    integer ( KDI ) :: &
      iV, &
      nV  

    nV = size ( F_I )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( AP_I ( iV ) + AM_I ( iV ) > 0.0_KDR ) then
        F_I ( iV ) &
          =  (    AP_I ( iV ) * F_IL ( iV ) &
               +  AM_I ( iV ) * F_IR ( iV ) &
               -  DF_I ( iV ) * AP_I ( iV ) * AM_I ( iV ) &
                  * ( U_IR ( iV ) - U_IL ( iV ) ) ) &
             /  ( AP_I ( iV ) + AM_I ( iV ) )
      else
        F_I ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP end parallel do

  end subroutine ComputeFluxesKernel


  subroutine ComputeIncrement_CSL_Kernel &
               ( dU, F_I, VJ_I, VJ, dX, dT, iD, oV )
    
    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      dU
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_I, &
      VJ_I, &
      VJ, &
      dX
    real ( KDR ), intent ( in ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV   

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV

!    dU  =  dU  +  dT * ( VJ_I * F_I  &
!                         -  cshift ( VJ_I * F_I, shift = +1, dim = iD ) ) &
!                       / ( VJ * dX )

    lV = 1
    where ( shape ( dU ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( dU ) > 1 )
      uV = shape ( dU ) - oV
    end where
      
    iaS = 0
    iaS ( iD ) = +1
      
    !$OMP parallel do private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          dU ( iV, jV, kV ) &
            = dU ( iV, jV, kV )  &
              +  dT * (    VJ_I ( iV, jV, kV ) &
                           *  F_I ( iV, jV, kV ) &
                        -  VJ_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                           *  F_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) &
                      / ( VJ ( iV, jV, kV ) * dX ( iV, jV, kV ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
    
  end subroutine ComputeIncrement_CSL_Kernel


  subroutine RecordBoundaryFluence_CSL &
               ( BF, CSL, F_I, VJ_I, Weight_RK, dT, iD, iC )

    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      BF
    type ( Chart_SLD_Form ), intent ( in ) :: &
      CSL
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_I, &
      VJ_I
    real ( KDR ), intent ( in ) :: &
      Weight_RK, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iD, &  !-- iDimension
      iC     !-- iConserved

    integer ( KDI ) :: &
      jD, kD   !-- jDimension, kDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oB, & !-- oBoundary
      nB    !-- nBoundary
    
    jD = mod ( iD, 3 ) + 1
    kD = mod ( jD, 3 ) + 1
    
    nB ( iD ) = 1
    nB ( jD ) = CSL % nCellsBrick ( jD )
    nB ( kD ) = CSL % nCellsBrick ( kD )

    if ( CSL % iaBrick ( iD ) == 1 ) then
      associate ( iCI => CSL % Atlas % Connectivity % iaInner ( iD ) )
      associate ( BF_Inner => BF ( iC, iCI ) % Value )
      oB = CSL % nGhostLayers
      call RecordBoundaryFluence_CSL_Kernel &
             ( BF_Inner, F_I, VJ_I, Weight_RK * dT, nB, oB )
      end associate !-- BF_Inner
      end associate !-- iCI
    end if !-- iaBrick ( iD ) == 1

    if ( CSL % iaBrick ( iD ) == CSL % nBricks ( iD ) ) then
      associate ( iCO => CSL % Atlas % Connectivity % iaOuter ( iD ) )
      associate ( BF_Outer => BF ( iC, iCO ) % Value )
      oB        = CSL % nGhostLayers
      oB ( iD ) =  oB ( iD ) + CSL % nCellsBrick ( iD )
      call RecordBoundaryFluence_CSL_Kernel &
             ( BF_Outer, F_I, VJ_I, Weight_RK * dT, nB, oB )
      end associate !-- BF_Outer
      end associate !-- iCO
    end if !-- iaBrick ( iD ) == nBricks ( iD )

  end subroutine RecordBoundaryFluence_CSL


  subroutine RecordBoundaryFluence_CSL_Kernel ( BF, F, VJ, Factor, nB, oB )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      BF
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F, &
      VJ
    real ( KDR ), intent ( in ) :: &
      Factor
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nB, &
      oB

    integer ( KDI ) :: &
      iV, jV, kV

    !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
    do kV = 1, nB ( 3 )
      do jV = 1, nB ( 2 )
        do iV = 1, nB ( 1 )
          BF ( iV, jV, kV ) &
            =  BF ( iV, jV, kV ) &
               +  Factor &
                  *  VJ ( oB ( 1 ) + iV, oB ( 2 ) + jV, oB ( 3 ) + kV ) &
                  *   F ( oB ( 1 ) + iV, oB ( 2 ) + jV, oB ( 3 ) + kV )
        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do

  end subroutine RecordBoundaryFluence_CSL_Kernel


end module IncrementDivergence_FV__Form
