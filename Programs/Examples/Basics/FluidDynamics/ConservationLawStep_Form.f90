module ConservationLawStep_Form

  use Basics
  use ConservedFields_Template

  implicit none
  private

  type, public :: ConservationLawStepForm
    integer ( KDI ) :: &
      ALPHA_PLUS  = 1, &
      ALPHA_MINUS = 2, &
      N_MODIFIED_SPEEDS = 2, &
      iTimerCommunication = 0, &
      iTimerRKStep = 0, &
      iTimerDifference = 0, &
      iTimerReconstruction = 0, &
      iTimerRiemannSolverInput = 0, &
      iTimerRawFluxes = 0, &
      iTimerFluxes = 0, &
      iTimerPrimitive = 0, &
      iTimerConserved = 0, &
      iTimerAuxiliary = 0, &
      iTimerUpdate = 0, &
      iTimerBoundaryCondition = 0, &
      iTimerDataTransferDevice = 0, &
      iTimerDataTransferHost = 0
    real ( KDR ) :: &
      LimiterParameter
    type ( StorageForm ) :: &
      Old, &
      Primitive, &
      DifferenceLeft, &
      DifferenceRight, &
      ReconstructionInner, &
      ReconstructionOuter, &
      ModifiedSpeedsInner, &
      ModifiedSpeedsOuter, &
      RawFluxInner, &
      RawFluxOuter, &
      FluxInner, &
      FluxOuter, &
      Update
    class ( ConservedFieldsTemplate ), pointer :: &
      ConservedFields => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Solve
    procedure, private, pass :: &
      ComputeUpdate
    procedure, private, pass :: &
      ComputeDifferences
    procedure, private, pass :: &
      ComputeReconstruction
    procedure, private, pass :: &
      ComputeFluxes
    procedure, private, pass :: &
      PrepareTimers
  end type ConservationLawStepForm
  
    private :: &
      ComputeDifferencesKernel, &
      ComputeReconstructionKernel, &
      ComputeFluxesKernel, &
      AddUpdateKernel, &
      ComputeUpdateKernel, &
      CombineUpdatesKernel
      
      
    interface
    
      module subroutine ComputeDifferencesKernel &
                   ( V, oV, iD, dV_Left, dV_Right, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          V
        integer ( KDI ), intent ( in ) :: &
          oV, &
          iD
        real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
          dV_Left, dV_Right
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeDifferencesKernel

      module subroutine ComputeReconstructionKernel &
                   ( V, dV_Left, dV_Right, Theta, V_Inner, V_Outer, &
                     UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          V, &
          dV_Left, dV_Right
        real ( KDR ), intent ( in ) :: &
          Theta
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          V_Inner, V_Outer
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeReconstructionKernel

      module subroutine ComputeFluxesKernel &
                   ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, U_I, U_O, oV, iD, &
                     F_I, F_O, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          AP_I, AP_O, &
          AM_I, AM_O, &
          RF_I, RF_O, &
          U_I, U_O
        integer ( KDI ), intent ( in ) :: &
          oV, &
          iD
        real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
          F_I, F_O
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeFluxesKernel

      module subroutine ComputeUpdateKernel &
                   ( dU, F_I, F_O, V, A, dT, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          dU
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          F_I, F_O
        real ( KDR ), intent ( in ) :: &
          V, &
          A, &
          dT
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeUpdateKernel
      
      module subroutine AddUpdateKernel ( O, U, C, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          O, &
          U
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          C
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine AddUpdateKernel
      
      module subroutine CombineUpdatesKernel ( C, O, U, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          C 
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          O, &
          U
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine CombineUpdatesKernel
      
    end interface

contains


  subroutine Initialize ( CLS, CF )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    class ( ConservedFieldsTemplate ), intent ( in ), target :: &
      CF

    integer ( KDI ) :: &
      iV  !-- iVariable
      
    call Show ( 'Initializing a ConservationLawStep', CONSOLE % INFO_3 )

    CLS % ConservedFields => CF

    associate ( DM => CF % DistributedMesh )
    associate ( nCells => DM % nProperCells + DM % nGhostCells )
    
    !-- Old
    
    call CLS % Old % Initialize ( [ nCells, CF % N_CONSERVED ] )
    if ( CF % AllocatedDevice ) &
      call CLS % Old % AllocateDevice ( )
    
    !-- Primitive, an overlay storage
    call CLS % Primitive % Initialize &
           ( CF, iaSelectedOption = CF % iaPrimitive )

    !-- Differences

    call CLS % DifferenceLeft % Initialize &
           ( [ nCells, CF % N_PRIMITIVE ], NameOption = 'DifferenceLeft', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaPrimitive ( iV ) ), &
                     iV = 1, CF % N_PRIMITIVE ) ] )
    call CLS % DifferenceRight % Initialize &
           ( [ nCells, CF % N_PRIMITIVE ], NameOption = 'DifferenceRight', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaPrimitive ( iV ) ), &
                     iV = 1, CF % N_PRIMITIVE ) ] )
    
    if ( CF % AllocatedDevice ) then
      call CLS % DifferenceLeft  % AllocateDevice ( )
      call CLS % DifferenceRight % AllocateDevice ( )
    end if
    call Clear ( CLS % DifferenceLeft % Value, &
                 UseDeviceOption = CF % AllocatedDevice )
    call Clear ( CLS % DifferenceRight % Value, &
                 UseDeviceOption = CF % AllocatedDevice )

    !-- Reconstruction

    CLS % LimiterParameter = 1.4_KDR
    call PROGRAM_HEADER % GetParameter &
           ( CLS % LimiterParameter, 'LimiterParameter' )

    call CLS % ReconstructionInner % Initialize &
           ( [ nCells, CF % nVariables ], NameOption = 'ReconstructionInner', &
             VariableOption = CF % Variable )
    call CLS % ReconstructionOuter % Initialize &
           ( [ nCells, CF % nVariables ], NameOption = 'ReconstructionOuter', &
             VariableOption = CF % Variable )
    
    if ( CF % AllocatedDevice ) then
      call CLS % ReconstructionInner % AllocateDevice ( )
      call CLS % ReconstructionOuter % AllocateDevice ( )
    end if 
    call Clear ( CLS % ReconstructionInner % Value, &
                 UseDeviceOption = CF % AllocatedDevice ) 
    call Clear ( CLS % ReconstructionOuter % Value, &
                 UseDeviceOption = CF % AllocatedDevice ) 

    !-- Riemann solver

    call CLS % ModifiedSpeedsInner % Initialize &
           ( [ nCells, CLS % N_MODIFIED_SPEEDS ], &
             NameOption = 'ModifiedSpeedsInner', &
             VariableOption &
               = [ 'AlphaPlus                      ', &
                   'AlphaMinus                     ' ] )
    call CLS % ModifiedSpeedsOuter % Initialize &
           ( [ nCells, CLS % N_MODIFIED_SPEEDS ], &
             NameOption = 'ModifiedSpeedsInner', &
             VariableOption &
               = [ 'AlphaPlus                      ', &
                   'AlphaMinus                     ' ] )
    
    if ( CF % AllocatedDevice ) then
      call CLS % ModifiedSpeedsInner % AllocateDevice ( )
      call CLS % ModifiedSpeedsOuter % AllocateDevice ( )
    end if
    call Clear ( CLS % ModifiedSpeedsInner % Value, &
                 UseDeviceOption = CF % AllocatedDevice )
    call Clear ( CLS % ModifiedSpeedsOuter % Value, &
                 UseDeviceOption = CF % AllocatedDevice )

    call CLS % RawFluxInner % Initialize &
           ( [ nCells, CF % N_CONSERVED ], NameOption = 'RawFluxInner', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaConserved ( iV ) ), &
                     iV = 1, CF % N_CONSERVED ) ] )
    call CLS % RawFluxOuter % Initialize &
           ( [ nCells, CF % N_CONSERVED ], NameOption = 'RawFluxOuter', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaConserved ( iV ) ), &
                     iV = 1, CF % N_CONSERVED ) ] )
    if ( CF % AllocatedDevice ) then
      call CLS % RawFluxInner % AllocateDevice ( )
      call CLS % RawFluxOuter % AllocateDevice ( )
    end if
    call Clear ( CLS % RawFluxInner % Value, &
                 UseDeviceOption = CF % AllocatedDevice )
    call Clear ( CLS % RawFluxOuter % Value, &
                 UseDeviceOption = CF % AllocatedDevice )

    call CLS % FluxInner % Initialize &
           ( [ nCells, CF % N_CONSERVED ], NameOption = 'FluxInner', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaConserved ( iV ) ), &
                     iV = 1, CF % N_CONSERVED ) ] )
    call CLS % FluxOuter % Initialize &
           ( [ nCells, CF % N_CONSERVED ], NameOption = 'FluxOuter', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaConserved ( iV ) ), &
                     iV = 1, CF % N_CONSERVED ) ] )
    if ( CF % AllocatedDevice ) then
      call CLS % FluxInner % AllocateDevice ( )
      call CLS % FluxOuter % AllocateDevice ( )
    end if
    call Clear ( CLS % FluxInner % Value, &
                 UseDeviceOption = CF % AllocatedDevice )
    call Clear ( CLS % FluxOuter % Value, &
                 UseDeviceOption = CF % AllocatedDevice )
    
    !-- Update

    call CLS % Update % Initialize &
           ( [ nCells, CF % N_CONSERVED ], NameOption = 'Update', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaConserved ( iV ) ), &
                     iV = 1, CF % N_CONSERVED ) ] )
    if ( CF % AllocatedDevice ) &
      call CLS % Update % AllocateDevice ( )

    end associate !-- nCells
    
    end associate !-- DM
    
  end subroutine Initialize


  subroutine Solve ( CLS, TimeStep )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iV  !-- iVariable
    type ( TimerForm ), pointer :: &
      T_RK, &
      T_P, &
      T_A, &
      T_DT_H, &
      T_C, &
      T_DT_D

    associate &
      ( CF => CLS % ConservedFields )
    associate &
      ( DM          => CF % DistributedMesh, &
        Current     => CF, &
        Old         => CLS % Old, &
        Primitive   => CLS % Primitive, &
        Update      => CLS % Update, &
        iaC  => CF % iaConserved )
    
    call CLS % PrepareTimers ( )
    
    call Show ( 'Preparing Step', CONSOLE % INFO_4 )
    
    T_RK  =>  PROGRAM_HEADER % Timer &
                ( CLS % iTimerRKStep, 'RK Step', Level = 2 )
    call T_RK % Start ( )
    do iV = 1, Current % N_CONSERVED
      call Copy ( Current % Value ( :, iaC ( iV ) ), Old % Value ( :, iV ), &
                  UseDeviceOption = Current % AllocatedDevice )
    end do
    call T_RK % Stop ( )
    

    !-- Substep 1
    
    call Show ( 'Solving Substep 1', CONSOLE % INFO_4 )

    call CLS % ComputeUpdate ( TimeStep ) !-- K1 = dT * RHS
    
    call Show ( 'Adding Update', CONSOLE % INFO_5 )
    call T_RK % Start ( )
    do iV = 1, Current % N_CONSERVED
      !-- Current = Old + K1
    !  Current % Value ( :, iaC ( iV ) ) &
    !    = Old % Value ( :, iV ) + Update % Value ( :, iV )
      call AddUpdateKernel &
             ( Old % Value ( :, iV ), Update % Value ( :, iV ), &
               Current % Value ( :, iaC ( iV ) ), &
               UseDeviceOption = Current % AllocatedDevice )
    end do
    
    call T_RK % Stop ( )
    
    call Show ( 'Computing Fluid', CONSOLE % INFO_5 )

    T_P  =>  PROGRAM_HEADER % Timer &
               ( CLS % iTimerPrimitive, 'Compute Primitive', Level = 2 )
    call T_P % Start ( )
    call Current % ComputePrimitive ( Current % Value )
    call T_P % Stop ( )
    
    T_A  =>  PROGRAM_HEADER % Timer &
               ( CLS % iTimerAuxiliary, 'Compute Auxiliary', Level = 2 )
    call T_A % Start ( )
    call Current % ComputeAuxiliary ( Current % Value )
    call T_A % Stop ( )
    
     T_DT_H  =>  PROGRAM_HEADER % Timer &
                  ( CLS % iTimerDataTransferHost, 'DataTransfer to Host', &
                    Level = 2 )
    if ( .not. DM % DevicesCommunicate ) then
      call T_DT_H % Start ( )
      call Primitive % UpdateHost ( ) 
      call T_DT_H % Stop ( )
    end if

    T_C  =>  PROGRAM_HEADER % Timer &
               ( CLS % iTimerCommunication, 'Communication', Level = 2 )
    call T_C % Start ( )
    call DM % StartGhostExchange ( )
    call DM % FinishGhostExchange ( )
    call T_C % Stop ( )
    
    T_DT_D  =>  PROGRAM_HEADER % Timer &
                  ( CLS % iTimerDataTransferDevice, 'DataTransfer to Device', &
                    Level = 2 )
    if ( .not. DM % DevicesCommunicate ) then
      call T_DT_D % Start ( )
      call Primitive % UpdateDevice ( )
      call T_DT_D % Stop ( )
    end if
    
    !-- Substep 2
    
    call Show ( 'Solving Substep 2', CONSOLE % INFO_4 )

    call CLS % ComputeUpdate ( TimeStep ) !-- K2 = dT * RHS
    
    call Show ( 'Combining Updates', CONSOLE % INFO_5 )
    call T_RK % Start ( )
    do iV = 1, Current % N_CONSERVED
      !-- Current = Old + 0.5 * ( k1 + k2 )
      !           = 0.5 Old + 0.5 ( Old + k1 + k2 )            
      !Current % Value ( :, iaC ( iV ) ) &
      !  = Current % Value ( :, iaC ( iV ) ) + Update % Value ( :, iV )
      !Current % Value ( :, iaC ( iV ) ) &
      !  = 0.5_KDR * ( Old % Value ( :, iV ) &
      !                + Current % Value ( :, iaC ( iV ) ) )
      call CombineUpdatesKernel &
             ( Current % Value ( :, iaC ( iV ) ), &
               Old % Value ( :, iV ), Update % Value ( :, iV ), &
               UseDeviceOption = Current % AllocatedDevice )
    end do
    call T_RK % Stop ( )
    
    call Show ( 'Computing Fluid', CONSOLE % INFO_5 )

    call T_P % Start ( )
    call Current % ComputePrimitive ( Current % Value )
    call T_P % Stop ( )
    
    call T_A % Start ( )
    call Current % ComputeAuxiliary ( Current % Value )
    call T_A % Stop ( )
    
    if ( .not. DM % DevicesCommunicate ) then
      call T_DT_H % Start ( )
      call Primitive % UpdateHost ( ) 
      call T_DT_H % Stop ( )
    end if
    
    call T_C % Start ( )
    call DM % StartGhostExchange ( )
    call DM % FinishGhostExchange ( )
    call T_C % Stop ( )
    
    if ( .not. DM % DevicesCommunicate ) then
      call T_DT_D % Start ( )
      call Primitive % UpdateDevice ( ) 
      call T_DT_D % Stop ( )
    end if
    
    end associate !-- DM, etc.
    end associate !-- CF

  end subroutine Solve


  subroutine ComputeUpdate ( CLS, TimeStep )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      iV
    type ( TimerForm ), pointer :: &
      T_U

    call Show ( 'Computing Update', CONSOLE % INFO_5 )
    
    associate ( CF => CLS % ConservedFields )
    associate ( DM => CF % DistributedMesh )

    call Clear ( CLS % Update % Value, UseDeviceOption = CF % AllocatedDevice )
    
    do iD = 1, DM % nDimensions

      call Show ( iD, 'iD', CONSOLE % INFO_5 )

      call CLS % ComputeDifferences ( iD )
      call CLS % ComputeReconstruction ( )
      call CLS % ComputeFluxes ( iD )
      
      T_U  =>  PROGRAM_HEADER % Timer &
                 ( CLS % iTimerUpdate, 'Update', Level = 2 )
      call T_U % Start ( )
      do iV = 1, CF % N_CONSERVED
        call ComputeUpdateKernel &
               ( CLS % Update % Value ( :, iV ), &
                 CLS % FluxInner % Value ( :, iV ), &
                 CLS % FluxOuter % Value ( :, iV ), DM % CellVolume, &
                 DM % CellArea ( iD ), TimeStep, &
                 UseDeviceOption = CF % AllocatedDevice )
      end do
      call T_U % Stop ( )

    end do
    
    end associate !-- DM
    end associate !-- CF

  end subroutine ComputeUpdate


  subroutine ComputeDifferences ( CLS, iD )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    integer ( KDI ), intent ( in ) :: &
      iD

    integer ( KDI ) :: &
      iP  !-- iPrimitive
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V, &
      dV_Left, &
      dV_Right
    type ( TimerForm ), pointer :: &
      T_BC, &
      T_D

    associate ( CF => CLS % ConservedFields )
    associate ( DM => CF % DistributedMesh )
    
    T_BC  =>  PROGRAM_HEADER % Timer &
                ( CLS % iTimerBoundaryCondition, 'ApplyBoundaryConditions', &
                  Level = 2 )
    call T_BC % Start ( )
    call CF % ApplyBoundaryConditions &
           ( CF % Value, CF % Value, iD, iBoundary = -1, &
             PrimitiveOnlyOption = .true. )
    call CF % ApplyBoundaryConditions &
           ( CF % Value, CF % Value, iD, iBoundary = +1, &
             PrimitiveOnlyOption = .true. )
    call T_BC % Stop ( )
           
    T_D  =>  PROGRAM_HEADER % Timer &
                ( CLS % iTimerDifference, 'Difference', &
                  Level = 2 )
    do iP = 1, CF % N_PRIMITIVE
      call DM % SetVariablePointer &
             ( CF % Value ( :, CF % iaPrimitive ( iP ) ), V )
      call DM % SetVariablePointer &
             ( CLS % DifferenceLeft % Value ( :, iP ), dV_Left )
      call DM % SetVariablePointer &
             ( CLS % DifferenceRight % Value ( :, iP ), dV_Right )
      call T_D % Start ( )
      call ComputeDifferencesKernel &
             ( V, DM % nGhostLayers ( iD ), iD, dV_Left, dV_Right, &
               UseDeviceOption = CF % AllocatedDevice )
      call T_D % Stop ( )
    end do
    
    nullify ( V, dV_Left, dV_Right )
    
    end associate !-- DM
    end associate !-- CF

  end subroutine ComputeDifferences


  subroutine ComputeReconstruction ( CLS )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS

    integer ( KDI ) :: &
      iP  !-- iPrimitive
    type ( TimerForm ), pointer :: &
      T_R, &
      T_A, &
      T_C

    associate ( CF => CLS % ConservedFields )
    associate ( iaP => CF % iaPrimitive )
    
    T_R  =>  PROGRAM_HEADER % Timer &
               ( CLS % iTimerReconstruction, 'Reconstruction', Level = 2 )
    do iP = 1, CF % N_PRIMITIVE
      call T_R % Start ( )
      call ComputeReconstructionKernel &
             ( CF % Value ( :, iaP ( iP ) ), &
               CLS % DifferenceLeft % Value ( :, iP ), &
               CLS % DifferenceRight % Value ( :, iP ), &
               CLS % LimiterParameter, &
               CLS % ReconstructionInner % Value ( :, iaP ( iP ) ), &
               CLS % ReconstructionOuter % Value ( :, iaP ( iP ) ), &
               UseDeviceOption = CF % AllocatedDevice )
      call T_R % Stop ( )
    end do

    T_A  =>  PROGRAM_HEADER % Timer &
               ( CLS % iTimerAuxiliary, 'ComputeAuxiliary', Level = 2 )
    call T_A % Start ( )
    call CF % ComputeAuxiliary &
           ( CLS % ReconstructionInner % Value )
    call CF % ComputeAuxiliary &
           ( CLS % ReconstructionOuter % Value )             
    call T_A % Stop ( )
    
    T_C  =>  PROGRAM_HEADER % Timer &
               ( CLS % iTimerConserved, 'ComputeConserved', Level = 2 )
    call T_C % Start ( )
    call CF % ComputeConserved &
           ( CLS % ReconstructionInner % Value )
    call CF % ComputeConserved &
           ( CLS % ReconstructionOuter % Value )             
    call T_C % Stop ( )
    
    end associate !-- iaP
    end associate !-- CF

  end subroutine ComputeReconstruction


  subroutine ComputeFluxes ( CLS, iDimension )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iV  !-- iVariable
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      AP_I, AP_O, &
      AM_I, AM_O, &
      RF_I, RF_O, &
      U_I, U_O, &
      F_I, F_O
    type ( TimerForm ), pointer :: &
      T_BC, &
      T_F, &
      T_RF, &
      T_RSI

    associate ( CF => CLS % ConservedFields )
    associate &
      ( DM  => CF % DistributedMesh, &
        iaC => CF % iaConserved )

    T_BC  =>  PROGRAM_HEADER % Timer &
                ( CLS % iTimerBoundaryCondition, 'ApplyBoundaryConditions', &
                  Level = 2 )
    call T_BC % Start ( )
    call CF % ApplyBoundaryConditions &
           ( CLS % ReconstructionOuter % Value, &
             CLS % ReconstructionInner % Value, iDimension, iBoundary = -1 )
             
    call CF % ApplyBoundaryConditions &
           ( CLS % ReconstructionInner % Value, &
             CLS % ReconstructionOuter % Value, iDimension, iBoundary = +1 )
    call T_BC % Stop ( )

    T_RSI  =>  PROGRAM_HEADER % Timer &
                 ( CLS % iTimerRiemannSolverInput, 'RiemannSolverInput', &
                   Level = 2 )
    call T_RSI % Start ( )
    call CF % ComputeRiemannSolverInput &
           ( CLS, CLS % ReconstructionInner % Value, &
             CLS % ReconstructionOuter % Value, iDimension )
    call T_RSI % Stop ( )

    T_RF  =>  PROGRAM_HEADER % Timer &
                ( CLS % iTimerRawFluxes, 'RawFluxes', &
                  Level = 2 )
    call T_RF % Start ( )
    call CF % ComputeRawFluxes &
           ( CLS % RawFluxInner % Value, CLS % ReconstructionInner % Value, &
             iDimension )
    call CF % ComputeRawFluxes &
           ( CLS % RawFluxOuter % Value, CLS % ReconstructionOuter % Value, &
             iDimension )
    call T_RF % Stop ( )
    
    call DM % SetVariablePointer &
           ( CLS % ModifiedSpeedsInner % Value ( :, CLS % ALPHA_PLUS ), AP_I )
    call DM % SetVariablePointer &
           ( CLS % ModifiedSpeedsOuter % Value ( :, CLS % ALPHA_PLUS ), AP_O )
    call DM % SetVariablePointer &
           ( CLS % ModifiedSpeedsInner % Value ( :, CLS % ALPHA_MINUS ), AM_I )
    call DM % SetVariablePointer &
           ( CLS % ModifiedSpeedsOuter % Value ( :, CLS % ALPHA_MINUS ), AM_O )

    T_F  =>  PROGRAM_HEADER % Timer &
               ( CLS % iTimerFluxes, 'Fluxes', &
                 Level = 2 )
    do iV = 1, CF % N_CONSERVED
      call DM % SetVariablePointer ( CLS % RawFluxInner % Value ( :, iV ), RF_I)
      call DM % SetVariablePointer ( CLS % RawFluxOuter % Value ( :, iV ), RF_O)
      call DM % SetVariablePointer ( CLS % FluxInner % Value ( :, iV ), F_I )
      call DM % SetVariablePointer ( CLS % FluxOuter % Value ( :, iV ), F_O )
      call DM % SetVariablePointer &
             ( CLS % ReconstructionInner % Value ( :, iaC ( iV ) ), U_I )
      call DM % SetVariablePointer &
             ( CLS % ReconstructionOuter % Value ( :, iaC ( iV ) ), U_O )
      call T_F % Start ( )
      call ComputeFluxesKernel &
             ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, U_I, U_O, &
               DM % nGhostLayers ( iDimension ), iDimension, F_I, F_O, &
               UseDeviceOption = CF % AllocatedDevice )
      call T_F % Stop ( )
    end do
    
    nullify ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, U_I, U_O, F_I, F_O )
    end associate !-- DM, iaC
    end associate !-- CF

  end subroutine ComputeFluxes
  
  
  subroutine PrepareTimers ( CLS )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
      
    type ( TimerForm ), pointer :: &
      T
      
    associate ( DM => CLS % ConservedFields % DistributedMesh )
    
    !-- Force Timers creation in certain order
    
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerCommunication, 'Communication', Level = 2 )
    
    if ( CONSOLE % Verbosity >= CONSOLE % INFO_3 ) then
      T => PROGRAM_HEADER % Timer &
             ( DM % iTimerPacking, 'Pack/Unpack', Level = 3 )
      T => PROGRAM_HEADER % Timer &
             ( DM % iTimerComm, 'Send/Recv', Level = 3 )
    end if
    
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerRKStep, 'RK Step', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerUpdate, 'Update', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerDifference, 'Difference', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerReconstruction, 'Reconstruction', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerRiemannSolverInput, 'RiemannSolverInput', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerRawFluxes, 'RawFluxes', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerFluxes, 'Fluxes', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerPrimitive, 'ComputePrimitive', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerConserved, 'ComputeConserved', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerAuxiliary, 'ComputeAuxiliary', Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerBoundaryCondition, 'ApplyBoundaryConditions', &
              Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerDataTransferDevice, 'DataTransfer to Device', &
             Level = 2 )
    T => PROGRAM_HEADER % Timer &
           ( CLS % iTimerDataTransferHost, 'DataTransfer to Host', &
             Level = 2 )
    
    end associate !-- DM
      
  end subroutine PrepareTimers


end module ConservationLawStep_Form
