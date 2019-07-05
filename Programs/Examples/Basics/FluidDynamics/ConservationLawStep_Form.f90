module ConservationLawStep_Form

  use iso_c_binding
  use Basics
  use ConservedFields_Template
  use ConservationLawStep_Kernel

  implicit none
  private

  type, public :: ConservationLawStepForm
    integer ( KDI ) :: &
      ALPHA_PLUS  = 1, &
      ALPHA_MINUS = 2, &
      N_MODIFIED_SPEEDS = 2, &
      iTimerCommunication, &
      iTimerRKStep, &
      iTimerDifference, &
      iTimerReconstruction, &
      iTimerRiemannSolverInput, &
      iTimerRawFluxes, &
      iTimerFluxes, &
      iTimerPrimitive, &
      iTimerConserved, &
      iTimerAuxiliary, &
      iTimerUpdate, &
      iTimerBoundaryCondition, &
      iTimerDataTransferDevice, &
      iTimerDataTransferHost
    real ( KDR ) :: &
      LimiterParameter
    type ( StorageForm ) :: &
      Old, &
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
  end type ConservationLawStepForm

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
    call CLS % Old % AllocateDevice ( )

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

    call CLS % DifferenceLeft  % AllocateDevice ( )
    call CLS % DifferenceRight % AllocateDevice ( )
    call Clear ( CLS % DifferenceLeft % D_Selected, &
                 CLS % DifferenceLeft % Value )
    call Clear ( CLS % DifferenceRight % D_Selected, &
                 CLS % DifferenceRight % Value )

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
    
    call CLS % ReconstructionInner % AllocateDevice ( )
    call CLS % ReconstructionOuter % AllocateDevice ( )
    call Clear ( CLS % ReconstructionInner % D_Selected, &
                 CLS % ReconstructionInner % Value ) 
    call Clear ( CLS % ReconstructionOuter % D_Selected, &
                 CLS % ReconstructionOuter % Value ) 

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
    
    call CLS % ModifiedSpeedsInner % AllocateDevice ( )
    call CLS % ModifiedSpeedsOuter % AllocateDevice ( )
    call Clear ( CLS % ModifiedSpeedsInner % D_Selected, &
                 CLS % ModifiedSpeedsInner % Value )
    call Clear ( CLS % ModifiedSpeedsOuter % D_Selected, &
                 CLS % ModifiedSpeedsOuter % Value )

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
    call CLS % RawFluxInner % AllocateDevice ( )
    call CLS % RawFluxOuter % AllocateDevice ( )
    call Clear ( CLS % RawFluxInner % D_Selected, &
                 CLS % RawFluxInner % Value )
    call Clear ( CLS % RawFluxOuter % D_Selected, &
                 CLS % RawFluxOuter % Value )

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
    call CLS % FluxInner % AllocateDevice ( )
    call CLS % FluxOuter % AllocateDevice ( )
    call Clear ( CLS % FluxInner % D_Selected, &
                 CLS % FluxInner % Value )
    call Clear ( CLS % FluxOuter % D_Selected, &
                 CLS % FluxOuter % Value )
    
    !-- Update

    call CLS % Update % Initialize &
           ( [ nCells, CF % N_CONSERVED ], NameOption = 'Update', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaConserved ( iV ) ), &
                     iV = 1, CF % N_CONSERVED ) ] )
    call CLS % Update % AllocateDevice ( )

    end associate !-- nCells
    end associate !-- DM
    
    call PROGRAM_HEADER % AddTimer &
           ( 'Communication', CLS % iTimerCommunication, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'RK Step', CLS % iTimerRKStep, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'Update', CLS % iTimerUpdate, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'Difference', CLS % iTimerDifference, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'Reconstruction', CLS % iTimerReconstruction, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'RiemannSolverInput', CLS % iTimerRiemannSolverInput, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'RawFluxes', CLS % iTimerRawFluxes, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'Fluxes', CLS % iTimerFluxes, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputePrimitive', CLS % iTimerPrimitive, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeConserved', CLS % iTimerConserved, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeAuxiliary', CLS % iTimerAuxiliary, Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ApplyBoundaryConditions', CLS % iTimerBoundaryCondition, &
             Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'DataTransfer to Device', CLS % iTimerDataTransferDevice, &
             Level = 2 )
    call PROGRAM_HEADER % AddTimer &
           ( 'DataTransfer to Host', CLS % iTimerDataTransferHost, &
             Level = 2 )
    
  end subroutine Initialize


  subroutine Solve ( CLS, TimeStep )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iV  !-- iVariable
    type ( StorageForm ) :: &
      Primitive!, &
      !Eigenspeed

    associate &
      ( CF => CLS % ConservedFields )
    associate &
      ( DM => CF % DistributedMesh, &
        Current => CF, &
        Old => CLS % Old, &
        Update  => CLS % Update, &
        iaC => CF % iaConserved, &
        T_C  => PROGRAM_HEADER % Timer ( CLS % iTimerCommunication ), &
        T_RK => PROGRAM_HEADER % Timer ( CLS % iTimerRKStep ), &
        T_A  => PROGRAM_HEADER % Timer ( CLS % iTimerAuxiliary ), &
        T_P  => PROGRAM_HEADER % Timer ( CLS % iTimerPrimitive ), &
        T_DT_D  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferDevice ), &
        T_DT_H  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferHost ) )

    call Show ( 'Preparing Step', CONSOLE % INFO_4 )    
    call Primitive % Initialize &
           ( Current, iaSelectedOption = Current % iaPrimitive )
    
    !call Eigenspeed % Initialize &
    !       ( Current, &
    !         iaSelectedOption = [ Current % FAST_EIGENSPEED_PLUS ( : ), &
    !                              Current % FAST_EIGENSPEED_MINUS ( : ) ] )

    call T_DT_D % Start ( )
    call Current % UpdateDevice ( )
    call T_DT_D % Stop ( )
    
    call T_RK % Start ( )
    do iV = 1, Current % N_CONSERVED
      call Copy ( Current % Value ( :, iaC ( iV ) ), &
                  Current % D_Selected ( iaC ( iV ) ), &
                  Old % D_Selected ( iV ), Old % Value ( :, iV ) )
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
               Old % D_Selected ( iV ), Update % D_Selected ( iV ), &
               Current % D_Selected ( iaC ( iV ) ), &
               Current % Value ( :, iaC ( iV ) ) )
    end do
    
    call T_RK % Stop ( )
    
    call Show ( 'Computing Fluid', CONSOLE % INFO_5 )

    call T_P % Start ( )
    call Current % ComputePrimitive ( Current % Value, Current % D_Selected )
    call T_P % Stop ( )
    
    call T_A % Start ( )
    call Current % ComputeAuxiliary ( Current % Value, Current % D_Selected )
    call T_A % Stop ( )
    
    call T_DT_H % Start ( )
    !call Primitive % UpdateHost ( ) 
    call Current % UpdateHost ( )
    call T_DT_H % Stop ( )
    
    call Show ( 'Ghost Exchange', CONSOLE % INFO_5 )

    call T_C % Start ( )
    call DM % StartGhostExchange ( )
    call DM % FinishGhostExchange ( )
    call T_C % Stop ( )
    
    call T_DT_D % Start ( )
    !call Primitive % UpdateDevice ( )
    call Current % UpdateDevice ( )
    call T_DT_D % Stop ( )
    
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
               Current % D_Selected ( iaC ( iV ) ), &
               Old % D_Selected ( iV ), Update % D_Selected ( iV ) )
    end do
    call T_RK % Stop ( )
    
    call Show ( 'Computing Fluid', CONSOLE % INFO_5 )

    call T_P % Start ( )
    call Current % ComputePrimitive ( Current % Value, Current % D_Selected )
    call T_P % Stop ( )
    
    call T_A % Start ( )
    call Current % ComputeAuxiliary ( Current % Value, Current % D_Selected )
    call T_A % Stop ( )
    
    call T_DT_H % Start ( )
    !call Primitive % UpdateHost ( ) 
    call Current % UpdateHost ( ) 
    call T_DT_H % Stop ( )
    
    call Show ( 'Ghost Exchange', CONSOLE % INFO_5 )

    call T_C % Start ( )
    call DM % StartGhostExchange ( )
    call DM % FinishGhostExchange ( )
    call T_C % Stop ( )
    
    call T_DT_H % Start ( )
    !call Primitive % UpdateHost ( ) 
    call Current % UpdateDevice ( ) 
    call T_DT_H % Stop ( )
    
    !call T_DT_H % Start ( )
    !call Eigenspeed % UpdateHost ( ) 
    !call T_DT_H % Stop ( )
    
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

    call Show ( 'Computing Update', CONSOLE % INFO_5 )

    associate ( CF => CLS % ConservedFields )
    associate ( DM => CF % DistributedMesh )
    associate ( T_U => PROGRAM_HEADER % Timer ( CLS % iTimerUpdate ) )

    call Clear ( CLS % Update % D_Selected, CLS % Update % Value )
    
    do iD = 1, DM % nDimensions

      call Show ( iD, 'iD', CONSOLE % INFO_5 )

      call CLS % ComputeDifferences ( iD )
      call CLS % ComputeReconstruction ( )
      call CLS % ComputeFluxes ( iD )
      
      call T_U % Start ( )
      do iV = 1, CF % N_CONSERVED
        call ComputeUpdateKernel &
               ( CLS % Update % Value ( :, iV ), &
                 CLS % FluxInner % Value ( :, iV ), &
                 CLS % FluxOuter % Value ( :, iV ), DM % CellVolume, &
                 DM % CellArea ( iD ), TimeStep, &
                 CLS % Update % D_Selected ( iV ), &
                 CLS % FluxInner % D_Selected ( iV ), &
                 CLS % FluxOuter % D_Selected ( iV ) )
      end do
      call T_U % Stop ( )

    end do
    
    end associate !-- T_DT_D, etc.
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

    associate ( CF => CLS % ConservedFields )
    associate ( DM => CF % DistributedMesh )
    associate &
      ( T_DT_D  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferDevice ), &
        T_DT_H  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferHost ), &
        T_D     => PROGRAM_HEADER % Timer ( CLS % iTimerDifference ), &
        T_BC    => PROGRAM_HEADER % Timer ( CLS % iTimerBoundaryCondition ) )
    
    call T_BC % Start ( )
    call CF % ApplyBoundaryConditions &
           ( CF % Value, CF % Value, iD, iBoundary = -1, &
             D_ExteriorValue = CF % D_Selected, &
             D_InteriorValue = CF % D_Selected )
    call CF % ApplyBoundaryConditions &
           ( CF % Value, CF % Value, iD, iBoundary = +1, &
             D_ExteriorValue = CF % D_Selected, &
             D_InteriorValue = CF % D_Selected )
    call T_BC % Stop ( )
           
    do iP = 1, CF % N_PRIMITIVE
      call DM % SetVariablePointer &
             ( CF % Value ( :, CF % iaPrimitive ( iP ) ), V )
      call DM % SetVariablePointer &
             ( CLS % DifferenceLeft % Value ( :, iP ), dV_Left )
      call DM % SetVariablePointer &
             ( CLS % DifferenceRight % Value ( :, iP ), dV_Right )
      call T_D % Start ( )
      call ComputeDifferencesKernel &
             ( V, DM % nGhostLayers ( iD ), iD, &
               CF % D_Selected ( CF % iaPrimitive ( iP ) ), &
               CLS % DifferenceLeft % D_Selected ( iP ), &
               CLS % DifferenceRight % D_Selected ( iP ), &
               dV_Left, dV_Right )
      call T_D % Stop ( )
    end do
    
    nullify ( V, dV_Left, dV_Right )
    
    end associate !-- T_DT
    end associate !-- DM
    end associate !-- CF

  end subroutine ComputeDifferences


  subroutine ComputeReconstruction ( CLS )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS

    integer ( KDI ) :: &
      iP  !-- iPrimitive

    associate ( CF => CLS % ConservedFields )
    associate ( iaP => CF % iaPrimitive )
    
    associate &
      ( T_DT_D  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferDevice ), &
        T_DT_H  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferHost ), &
        T_R     => PROGRAM_HEADER % Timer ( CLS % iTimerReconstruction ), &
        T_A     => PROGRAM_HEADER % Timer ( CLS % iTimerAuxiliary ), &
        T_C     => PROGRAM_HEADER % Timer ( CLS % iTimerConserved ) )
    
    do iP = 1, CF % N_PRIMITIVE
      call T_R % Start ( )
      call ComputeReconstructionKernel &
             ( CF % Value ( :, iaP ( iP ) ), &
               CLS % DifferenceLeft % Value ( :, iP ), &
               CLS % DifferenceRight % Value ( :, iP ), &
               CLS % LimiterParameter, &
               CF % D_Selected ( iaP ( iP ) ), &
               CLS % DifferenceLeft % D_Selected ( iP ), &
               CLS % DifferenceRight % D_Selected ( iP ), &
               CLS % ReconstructionInner % D_Selected ( iaP ( iP ) ), &
               CLS % ReconstructionOuter % D_Selected ( iaP ( iP ) ), &
               CLS % ReconstructionInner % Value ( :, iaP ( iP ) ), &
               CLS % ReconstructionOuter % Value ( :, iaP ( iP ) ) )
      call T_R % Stop ( )
    end do

    call T_A % Start ( )
    call CF % ComputeAuxiliary &
           ( CLS % ReconstructionInner % Value, &
             CLS % ReconstructionInner % D_Selected )
    call CF % ComputeAuxiliary &
           ( CLS % ReconstructionOuter % Value, &
             CLS % ReconstructionOuter % D_Selected )
    call T_A % Stop ( )
    
    call T_C % Start ( )
    call CF % ComputeConserved &
           ( CLS % ReconstructionInner % Value, &
             CLS % ReconstructionInner % D_Selected )
    call CF % ComputeConserved &
           ( CLS % ReconstructionOuter % Value, &
             CLS % ReconstructionOuter % D_Selected )
    call T_C % Stop ( )
    
    end associate !-- Timer
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

    associate ( CF => CLS % ConservedFields )
    associate &
      ( DM  => CF % DistributedMesh, &
        iaC => CF % iaConserved, &
        T_DT_D  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferDevice ), &
        T_DT_H  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferHost ), &
        T_RSI => PROGRAM_HEADER % Timer ( CLS % iTimerRiemannSolverInput ), &
        T_RF  => PROGRAM_HEADER % Timer ( CLS % iTimerRawFluxes ), &
        T_F   => PROGRAM_HEADER % Timer ( CLS % iTimerFluxes ), &
        T_BC  => PROGRAM_HEADER % Timer ( CLS % iTimerBoundaryCondition ) )

    call T_BC % Start ( )
    call CF % ApplyBoundaryConditions &
           ( CLS % ReconstructionOuter % Value, &
             CLS % ReconstructionInner % Value, iDimension, iBoundary = -1, &
             D_ExteriorValue = CLS % ReconstructionOuter % D_Selected, &
             D_InteriorValue = CLS % ReconstructionInner % D_Selected )
    call CF % ApplyBoundaryConditions &
           ( CLS % ReconstructionInner % Value, &
             CLS % ReconstructionOuter % Value, iDimension, iBoundary = +1, &
             D_ExteriorValue = CLS % ReconstructionInner % D_Selected, &
             D_InteriorValue = CLS % ReconstructionOuter % D_Selected )
    call T_BC % Stop ( )

    call T_RSI % Start ( )
    call CF % ComputeRiemannSolverInput &
           ( CLS, CLS % ReconstructionInner % Value, &
             CLS % ReconstructionOuter % Value, iDimension, &
             CLS % ReconstructionInner % D_Selected, &
             CLS % ReconstructionOuter % D_Selected )
    call T_RSI % Stop ( )

    call T_RF % Start ( )
    call CF % ComputeRawFluxes &
           ( CLS % RawFluxInner % Value, CLS % ReconstructionInner % Value, &
             iDimension, CLS % RawFluxInner % D_Selected, &
             CLS % ReconstructionInner % D_Selected )
    call CF % ComputeRawFluxes &
           ( CLS % RawFluxOuter % Value, CLS % ReconstructionOuter % Value, &
             iDimension, CLS % RawFluxOuter % D_Selected, &
             CLS % ReconstructionOuter % D_Selected )
    call T_RF % Stop ( )
    
    call DM % SetVariablePointer &
           ( CLS % ModifiedSpeedsInner % Value ( :, CLS % ALPHA_PLUS ), AP_I )
    call DM % SetVariablePointer &
           ( CLS % ModifiedSpeedsOuter % Value ( :, CLS % ALPHA_PLUS ), AP_O )
    call DM % SetVariablePointer &
           ( CLS % ModifiedSpeedsInner % Value ( :, CLS % ALPHA_MINUS ), AM_I )
    call DM % SetVariablePointer &
           ( CLS % ModifiedSpeedsOuter % Value ( :, CLS % ALPHA_MINUS ), AM_O )

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
               DM % nGhostLayers ( iDimension ), iDimension, &
               CLS % ModifiedSpeedsInner % D_Selected ( CLS % ALPHA_PLUS ), &
               CLS % ModifiedSpeedsOuter % D_Selected ( CLS % ALPHA_PLUS ), &
               CLS % ModifiedSpeedsInner % D_Selected ( CLS % ALPHA_MINUS ), &
               CLS % ModifiedSpeedsOuter % D_Selected ( CLS % ALPHA_MINUS ), &
               CLS % RawFluxInner % D_Selected ( iV ), &
               CLS % RawFluxOuter % D_Selected ( iV ), &
               CLS % ReconstructionInner % D_Selected ( iaC ( iV ) ), &
               CLS % ReconstructionOuter % D_Selected ( iaC ( iV ) ), &
               CLS % FluxInner % D_Selected ( iV ), &
               CLS % FluxOuter % D_Selected ( iV ), &
               F_I, F_O )
      call T_F % Stop ( )
    end do
    
    nullify ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, U_I, U_O, F_I, F_O )
    end associate !-- DM, iaC
    end associate !-- CF

  end subroutine ComputeFluxes


end module ConservationLawStep_Form
