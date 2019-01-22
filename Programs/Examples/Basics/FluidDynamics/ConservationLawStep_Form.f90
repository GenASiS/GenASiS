
#include "Preprocessor"

module ConservationLawStep_Form

  use iso_c_binding
  use Basics
  use ConservedFields_Template

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

    private :: &
      ComputeDifferencesKernel, &
      ComputeReconstructionKernel, &
      ComputeFluxesKernel, &
      AddUpdateKernel, &
      CombineUpdatesKernel

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


  subroutine ComputeDifferencesKernel &
               ( V, oV, iD, D_V, D_dV_Left, D_dV_Right, dV_Left, dV_Right )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      V
    integer ( KDI ), intent ( in ) :: &
      oV, &
      iD
    type ( c_ptr ), intent ( in ) :: &
      D_V, &
      D_dV_Left, D_dV_Right
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dV_Left, dV_Right
      
    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV
      
    call AssociateHost ( D_V, V )
    call AssociateHost ( D_dV_Left, dV_Left )
    call AssociateHost ( D_dV_Right, dV_Right )

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV

    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD ) - 1
    
!    dV_Left  = V - cshift ( V, shift = -1, dim = iD )    

    iaS = 0
    iaS ( iD ) = -1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
        
          !call Show ( [ iV, jV, kV ], '<<< iV, jV, kV: 1' )

          iaVS = [ iV, jV, kV ] + iaS

          dV_Left ( iV, jV, kV )  &
            =  V ( iV, jV, kV )  -  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
!    dV_Right = cshift ( V, shift = 1, dim = iD ) - V

    iaS = 0
    iaS ( iD ) = +1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
          
          !call Show ( [ iV, jV, kV ], '<<< iV, jV, kV: 2' )

          iaVS = [ iV, jV, kV ] + iaS

          dV_Right ( iV, jV, kV )  &
            =  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) - V ( iV, jV, kV )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( dV_Right )
    call DisassociateHost ( dV_Left )
    call DisassociateHost ( V )
    
  end subroutine ComputeDifferencesKernel


  subroutine ComputeReconstructionKernel &
               ( V, dV_Left, dV_Right, Theta, D_V, D_dV_Left, D_dV_Right, &
                 D_V_Inner, D_V_Outer, V_Inner, V_Outer )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V, &
      dV_Left, dV_Right
    real ( KDR ), intent ( in ) :: &
      Theta
    type ( c_ptr ), intent ( in ) :: &
      D_V, &
      D_dV_Left, D_dV_Right, &
      D_V_Inner, D_V_Outer
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      V_Inner, V_Outer

    !real ( KDR ), dimension ( size ( V ) ) :: &
    !  dV
    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      dV
      
    call AssociateHost ( D_V, V )
    call AssociateHost ( D_dV_Left, dV_Left )
    call AssociateHost ( D_dV_Right, dV_Right )
    call AssociateHost ( D_V_Inner, V_Inner )
    call AssociateHost ( D_V_Outer, V_Outer )

    !where ( dV_Left > 0.0_KDR .and. dV_Right > 0.0_KDR )
    !  dV = min ( Theta * dV_Left, Theta * dV_Right, &
    !               0.5_KDR * ( dV_Left + dV_Right ) )
    !elsewhere ( dV_Left < 0.0_KDR .and. dV_Right < 0.0_KDR )
    !  dV = max ( Theta * dV_Left, Theta * dV_Right, &
    !               0.5_KDR * ( dV_Left + dV_Right ) )      
    !elsewhere
    !  dV = 0.0_KDR
    !endwhere

    !V_Inner = V - 0.5_KDR * dV
    !V_Outer = V + 0.5_KDR * dV
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( dV )
    do iV = 1, size ( V )
      dV = ( sign ( 0.5_KDR, dV_Left ( iV ) ) &
             + sign ( 0.5_KDR, dV_Right ( iV ) ) ) &
             * min ( abs ( Theta * dV_Left ( iV ) ), &
                     abs ( Theta * dV_Right ( iV ) ), &
                     abs ( 0.5_KDR * ( dV_Left ( iV ) + dV_Right ( iV ) ) ) )
      V_Inner ( iV ) = V ( iV ) - 0.5_KDR * dV
      V_Outer ( iV ) = V ( iV ) + 0.5_KDR * dV
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( V_Outer )
    call DisassociateHost ( V_Inner )
    call DisassociateHost ( dV_Right )
    call DisassociateHost ( dV_Left )
    call DisassociateHost ( V )
    
  end subroutine ComputeReconstructionKernel


  subroutine ComputeFluxesKernel &
               ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, U_I, U_O, oV, iD, &
                 D_AP_I, D_AP_O, D_AM_I, D_AM_O, D_RF_I, D_RF_O, &
                 D_U_I, D_U_O, D_F_I, D_F_O, F_I, F_O )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      AP_I, AP_O, &
      AM_I, AM_O, &
      RF_I, RF_O, &
      U_I, U_O
    integer ( KDI ), intent ( in ) :: &
      oV, &
      iD
    type ( c_ptr ), intent ( in ) :: &
      D_AP_I, D_AP_O, &
      D_AM_I, D_AM_O, &
      D_RF_I, D_RF_O, &
      D_U_I, D_U_O, &
      D_F_I, D_F_O
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      F_I, F_O
      
    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV
      
    call AssociateHost ( D_AP_I, AP_I )
    call AssociateHost ( D_AP_O, AP_O )
    call AssociateHost ( D_AM_I, AM_I )
    call AssociateHost ( D_AM_O, AM_O )
    call AssociateHost ( D_RF_I, RF_I )
    call AssociateHost ( D_RF_O, RF_O )
    call AssociateHost ( D_U_I, U_I )
    call AssociateHost ( D_U_O, U_O )
    call AssociateHost ( D_F_I, F_I )
    call AssociateHost ( D_F_O, F_O )
            
    lV = 1
    where ( shape ( F_I ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV

    uV = 1
    where ( shape ( F_I ) > 1 )
      uV = shape ( F_I ) - oV
    end where
    uV ( iD ) = size ( F_I, dim = iD ) - 1
    
    !where ( AP_I + AM_I > 0.0_KDR )
    !  F_I = ( AP_I * cshift ( RF_O, shift = -1, dim = iD )  +  AM_I * RF_I &
    !          - AP_I * AM_I * ( U_I - cshift ( U_O, shift = -1, dim = iD ) ) ) &
    !        / ( AP_I + AM_I )
    !elsewhere
    !  F_I = 0.0_KDR
    !end where
    
    iaS = 0
    iaS ( iD ) = -1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
            
          iaVS = [ iV, jV, kV ] + iaS
          
          if ( AP_I ( iV, jV, kV ) + AM_I ( iV, jV, kV ) > 0.0_KDR ) then
            F_I ( iV, jV, kV ) &
              = ( AP_I ( iV, jV, kV ) &
                    * RF_O ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                  + AM_I ( iV, jV, kV ) * RF_I ( iV, jV, kV ) &
                  - AP_I ( iV, jV, kV ) * AM_I ( iV, jV, kV ) &
                    * ( U_I ( iV, jV, kV ) &
                        -  U_O ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) ) &
                / ( AP_I ( iV, jV, kV ) + AM_I ( iV, jV, kV ) )
          else
            F_I ( iV, jV, kV ) = 0.0_KDR
          end if

        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    !where ( AP_O + AM_O > 0.0_KDR )
    !  F_O = ( AP_O * RF_O  +  AM_O * cshift ( RF_I, shift = +1, dim = iD ) &
    !          - AP_O * AM_O * ( cshift ( U_I, shift = +1, dim = iD ) - U_O ) ) &
    !        / ( AP_O + AM_O )
    !elsewhere
    !  F_O = 0.0_KDR
    !end where
    
    iaS = 0
    iaS ( iD ) = +1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
            
          iaVS = [ iV, jV, kV ] + iaS
          
          if ( AP_O ( iV, jV, kV ) + AM_O ( iV, jV, kV ) > 0.0_KDR ) then
            F_O ( iV, jV, kV ) &
              = ( AP_O ( iV, jV, kV ) * RF_O ( iV, jV, kV ) &
                  +  AM_O ( iV, jV, kV ) &
                     * RF_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                  - AP_O ( iV, jV, kV ) * AM_O ( iV, jV, kV ) &
                    * ( U_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                        - U_O ( iV, jV, kV ) ) ) &
                / ( AP_O ( iV, jV, kV ) + AM_O ( iV, jV, kV ) )
          else
            F_O ( iV, jV, kV ) = 0.0_KDR
          end if
          
        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( F_O )
    call DisassociateHost ( F_I )
    call DisassociateHost ( U_O )
    call DisassociateHost ( U_I )
    call DisassociateHost ( RF_O )
    call DisassociateHost ( RF_I )
    call DisassociateHost ( AM_O )
    call DisassociateHost ( AM_I )
    call DisassociateHost ( AP_O )
    call DisassociateHost ( AP_I )
    
  end subroutine ComputeFluxesKernel


  subroutine ComputeUpdateKernel &
               ( dU, F_I, F_O, V, A, dT, D_dU, D_F_I, D_F_O )
    
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dU
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      F_I, F_O
    real ( KDR ), intent ( in ) :: &
      V, &
      A, &
      dT
    type ( c_ptr ), intent ( in ) :: &
      D_dU, &
      D_F_I, D_F_O
    
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_dU, dU )
    call AssociateHost ( D_F_I, F_I )
    call AssociateHost ( D_F_O, F_O )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( dU )
      dU ( iV ) = dU ( iV ) - dT * ( F_O ( iV ) - F_I ( iV ) ) * A / V
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( F_O )
    call DisassociateHost ( F_I )
    call DisassociateHost ( dU )

  end subroutine ComputeUpdateKernel
  
  
  subroutine AddUpdateKernel ( O, U, D_O, D_U, D_C, C )
    
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      O, &
      U
    type ( c_ptr ), intent ( in ) :: &
      D_O, &
      D_U, &
      D_C
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      C 
      
    integer ( KDI ) :: &
      iV, &
      nV
    
    call AssociateHost ( D_O, O )
    call AssociateHost ( D_U, U )
    call AssociateHost ( D_C, C )
    
    nV = size ( O )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, nV
      C ( iV ) = O ( iV ) + U ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( C )
    call DisassociateHost ( U )
    call DisassociateHost ( O )

  end subroutine AddUpdateKernel
  
  
  subroutine CombineUpdatesKernel ( C, O, U, D_C, D_O, D_U )
    
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      C 
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      O, &
      U
    type ( c_ptr ), intent ( in ) :: &
      D_C, &
      D_O, &
      D_U
      
    integer ( KDI ) :: &
      iV, &
      nV
    
    call AssociateHost ( D_C, C )    
    call AssociateHost ( D_O, O )
    call AssociateHost ( D_U, U )
    
    nV = size ( O )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, nV
      C ( iV ) = 0.5_KDR * ( O ( iV ) + ( C ( iV ) + U ( iV ) ) )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( U )
    call DisassociateHost ( O )
    call DisassociateHost ( C )

  end subroutine CombineUpdatesKernel
  

end module ConservationLawStep_Form
