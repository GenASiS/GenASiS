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
      iTimerDifference, &
      iTimerReconstruction, &
      iTimerRiemannSolverInput, &
      iTimerDataTransferDevice, &
      iTimerDataTransferHost
    real ( KDR ) :: &
      LimiterParameter
    type ( StorageForm ) :: &
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
      ComputeFluxesKernel

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

    !-- Update

    call CLS % Update % Initialize &
           ( [ nCells, CF % N_CONSERVED ], NameOption = 'Update', &
             VariableOption &
               = [ ( CF % Variable ( CF % iaConserved ( iV ) ), &
                     iV = 1, CF % N_CONSERVED ) ] )

    end associate !-- nCells
    end associate !-- DM
    
    call PROGRAM_HEADER % AddTimer &
           ( 'Difference', CLS % iTimerDifference, Level = 1 )
    call PROGRAM_HEADER % AddTimer &
           ( 'Reconstruction', CLS % iTimerReconstruction, Level = 1 )
    call PROGRAM_HEADER % AddTimer &
           ( 'RiemannSolverInput', CLS % iTimerRiemannSolverInput, Level = 1 )
    call PROGRAM_HEADER % AddTimer &
           ( 'DataTransfer to Device', CLS % iTimerDataTransferDevice, Level = 1 )
    call PROGRAM_HEADER % AddTimer &
           ( 'DataTransfer to Host', CLS % iTimerDataTransferHost, Level = 1 )
    
  end subroutine Initialize


  subroutine Solve ( CLS, TimeStep )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iV  !-- iVariable
    type ( StorageForm ) :: &
      Old, &
      Primitive

    associate &
      ( CF => CLS % ConservedFields )
    associate &
      ( DM => CF % DistributedMesh, &
        Current => CF, &
        Update  => CLS % Update, &
        iaC => CF % iaConserved )

    call Old % Initialize ( [ Current % nValues, Current % N_CONSERVED ] )
    call Primitive % Initialize &
           ( Current, iaSelectedOption = Current % iaPrimitive )

    do iV = 1, Current % N_CONSERVED
      call Copy ( Current % Value ( :, iaC ( iV ) ), Old % Value ( :, iV ) )
    end do

    !-- Substep 1
    !call Show ( '<<< Substep 1')

    call CLS % ComputeUpdate ( TimeStep ) !-- K1 = dT * RHS
    
    !call Show ( '<<< After Compute Update')
    
    do iV = 1, Current % N_CONSERVED
      !-- Current = Old + K1
      Current % Value ( :, iaC ( iV ) ) &
        = Old % Value ( :, iV ) + Update % Value ( :, iV )
    end do

    call Current % ComputePrimitive ( Current % Value )
    call DM % StartGhostExchange ( Primitive )
    call Current % ComputeAuxiliary ( Current % Value )
    call DM % FinishGhostExchange ( )

    !-- Substep 2

    call CLS % ComputeUpdate ( TimeStep ) !-- K2 = dT * RHS
    
    do iV = 1, Current % N_CONSERVED
      !-- Current = Old + 0.5 * ( k1 + k2 )
      !           = 0.5 Old + 0.5 ( Old + k1 + k2 )            
      Current % Value ( :, iaC ( iV ) ) &
        = Current % Value ( :, iaC ( iV ) ) + Update % Value ( :, iV )
      Current % Value ( :, iaC ( iV ) ) &
        = 0.5_KDR * ( Old % Value ( :, iV ) &
                      + Current % Value ( :, iaC ( iV ) ) )
    end do

    call Current % ComputePrimitive ( Current % Value )
    call DM % StartGhostExchange ( Primitive )
    call Current % ComputeAuxiliary ( Current % Value )
    call DM % FinishGhostExchange ( )
    
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

    associate ( CF => CLS % ConservedFields )
    associate ( DM => CF % DistributedMesh )

    call Clear ( CLS % Update % Value )
    
    !call Show ( '<<< ComputeUpdate: Dimension loop' )
    do iD = 1, DM % nDimensions

      call CLS % ComputeDifferences ( iD )
      call CLS % ComputeReconstruction ( )
      call CLS % ComputeFluxes ( iD )

      do iV = 1, CF % N_CONSERVED
        call ComputeUpdateKernel &
               ( CLS % Update % Value ( :, iV ), &
                 CLS % FluxInner % Value ( :, iV ), &
                 CLS % FluxOuter % Value ( :, iV ), DM % CellVolume, &
                 DM % CellArea ( iD ), TimeStep )
      end do

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

    associate ( CF => CLS % ConservedFields )
    associate ( DM  => CF % DistributedMesh )
    associate &
      ( T_DT_D  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferDevice ), &
        T_DT_H  => PROGRAM_HEADER % Timer ( CLS % iTimerDataTransferHost ), &
        T_D     => PROGRAM_HEADER % Timer ( CLS % iTimerDifference ) )
    
    !call Show ('<<< Compute Differences')

    call CF % ApplyBoundaryConditions &
           ( CF % Value, CF % Value, iD, iBoundary = -1 )
    call CF % ApplyBoundaryConditions &
           ( CF % Value, CF % Value, iD, iBoundary = +1 )
           
    !call Show ('<<< After ApplyBoundary')
    
    call T_DT_D % Start ( )
    call CF % UpdateDevice ( CF % iaPrimitive ( 1 ) )
    call T_DT_D % Stop ( )
    
    do iP = 1, CF % N_PRIMITIVE
      if ( iP < CF % N_PRIMITIVE ) then
        call T_DT_D % Start ( )
        call CF % UpdateDevice ( CF % iaPrimitive ( iP + 1 ) )
        call T_DT_D % Stop ( )
      end if
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
      
      !call T_DT_H % Start ( )
      !call CLS % DifferenceLeft % UpdateHost ( iP )
      !call CLS % DifferenceRight % UpdateHost ( iP )
      !call T_DT_H % Stop ( )
    end do
    
    !call Show ( '<<< After Differences' )
    
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
        T_R     => PROGRAM_HEADER % Timer ( CLS % iTimerReconstruction ) )
    
    
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
      
      call T_DT_H % Start ( )
      call CLS % ReconstructionInner % UpdateHost ( iaP ( iP ) )
      call CLS % ReconstructionOuter % UpdateHost ( iaP ( iP ) )
      call T_DT_H % Stop ( )
    end do

    call CF % ComputeAuxiliary ( CLS % ReconstructionInner % Value )
    call CF % ComputeAuxiliary ( CLS % ReconstructionOuter % Value )
    call CF % ComputeConserved ( CLS % ReconstructionInner % Value )
    call CF % ComputeConserved ( CLS % ReconstructionOuter % Value )
    
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
        T_RSI => PROGRAM_HEADER % Timer ( CLS % iTimerRiemannSolverInput ) )

    call CF % ApplyBoundaryConditions &
           ( CLS % ReconstructionOuter % Value, &
             CLS % ReconstructionInner % Value, iDimension, iBoundary = -1 )
    call CF % ApplyBoundaryConditions &
           ( CLS % ReconstructionInner % Value, &
             CLS % ReconstructionOuter % Value, iDimension, iBoundary = +1 )

    call T_DT_D % Start ( )
    call CLS % ReconstructionInner % UpdateDevice ( )
    call CLS % ReconstructionOuter % UpdateDevice ( )
    call T_DT_D % Stop ( )
    
    call T_RSI % Start ( )
    call CF % ComputeRiemannSolverInput &
           ( CLS, CLS % ReconstructionInner % Value, &
             CLS % ReconstructionOuter % Value, iDimension, &
             CLS % ReconstructionInner % D_Selected, &
             CLS % ReconstructionOuter % D_Selected )
    call T_RSI % Stop ( )

    call CF % ComputeRawFluxes &
           ( CLS % RawFluxInner % Value, CLS % ReconstructionInner % Value, &
             iDimension )
    call CF % ComputeRawFluxes &
           ( CLS % RawFluxOuter % Value, CLS % ReconstructionOuter % Value, &
             iDimension )

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
      call ComputeFluxesKernel &
             ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, U_I, U_O, iDimension, &
               F_I, F_O )
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
      
    call AssociateDevice ( D_V, V )
    call AssociateDevice ( D_dV_Left, dV_Left )
    call AssociateDevice ( D_dV_Right, dV_Right )

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
    
    !call Show ( '<<< compute differences kernel' )

!    dV_Left  = V - cshift ( V, shift = -1, dim = iD )    

    iaS = 0
    iaS ( iD ) = -1
    
    !-- $OMP parallel do private ( iV, jV, kV, iaVS )
    
    !$OMP  target teams distribute parallel do collapse ( 3 ) &
    !$OMP& schedule ( static, 1 ) private ( iaVS )
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
    !$OMP end target teams distribute parallel do
    
    !-- $OMP end parallel do
    

!    dV_Right = cshift ( V, shift = 1, dim = iD ) - V

    iaS = 0
    iaS ( iD ) = +1
    
    !-- $OMP parallel do private ( iV, jV, kV, iaVS )
    
    !$OMP  target teams distribute parallel do collapse ( 3 ) &
    !$OMP& schedule ( static, 1 ) private ( iaVS )
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
    !$OMP end target teams distribute parallel do
    
    !-- $OMP end parallel do

    call DisassociateDevice ( dV_Right )
    call DisassociateDevice ( dV_Left )
    call DisassociateDevice ( V )
    
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
      
    call AssociateDevice ( D_V, V )
    call AssociateDevice ( D_dV_Left, dV_Left )
    call AssociateDevice ( D_dV_Right, dV_Right )
    call AssociateDevice ( D_V_Inner, V_Inner )
    call AssociateDevice ( D_V_Outer, V_Outer )

      
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
    
    !$OMP  target teams distribute parallel do &
    !$OMP& schedule ( static, 1 ) private ( dV )
    do iV = 1, size ( V )
      dV = ( sign ( 0.5_KDR, dV_Left ( iV ) ) &
             + sign ( 0.5_KDR, dV_Right ( iV ) ) ) &
             * min ( abs ( Theta * dV_Left ( iV ) ), &
                     abs ( Theta * dV_Right ( iV ) ), &
                     abs ( 0.5_KDR * ( dV_Left ( iV ) + dV_Right ( iV ) ) ) )
      V_Inner ( iV ) = V ( iV ) - 0.5_KDR * dV
      V_Outer ( iV ) = V ( iV ) + 0.5_KDR * dV
    end do
    !$OMP end target teams distribute parallel do
    
    call DisassociateDevice ( V_Outer )
    call DisassociateDevice ( V_Inner )
    call DisassociateDevice ( dV_Right )
    call DisassociateDevice ( dV_Left )
    call DisassociateDevice ( V )
    
  end subroutine ComputeReconstructionKernel


  subroutine ComputeFluxesKernel &
               ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, U_I, U_O, iD, F_I, F_O )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      AP_I, AP_O, &
      AM_I, AM_O, &
      RF_I, RF_O, &
      U_I, U_O
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      F_I, F_O

    where ( AP_I + AM_I > 0.0_KDR )
      F_I = ( AP_I * cshift ( RF_O, shift = -1, dim = iD )  +  AM_I * RF_I &
              - AP_I * AM_I * ( U_I - cshift ( U_O, shift = -1, dim = iD ) ) ) &
            / ( AP_I + AM_I )
    elsewhere
      F_I = 0.0_KDR
    end where

    where ( AP_O + AM_O > 0.0_KDR )
      F_O = ( AP_O * RF_O  +  AM_O * cshift ( RF_I, shift = +1, dim = iD ) &
              - AP_O * AM_O * ( cshift ( U_I, shift = +1, dim = iD ) - U_O ) ) &
            / ( AP_O + AM_O )
    elsewhere
      F_O = 0.0_KDR
    end where

  end subroutine ComputeFluxesKernel


  subroutine ComputeUpdateKernel ( dU, F_I, F_O, V, A, dT )
    
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dU
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      F_I, F_O
    real ( KDR ), intent ( in ) :: &
      V, &
      A, &
      dT

    dU = dU - dT * ( F_O - F_I ) * A / V
 
  end subroutine ComputeUpdateKernel


end module ConservationLawStep_Form
