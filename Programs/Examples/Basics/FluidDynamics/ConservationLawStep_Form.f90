module ConservationLawStep_Form

  use Basics
  use ConservedFields_Template

  implicit none
  private

  type, public :: ConservationLawStepForm
    integer ( KDI ) :: &
      ALPHA_PLUS  = 1, &
      ALPHA_MINUS = 2, &
      N_MODIFIED_SPEEDS = 2
    real ( KDR ) :: &
      LimiterParameter
    type ( VariableGroupForm ) :: &
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

  end subroutine Initialize


  subroutine Solve ( CLS, TimeStep )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iV  !-- iVariable
    type ( VariableGroupForm ) :: &
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

    call CLS % ComputeUpdate ( TimeStep ) !-- K1 = dT * RHS
    
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


  subroutine ComputeDifferences ( CLS, iDimension )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iV  !-- iVariable
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V, &
      dV_Left, &
      dV_Right

    associate ( CF => CLS % ConservedFields )
    associate ( DM  => CF % DistributedMesh )

    call CF % ApplyBoundaryConditions &
           ( CF % Value, CF % Value, iDimension, iBoundary = -1 )
    call CF % ApplyBoundaryConditions &
           ( CF % Value, CF % Value, iDimension, iBoundary = +1 )

    do iV = 1, CF % N_PRIMITIVE
      call DM % SetVariablePointer &
             ( CF % Value ( :, CF % iaPrimitive ( iV ) ), V )
      call DM % SetVariablePointer &
             ( CLS % DifferenceLeft % Value ( :, iV ), dV_Left )
      call DM % SetVariablePointer &
             ( CLS % DifferenceRight % Value ( :, iV ), dV_Right )
      call ComputeDifferencesKernel ( V, iDimension, dV_Left, dV_Right )
    end do

    nullify ( V, dV_Left, dV_Right )
    end associate !-- DM
    end associate !-- CF

  end subroutine ComputeDifferences


  subroutine ComputeReconstruction ( CLS )

    class ( ConservationLawStepForm ), intent ( inout ) :: &
      CLS

    integer ( KDI ) :: &
      iV  !-- iVariable

    associate ( CF => CLS % ConservedFields )
    associate ( iaP => CF % iaPrimitive )

    do iV = 1, CF % N_PRIMITIVE
      call ComputeReconstructionKernel &
             ( CF % Value ( :, iaP ( iV ) ), &
               CLS % DifferenceLeft % Value ( :, iV ), &
               CLS % DifferenceRight % Value ( :, iV ), &
               CLS % LimiterParameter, &
               CLS % ReconstructionInner % Value ( :, iaP ( iV ) ), &
               CLS % ReconstructionOuter % Value ( :, iaP ( iV ) ) )
    end do

    call CF % ComputeAuxiliary ( CLS % ReconstructionInner % Value )
    call CF % ComputeAuxiliary ( CLS % ReconstructionOuter % Value )
    call CF % ComputeConserved ( CLS % ReconstructionInner % Value )
    call CF % ComputeConserved ( CLS % ReconstructionOuter % Value )

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
        iaC => CF % iaConserved )

    call CF % ApplyBoundaryConditions &
           ( CLS % ReconstructionOuter % Value, &
             CLS % ReconstructionInner % Value, iDimension, iBoundary = -1 )
    call CF % ApplyBoundaryConditions &
           ( CLS % ReconstructionInner % Value, &
             CLS % ReconstructionOuter % Value, iDimension, iBoundary = +1 )

    call CF % ComputeRiemannSolverInput &
           ( CLS, CLS % ReconstructionInner % Value, &
             CLS % ReconstructionOuter % Value, iDimension )

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


  subroutine ComputeDifferencesKernel ( V, iDimension, dV_Left, dV_Right )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      V
    integer ( KDI ), intent ( in ) :: &
      iDimension
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dV_Left, &
      dV_Right

    dV_Left  = V - cshift ( V, shift = -1, dim = iDimension )
    dV_Right = cshift ( V, shift = 1, dim = iDimension ) - V

  end subroutine ComputeDifferencesKernel


  subroutine ComputeReconstructionKernel &
               ( V, dV_Left, dV_Right, Theta, V_Inner, V_Outer )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V, &
      dV_Left, &
      dV_Right
    real ( KDR ), intent ( in ) :: &
      Theta
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      V_Inner, &
      V_Outer

    real ( KDR ), dimension ( size ( V ) ) :: &
      dV

    where ( dV_Left > 0.0_KDR .and. dV_Right > 0.0_KDR )
      dV = min ( Theta * dV_Left, Theta * dV_Right, &
                   0.5_KDR * ( dV_Left + dV_Right ) )
    elsewhere ( dV_Left < 0.0_KDR .and. dV_Right < 0.0_KDR )
      dV = max ( Theta * dV_Left, Theta * dV_Right, &
                   0.5_KDR * ( dV_Left + dV_Right ) )      
    elsewhere
      dV = 0.0_KDR
    endwhere

    V_Inner = V - 0.5_KDR * dV
    V_Outer = V + 0.5_KDR * dV

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
