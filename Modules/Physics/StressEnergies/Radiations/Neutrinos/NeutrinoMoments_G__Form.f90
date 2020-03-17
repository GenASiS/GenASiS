module NeutrinoMoments_G__Form

  !-- NeutrinoMoments_Grey__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics
  use Photons

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_NM = 1, &
      N_CONSERVED_NM = 1, &
      N_FIELDS_NM    = 6, &
      N_VECTORS_NM   = 0

  type, public, extends ( PhotonMoments_G_Form ) :: NeutrinoMoments_G_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_NM          = N_PRIMITIVE_NM, &
      N_CONSERVED_NM          = N_CONSERVED_NM, &
      N_FIELDS_NM             = N_FIELDS_NM, &
      N_VECTORS_NM            = N_VECTORS_NM, &
      COMOVING_NUMBER         = 0, &
      CONSERVED_NUMBER        = 0, &
      DEGENERACY_PARAMETER    = 0, &
      DEGENERACY_PARAMETER_EQ = 0, &
      ENERGY_AVERAGE          = 0, &
      OCCUPANCY_AVERAGE       = 0!, &
  !     BETA_EQUILIBRIUM        = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetOutput
    final :: &
      Finalize
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, pass ( C ) :: &
      ComputeDiffusionFactor_HLL
    procedure, private, pass ( NM ) :: &
      ComputeSpectralParameters_NM
    generic, public :: &
      ComputeSpectralParameters => ComputeSpectralParameters_NM
  end type NeutrinoMoments_G_Form

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeConservedNumber, &
      ComputeComovingNumber, &
      ComputeRawFluxesKernel, &
      EvaluateZeroEta

contains


  subroutine InitializeAllocate_RM &
               ( RM, RadiationType, MomentsType, RiemannSolverType, &
                 ReconstructedType, UseLimiter, Units, LimiterParameter, &
                 nValues, VariableOption, VectorOption, NameOption, &
                 ClearOption, UnitOption, VectorIndicesOption )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      RM
    character ( * ), intent ( in ) :: &
      RadiationType, &
      MomentsType, &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
    real ( KDR ), intent ( in ) :: &
      LimiterParameter
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit

    call InitializeBasics &
           ( RM, Variable, VariableUnit, VariableOption, UnitOption )

    call SetUnits ( VariableUnit, RM, Units )

    call RM % PhotonMoments_G_Form % Initialize &
           ( RadiationType, MomentsType, RiemannSolverType, ReconstructedType, &
             UseLimiter, Units, LimiterParameter, nValues, &
             VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeAllocate_RM


  subroutine SetPrimitiveConserved ( C )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_NM ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_NM ) :: &
      ConservedName

    oP = C % N_PRIMITIVE_TEMPLATE + C % N_PRIMITIVE_PM + C % N_PRIMITIVE_RM
    oC = C % N_CONSERVED_TEMPLATE + C % N_CONSERVED_PM + C % N_CONSERVED_RM

    call Show( oP, 'oP' )
    call Show( oC, 'oC' )

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_NM
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    call Show( C % COMOVING_NUMBER, 'COMOVING_NUMBER' )
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_NM ) &
      = [ C % COMOVING_NUMBER ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_NM
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_NM ) &
      = [ C % CONSERVED_NUMBER ]

    call Show( C % iaConserved, 'iaConserved' )
    call Show( C % iaPrimitive, 'iaPrimitive' )

    do iF = 1, C % N_PRIMITIVE_NM
      call Show( oP+iF, 'oP+iF' )
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_NM
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )
    
    call C % RadiationMomentsForm % SetPrimitiveConserved ( )

  end subroutine SetPrimitiveConserved


  subroutine SetOutput ( RM, Output )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      RM
    type ( StorageForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( RM % COMOVING_MOMENTUM_U )

    call Output % Initialize &
           ( RM, iaSelectedOption = [ RM % COMOVING_ENERGY, &
                                      RM % COMOVING_MOMENTUM_U, &
                                      RM % FLUID_VELOCITY_U, &
                                      RM % FLUX_FACTOR, &
                                      RM % STRESS_FACTOR, &
                                      RM % TEMPERATURE_PARAMETER, &
                                      RM % TEMPERATURE_PARAMETER_EQ, &
                                      RM % COMOVING_NUMBER, &
                                      RM % DEGENERACY_PARAMETER, &
                                      RM % DEGENERACY_PARAMETER_EQ ], &
!                                      RM % ENERGY_AVERAGE, &
!                                      RM % OCCUPANCY_AVERAGE, &
!                                      RM % BETA_EQUILIBRIUM, &
!                                      RM % DIFFUSION_FACTOR_E ], &
             VectorOption = [ 'ComovingMomentum_U' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  impure elemental subroutine Finalize ( NM )

    type ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      NM

    !-- Trigger finalization in parent

  end subroutine Finalize


  subroutine ComputeFromPrimitiveCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
    real ( KDR ), dimension ( :, : ), pointer :: &
      RMV
      
    call C % RadiationMomentsForm % ComputeFromPrimitiveCommon &
           ( Value_C, G, Value_G, nValuesOption, oValueOption )

    RMV => Value_C
    associate ( GV => Value_G )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( RMV, dim = 1 )
    end if
    
    associate &
      ( M_DD_22 => GV ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => GV ( oV + 1 : oV + nV, G % METRIC_DD_33 ) )
    associate &
      ( J       => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY ), &
        H_1     => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 1 ) ), &
        H_2     => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 2 ) ), &
        H_3     => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 3 ) ), &
        T       => RMV ( oV + 1 : oV + nV, C % TEMPERATURE_PARAMETER ), &
        T_EQ    => RMV ( oV + 1 : oV + nV, C % TEMPERATURE_PARAMETER_EQ ), &
        N       => RMV ( oV + 1 : oV + nV, C % COMOVING_NUMBER ), &
        D       => RMV ( oV + 1 : oV + nV, C % CONSERVED_NUMBER ), &
        Eta     => RMV ( oV + 1 : oV + nV, C % DEGENERACY_PARAMETER ), &
        Eta_EQ  => RMV ( oV + 1 : oV + nV, C % DEGENERACY_PARAMETER_EQ ), &
        E_Ave   => RMV ( oV + 1 : oV + nV, C % ENERGY_AVERAGE ), &
        F_Ave   => RMV ( oV + 1 : oV + nV, C % OCCUPANCY_AVERAGE ), &
        V_1   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 1 ) ), &
        V_2   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 2 ) ), &
        V_3   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 2 ) ) )

    call ComputeConservedNumber &
           ( D, N, H_1, H_2, H_3, E_Ave, M_DD_22, M_DD_33, V_1, V_2, V_3 )

    if ( associated ( C % Interactions ) ) then
      call C % Interactions % ComputeEquilibriumParameters ( T_EQ, Eta_EQ )
      associate ( I => C % Interactions )
      associate &
        ( J_EQ => I % Value ( :, I % EQUILIBRIUM_J ), &
          N_EQ => I % Value ( :, I % EQUILIBRIUM_N ) )
      call C % ComputeSpectralParameters &
             ( T, Eta, J, N, E_Ave_Option = E_Ave, F_Ave_Option = F_Ave, &
               J_EQ_Option = J_EQ, N_EQ_Option = N_EQ, T_EQ_Option = T_EQ, &
               Eta_EQ_Option = Eta_EQ )
      end associate !-- J_EQ, etc.
      end associate !-- I
    else
      call C % ComputeSpectralParameters ( T, Eta, J, N )
    end if  !-- associated C % Interactions

    end associate !-- J, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- GV
    nullify ( RMV )

  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
    real ( KDR ), dimension ( :, : ), pointer :: &
      RMV
      
    call C % RadiationMomentsForm % ComputeFromConservedCommon &
           ( Value_C, G, Value_G, nValuesOption, oValueOption )

    RMV => Value_C
    associate ( GV => Value_G )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( RMV, dim = 1 )
    end if
    
    associate &
      ( M_DD_22 => GV ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => GV ( oV + 1 : oV + nV, G % METRIC_DD_33 ) )
    associate &
      ( J       => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY ), &
        H_1     => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 1 ) ), &
        H_2     => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 2 ) ), &
        H_3     => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 3 ) ), &
        T       => RMV ( oV + 1 : oV + nV, C % TEMPERATURE_PARAMETER ), &
        T_EQ    => RMV ( oV + 1 : oV + nV, C % TEMPERATURE_PARAMETER_EQ ), &
        N       => RMV ( oV + 1 : oV + nV, C % COMOVING_NUMBER ), &
        D       => RMV ( oV + 1 : oV + nV, C % CONSERVED_NUMBER ), &
        Eta     => RMV ( oV + 1 : oV + nV, C % DEGENERACY_PARAMETER ), &
        Eta_EQ  => RMV ( oV + 1 : oV + nV, C % DEGENERACY_PARAMETER_EQ ), &
        E_Ave   => RMV ( oV + 1 : oV + nV, C % ENERGY_AVERAGE ), &
        F_Ave   => RMV ( oV + 1 : oV + nV, C % OCCUPANCY_AVERAGE ), & 
        V_1   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 1 ) ), &
        V_2   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 2 ) ), &
        V_3   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 2 ) ) )

    call ComputeComovingNumber &
           ( N, D, H_1, H_2, H_3, E_Ave, M_DD_22, M_DD_33, V_1, V_2, V_3 )

    if ( associated ( C % Interactions ) ) then
      call C % Interactions % ComputeEquilibriumParameters ( T_EQ, Eta_EQ )
      associate ( I => C % Interactions )
      associate &
        ( J_EQ => I % Value ( :, I % EQUILIBRIUM_J ), &
          N_EQ => I % Value ( :, I % EQUILIBRIUM_N ) )
      call C % ComputeSpectralParameters &
             ( T, Eta, J, N, E_Ave_Option = E_Ave, F_Ave_Option = F_Ave, &
               J_EQ_Option = J_EQ, N_EQ_Option = N_EQ, T_EQ_Option = T_EQ, &
               Eta_EQ_Option = Eta_EQ )
      end associate !-- J_EQ, etc.
      end associate !-- I
    else
      call C % ComputeSpectralParameters ( T, Eta, J, N )
    end if  !-- associated C % Interactions

    end associate !-- J, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- GV
    nullify ( RMV )

  end subroutine ComputeFromConservedCommon


  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Value_C, Value_G, iDimension, &
                 nValuesOption, oValueOption )
    
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_C, &
      Value_G
    integer ( KDI ), intent ( in ) :: &
      iDimension
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      iNumber
    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues

    call C % PhotonMoments_G_Form % ComputeRawFluxes &
           ( RawFlux, G, Value_C, Value_G, iDimension, &
             nValuesOption, oValueOption )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( Value_C, dim = 1 )
    end if
    
    call Search ( C % iaConserved, C % CONSERVED_NUMBER, iNumber )
    
    associate &
      ( F_D   => RawFlux ( oV + 1 : oV + nV, iNumber ), &
        J     => Value_C ( oV + 1 : oV + nV, C % COMOVING_ENERGY ), &
        H_Dim => Value_C ( oV + 1 : oV + nV, &
                           C % COMOVING_MOMENTUM_U ( iDimension ) ), &
        N     => Value_C ( oV + 1 : oV + nV, C % COMOVING_NUMBER ), &
        V_Dim => Value_C ( oV + 1 : oV + nV, &
                           C % FLUID_VELOCITY_U ( iDimension ) ) )

    call ComputeRawFluxesKernel ( F_D, J, H_Dim, N, V_Dim )

    end associate !-- F_E, etc.

  end subroutine ComputeRawFluxes
  

  subroutine ComputeDiffusionFactor_HLL ( DF_I, Grid, C, iDimension )

    type ( StorageForm ), intent ( inout ) :: &
      DF_I
    class ( * ), intent ( in ), target :: &
      Grid
    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iEnergy, &
      iNumber

    if ( .not. associated ( C % Interactions ) ) &
      return

    call C % PhotonMoments_G_Form % ComputeDiffusionFactor_HLL &
           ( DF_I, Grid, iDimension )

    call Search ( C % iaConserved, C % CONSERVED_ENERGY, iEnergy )
    call Search ( C % iaConserved, C % CONSERVED_NUMBER, iNumber )

    call Copy ( DF_I % Value ( :, iEnergy ), DF_I % Value ( :, iNumber ) )

  end subroutine ComputeDiffusionFactor_HLL


  subroutine ComputeSpectralParameters_NM &
               ( T, Eta, NM, J, N, E_Ave_Option, F_Ave_Option, &
                 J_EQ_Option, N_EQ_Option, T_EQ_Option, Eta_EQ_Option )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      Eta
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      NM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      N
    real ( KDR ), dimension ( : ), intent ( inout ), optional :: &
      E_Ave_Option, &
      F_Ave_Option, &
      J_EQ_Option, &
      N_EQ_Option
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      T_EQ_Option, &
      Eta_EQ_Option

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      OnePlusEpsilon, &
      MomentRatio, &
      Factor_ND, Factor_ED_1, Factor_ED_2, &
      Factor_J_N, &
      Eta_ND, Eta_ED, Eta_Sol, & !-- NonDegenerate, ExtremeDegenerate, Solution
      Fermi_2, Fermi_3, &
      Fermi_2_EQ, Fermi_3_EQ, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta
    type ( RootFinderForm ) :: &
      RF

    call RF % Initialize ( MomentRatio )
    RF % EvaluateZero  =>  EvaluateZeroEta

    nValues = size ( T )

    OnePlusEpsilon  =  1.0_KDR  +  10.0_KDR * epsilon ( 0.0_KDR )

    associate &
      ( Pi        => CONSTANT % PI, &
        TwoPi     => 2.0_KDR  *  CONSTANT % PI, &
        FourPi    => 4.0_KDR  *  CONSTANT % PI, &
        SixPi_2   => 6.0_KDR  *  CONSTANT % PI ** 2, &
        EightPi_2 => 8.0_KDR  *  CONSTANT % PI ** 2 )
    Factor_ND    =  Pi ** ( - 2.0_KDR / 3.0_KDR )  /  3.0_KDR
    Factor_ED_1  =  EightPi_2 ** ( 1.0_KDR / 4.0_KDR )  &
                    /  SixPi_2 ** ( 1.0_KDR / 3.0_KDR )
    Factor_ED_2  =  6.0_KDR  /  Pi ** 2
    Factor_J_N   =  FourPi  /  TwoPi ** 3
    end associate !-- TwoPi, etc.

!call Show ( minval ( pack ( J / N ** (4./3.), mask = N > 0.0_KDR ) ), '>>> min MomentRatio' )

    !$OMP parallel do &
    !$OMP   private ( iV, MomentRatio, Eta_ND, Eta_ED, Eta_Sol, &
    !$OMP             Fermi_2, Fermi_3, &
    !$OMP             fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
    do iV = 1, nValues

!call Show ( iV, '>>> iV' )
!call Show ( N ( iV ) ** ( 4./3. ), '>>> N^(4/3)' )
!call Show ( J ( iV ), '>>> J' )
      if ( N ( iV ) > 0.0_KDR .and. J ( iV ) > 0.0_KDR ) then
        MomentRatio  =  J ( iV ) ** ( 1.0_KDR / 4.0_KDR )  &
                        *  N ( iV ) ** ( - 1.0_KDR / 3.0_KDR )
        MomentRatio  =  max ( MomentRatio, OnePlusEpsilon / Factor_ED_1 )
        Eta_ND  =  - 3.0_KDR  * log ( Factor_ND  *  MomentRatio ** 4 )
        Eta_ED  =  ( Factor_ED_2 * ( Factor_ED_1 * MomentRatio - 1.0_KDR ) ) &
                   ** ( - 0.5_KDR )
!call Show ( MomentRatio, '>>> MomentRatio' )
!call Show ( Eta_ND, '>>> Eta_ND' )
        if ( Eta_ND < -10.0_KDR ) then
          Eta ( iV )  =  Eta_ND
!call Show ( '>>> Eta = Eta_ND' )
        else if ( Eta_ED > 50.0_KDR ) then
          Eta ( iV )  =  Eta_ED
        else
!call Show ( Eta ( iV ), '>>> Eta pre' )
!          call RF % Solve ( Eta ( iV ), 1.1 * Eta ( iV ), Eta ( iV ) )
          Eta_Sol = max ( Eta_ND, Eta ( iV ) )
          call RF % Solve ( 0.99 * Eta_Sol, Eta_Sol, Eta_Sol )
!call Show ( RF % Success, '>>> Success' )
!call Show ( RF % SolutionAccuracy, '>>> SolutionAccuracy' )
!call Show ( RF % nIterations, '>>> nIterations' )
          if ( RF % Success ) then
            Eta ( iV )  =  Eta_Sol
          else
call Show ( '>>> Eta fail', CONSOLE % ERROR )
call Show ( NM % Name, '>>> Species', CONSOLE % ERROR )
call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
call Show ( iV, '>>> iV', CONSOLE % ERROR )
call Show ( J ( iV ), '>>> J', CONSOLE % ERROR )
call Show ( N ( iV ), '>>> N', CONSOLE % ERROR )
call Show ( MomentRatio, '>>> MomentRatio', CONSOLE % ERROR )
call Show ( 1.0_KDR / Factor_ED_1, '>>> MomentRatioMin', CONSOLE % ERROR )
call Show ( Eta ( iV ), '>>> Eta previous', CONSOLE % ERROR )
call Show ( Eta_ND, '>>> Eta_ND', CONSOLE % ERROR )
call Show ( Eta_ED, '>>> Eta_ED', CONSOLE % ERROR )
call Show ( Eta_Sol, '>>> Eta_Sol', CONSOLE % ERROR )
call Show ( RF % SolutionAccuracy, '>>> SolutionAccuracy', CONSOLE % ERROR )
            if ( Eta_ND < 1.0_KDR ) then
call Show ( Eta_ND, '>>> Falling back to Eta_ND', CONSOLE % ERROR )
              Eta ( iV )  =  Eta_ND
            else
call Show ( Eta_ED, '>>> Falling back to Eta_ED', CONSOLE % ERROR )
              Eta ( iV )  =  Eta_ED
            end if
!call PROGRAM_HEADER % Abort ( )
          end if
        end if
!call Show ( Eta ( iV ), '>>> Eta post' )
      else
        Eta ( iV )  =  log ( tiny ( 0.0_KDR ) )  *  10.0_KDR
      end if  !-- N > 0 and J > 0

!      if ( trim ( NM % Type ) == 'NEUTRINOS_E_NU' ) then
!        Eta ( iV )  =  min ( Eta ( iV ), Eta_EQ ( iV ) )
!      end if

      call DFERMI ( 2.0_KDR, Eta ( iV ), 0.0_KDR, Fermi_2, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 3.0_KDR, Eta ( iV ), 0.0_KDR, Fermi_3, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      if ( N ( iV ) > 0.0_KDR .and. J ( iV ) > 0.0_KDR ) then
        T  ( iV )     =  J ( iV )  /  N ( iV )  *  Fermi_2 / Fermi_3
      else
        T ( iV )      =  0.0_KDR
      end if

    end do !-- iV
    !$OMP end parallel do

    if ( present ( E_Ave_Option ) .and. present ( F_Ave_Option ) &
         .and. present ( J_EQ_Option ) .and. present ( N_EQ_Option ) &
         .and. present ( T_EQ_Option ) .and. present ( Eta_EQ_Option ) ) &
    then

      !$OMP parallel do &
      !$OMP   private ( iV, Fermi_2_EQ, Fermi_3_EQ, &
      !$OMP             fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nValues

        call DFERMI ( 2.0_KDR, Eta_EQ_Option ( iV ), 0.0_KDR, Fermi_2_EQ, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 3.0_KDR, Eta_EQ_Option ( iV ), 0.0_KDR, Fermi_3_EQ, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

        if ( N ( iV ) > 0.0_KDR .and. J ( iV ) > 0.0_KDR ) then
          E_Ave_Option ( iV )  =  J ( iV )  /  N ( iV )
          F_Ave_Option ( iV )  &
            =  1.0_KDR &
               /  ( exp ( E_Ave_Option ( iV ) / T ( iV )  &
                          -  Eta ( iV ) )  +  1.0_KDR )
        else
          E_Ave_Option ( iV )  =  0.0_KDR
          F_Ave_Option ( iV )  =  0.0_KDR
        end if

        N_EQ_Option ( iV )  &
          =  Factor_J_N  *  T_EQ_Option ( iV ) ** 3  *  Fermi_2_EQ
        J_EQ_Option ( iV )  &
          =  Factor_J_N  *  T_EQ_Option ( iV ) ** 4  *  Fermi_3_EQ

      end do !-- iV
      !$OMP end parallel do

    end if !-- options present

  end subroutine ComputeSpectralParameters_NM


  subroutine InitializeBasics &
               ( NM, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      NM
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    type ( MeasuredValueForm ), dimension ( : ), optional, intent ( in ) :: &
      VariableUnitOption

    integer ( KDI ) :: &
      oF, &  !-- oField
      oV     !-- oVector

    if ( NM % Type == '' ) &
      NM % Type = 'a NeutrinoMoments_G'

    !-- variable indices

    oF = NM % N_FIELDS_TEMPLATE + NM % N_FIELDS_RM + NM % N_FIELDS_PM
    if ( NM % N_FIELDS == 0 ) &
      NM % N_FIELDS = oF + NM % N_FIELDS_NM 

    NM % COMOVING_NUMBER          =  oF + 1
    NM % CONSERVED_NUMBER         =  oF + 2
    NM % DEGENERACY_PARAMETER     =  oF + 3
    NM % DEGENERACY_PARAMETER_EQ  =  oF + 4
    NM % ENERGY_AVERAGE           =  oF + 5
    NM % OCCUPANCY_AVERAGE        =  oF + 6
    ! NM % BETA_EQUILIBRIUM         =  oF + 7

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( NM % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + NM % N_FIELDS_NM ) &
      = [ 'ComovingNumber        ', &
          'ConservedNumber       ', &
          'DegeneracyParameter   ', &
          'DegeneracyParameter_EQ', &
          'EnergyAverage         ', &
          'OccupancyAverage      ' ]!, &
!          'BetaEquilibrium       ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( NM % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = NM % N_VECTORS_TEMPLATE + NM % N_VECTORS_RM + NM % N_VECTORS_PM
    if ( NM % N_VECTORS == 0 ) &
      NM % N_VECTORS = oV + NM % N_VECTORS_NM

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, NM, Units )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      NM
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units

    VariableUnit ( NM % COMOVING_NUMBER ) &
      =  Units % EnergyDensity / Units % Temperature
    VariableUnit ( NM % CONSERVED_NUMBER ) &
      =  Units % EnergyDensity / Units % Temperature

   VariableUnit ( NM % ENERGY_AVERAGE )  =  Units % Temperature

  end subroutine SetUnits


  subroutine ComputeConservedNumber &
               ( D, N, H_1, H_2, H_3, E_Ave, M_DD_22, M_DD_33, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      D
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      H_1, H_2, H_3, &
      E_Ave, &
      M_DD_22, M_DD_33, &
      V_1, V_2, V_3

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( D )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV )  >  0.0_KDR ) then
        D ( iV )  =  N ( iV )
        if ( E_Ave ( iV ) > 0.0_KDR ) &
          D ( iV )  =  D ( iV )  &
                       +  (                       V_1 ( iV )  *  H_1 ( iV )  &
                            +  M_DD_22 ( iV )  *  V_2 ( iV )  *  H_2 ( iV )  &
                            +  M_DD_33 ( iV )  *  V_3 ( iV )  *  H_3 ( iV ) ) &
                          /  E_Ave ( iV )
      else
        D ( iV )  =  0.0_KDR
        N ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeConservedNumber


  subroutine ComputeComovingNumber &
               ( N, D, H_1, H_2, H_3, E_Ave, M_DD_22, M_DD_33, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      D
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      H_1, H_2, H_3, &
      E_Ave, &
      M_DD_22, M_DD_33, &
      V_1, V_2, V_3

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( N )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( D ( iV )  >  0.0_KDR ) then
        N ( iV )  =  D ( iV )
        if ( E_Ave ( iV ) > 0.0_KDR ) &
          N ( iV )  =  N ( iV )  &
                       -  (                       V_1 ( iV )  *  H_1 ( iV )  &
                            +  M_DD_22 ( iV )  *  V_2 ( iV )  *  H_2 ( iV )  &
                            +  M_DD_33 ( iV )  *  V_3 ( iV )  *  H_3 ( iV ) ) &
                          /  E_Ave ( iV )
      else
        N ( iV )  =  0.0_KDR
        D ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeComovingNumber


  subroutine ComputeRawFluxesKernel ( F_D, J, H_Dim, N, V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_D
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      H_Dim, &
      N, &
      V_Dim

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      ff_D

    nValues = size ( F_D )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      F_D ( iV )  =  N ( iV )  *  V_Dim ( iV )
      if ( J ( iV ) > 0.0_KDR ) &
        F_D ( iV )  =  F_D ( iV )  &
                       +  ( H_Dim ( iV )  /  J ( iV ) )  *  N ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeRawFluxesKernel


  subroutine EvaluateZeroEta ( Parameters, Input, Result )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      Input
    real ( KDR ), intent ( out ) :: &
      Result

    real ( KDR ) :: &
      Pi, &
      Factor, &
      X, &
      Fermi_2, Fermi_3, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta

    Pi      =  CONSTANT % PI
    Factor  =  ( 2 * Pi ** 2 ) ** ( 1. / 12. )

    X  =  min ( Input, log ( huge ( 1.0_KDR ) ) - 10.0_KDR )

    call DFERMI ( 2.0_KDR, X, 0.0_KDR, Fermi_2, &
                  fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    call DFERMI ( 3.0_KDR, X, 0.0_KDR, Fermi_3, &
                  fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

    select type ( P => Parameters )
    type is ( real ( KDR ) )    
    
    Fermi_2  =  max ( Fermi_2, tiny ( 0.0_KDR ) )

    Result  =  Factor  *  Fermi_3 ** ( 1.0_KDR / 4.0_KDR )  &
                  *  Fermi_2 ** ( - 1.0_KDR / 3.0_KDR ) &
               -  P

    end select
  
  end subroutine EvaluateZeroEta

  
end module NeutrinoMoments_G__Form
