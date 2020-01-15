module Fluid_P_I__Form

  !-- Fluid_Perfect_Ideal__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use Fluid_P__Template
  use FluidFeatures_P__Form
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_IDEAL = 0, &
      N_CONSERVED_IDEAL = 0, &
      N_FIELDS_IDEAL    = 0, &
      N_VECTORS_IDEAL   = 0

  type, public, extends ( Fluid_P_Template ) :: Fluid_P_I_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_IDEAL = N_PRIMITIVE_IDEAL, &
      N_CONSERVED_IDEAL = N_CONSERVED_IDEAL, &
      N_FIELDS_IDEAL    = N_FIELDS_IDEAL, &
      N_VECTORS_IDEAL   = N_VECTORS_IDEAL
    real ( KDR ) :: &
      BoltzmannConstant, &
      AdiabaticIndex, &
      MeanMolecularWeight, &
      SpecificHeatVolume, &  !-- per baryon
      FiducialBaryonDensity, &
      FiducialPressure
  contains
    procedure, public, pass :: &
      Initialize_P_I
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetAdiabaticIndex
    procedure, public, pass :: &
      SetMeanMolecularWeight
    procedure, public, pass :: &
      SetSpecificHeatVolume
    procedure, public, pass :: &
      SetFiducialParameters
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass ( C ) :: &
      ComputeFromTemperature
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, nopass :: &
      Apply_EOS_I_T_Kernel
    procedure, public, nopass :: &
      Apply_EOS_I_SB_E_Kernel
    procedure, public, nopass :: &
      Apply_EOS_I_E_Kernel
  end type Fluid_P_I_Form

    private :: &
      InitializeBasics


contains


  subroutine Initialize_P_I &
               ( F, FluidType, RiemannSolverType, ReconstructedType, &
                 UseEntropy, UseLimiter, Units, BaryonMassReference, &
                 LimiterParameter, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    character ( * ), intent ( in ) :: &
      FluidType, &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseEntropy, &
      UseLimiter
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
    real ( KDR ), intent ( in ) :: &
      BaryonMassReference, &
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
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), &
      optional :: &
        VectorIndicesOption

    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit

    call InitializeBasics &
           ( F, Variable, VariableUnit, VariableOption, UnitOption )

    call F % InitializeTemplate_P &
           ( FluidType, RiemannSolverType, ReconstructedType, UseEntropy, &
             UseLimiter, Units, BaryonMassReference, LimiterParameter, &
             nValues, VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    associate &
      ( k     => F % BoltzmannConstant, &
        gamma => F % AdiabaticIndex, &
        mu    => F % MeanMolecularWeight, &
        c_v   => F % SpecificHeatVolume, &
        n_0   => F % FiducialBaryonDensity, &
        p_0   => F % FiducialPressure )

    if ( Units % Temperature /= UNIT % IDENTITY ) then
      k = CONSTANT % BOLTZMANN
      call Show ( k, UNIT % JOULE / UNIT % KELVIN, 'BoltzmannConstant', &
                  F % IGNORABILITY )
    else !-- Dimensionless
      k = 1.0_KDR
      call Show ( k, 'BoltzmannConstant', F % IGNORABILITY )
    end if

    gamma  =  1.4_KDR
    mu     =  1.0_KDR
    c_v    =  k / ( mu * ( gamma - 1.0_KDR ) )
    call Show ( gamma, 'AdiabaticIndex', F % IGNORABILITY )
    call Show ( mu, 'MeanMolecularWeight', F % IGNORABILITY )
    if ( Units % Temperature /= UNIT % IDENTITY ) then
      call Show ( c_v, UNIT % BOLTZMANN, 'SpecificHeatVolume', &
                  F % IGNORABILITY )
    else
      call Show ( c_v, 'SpecificHeatVolume', F % IGNORABILITY )
    end if

    n_0  =  1.0_KDR
    p_0  =  1.0_KDR
    call Show ( n_0, Units % NumberDensity, 'FiducialBaryonDensity', &
                F % IGNORABILITY )
    call Show ( p_0, Units % EnergyDensity, 'FiducialPressure', &
                F % IGNORABILITY )

    end associate !-- amu, etc.

  end subroutine Initialize_P_I


  subroutine SetPrimitiveConserved ( C )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_IDEAL ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_IDEAL ) :: &
      ConservedName

    oP = C % N_PRIMITIVE_TEMPLATE + C % N_PRIMITIVE_DUST &
         + C % N_PRIMITIVE_PERFECT
    oC = C % N_CONSERVED_TEMPLATE + C % N_CONSERVED_DUST &
         + C % N_CONSERVED_PERFECT

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_IDEAL
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
!    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_IDEAL ) &
!      = [ ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_IDEAL
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
!    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_IDEAL ) &
!      = [ ]
    
    do iF = 1, C % N_PRIMITIVE_IDEAL
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_IDEAL
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )

    call C % SetPrimitiveConservedTemplate_P ( )

  end subroutine SetPrimitiveConserved


  subroutine SetAdiabaticIndex ( F, AdiabaticIndex )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      AdiabaticIndex

    F % AdiabaticIndex = AdiabaticIndex

    call Show ( 'Setting AdiabaticIndex of a Fluid_P_I', F % IGNORABILITY )
    call Show ( F % Name, 'Name', F % IGNORABILITY )
    call Show ( F % AdiabaticIndex, 'AdiabaticIndex', F % IGNORABILITY )

  end subroutine SetAdiabaticIndex


  subroutine SetMeanMolecularWeight ( F, MeanMolecularWeight )

    !-- Assumes AdiabaticIndex already set.

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      MeanMolecularWeight

    associate &
      ( k     => F % BoltzmannConstant, &
        gamma => F % AdiabaticIndex, &
        mu    => F % MeanMolecularWeight, &
        c_v   => F % SpecificHeatVolume )

    mu  =  MeanMolecularWeight

    c_v  =  k / ( mu * ( gamma - 1.0_KDR ) )

    call Show ( 'Setting MeanMolecularWeight of a Fluid_P_I', F % IGNORABILITY )
    call Show ( F % Name, 'Name', F % IGNORABILITY )
    call Show ( F % MeanMolecularWeight, 'MeanMolecularWeight', &
                F % IGNORABILITY )

    if ( F % Unit ( F % TEMPERATURE ) /= UNIT % IDENTITY ) then
      call Show ( c_v, UNIT % BOLTZMANN, 'SpecificHeatVolume', &
                  F % IGNORABILITY )
    else
      call Show ( c_v, 'SpecificHeatVolume', F % IGNORABILITY )
    end if

    end associate !-- amu, etc.

  end subroutine SetMeanMolecularWeight


  subroutine SetSpecificHeatVolume ( F, SpecificHeatVolume )

    !-- Assumes AdiabaticIndex already set.

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      SpecificHeatVolume

    associate &
      ( k     => F % BoltzmannConstant, &
        gamma => F % AdiabaticIndex, &
        mu    => F % MeanMolecularWeight, &
        c_v   => F % SpecificHeatVolume )

    c_v  =  SpecificHeatVolume

    mu  =  k / ( c_v * ( gamma - 1.0_KDR ) )

    call Show ( 'Setting SpecificHeatVolume of a Fluid_P_I', F % IGNORABILITY )
    call Show ( F % Name, 'Name', F % IGNORABILITY )
    if ( F % Unit ( F % TEMPERATURE ) /= UNIT % IDENTITY ) then
      call Show ( c_v, UNIT % BOLTZMANN, 'SpecificHeatVolume', &
                  F % IGNORABILITY )
    else
      call Show ( c_v, 'SpecificHeatVolume', F % IGNORABILITY )
    end if

    call Show ( F % MeanMolecularWeight, 'MeanMolecularWeight', &
                F % IGNORABILITY )

    end associate !-- amu, etc.

  end subroutine SetSpecificHeatVolume


  subroutine SetFiducialParameters &
               ( F, FiducialBaryonDensity, FiducialPressure )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      FiducialBaryonDensity, &
      FiducialPressure

    F % FiducialBaryonDensity = FiducialBaryonDensity
    F % FiducialPressure = FiducialPressure

    call Show ( 'Setting fiducial parameters of a Fluid_P_I', F % IGNORABILITY )
    call Show ( F % Name, 'Name', F % IGNORABILITY )
    call Show ( F % FiducialBaryonDensity, &
                F % Unit ( F % COMOVING_BARYON_DENSITY ), &
                'FiducialBaryonDensity', F % IGNORABILITY )
    call Show ( F % FiducialPressure, F % Unit ( F % PRESSURE ), &
                'FiducialPressure', F % IGNORABILITY )

  end subroutine SetFiducialParameters


  subroutine SetOutput ( F, Output )

    class ( Fluid_P_I_Form ), intent ( in ) :: &
      F
    type ( StorageForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( F % VELOCITY_U )
    call Output % Initialize &
           ( F, iaSelectedOption &
                  = [ F % COMOVING_BARYON_DENSITY, F % VELOCITY_U, &
                      F % PRESSURE, F % TEMPERATURE, F % MACH_NUMBER, &
                      F % ENTROPY_PER_BARYON ], &
             VectorOption = [ 'Velocity' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  subroutine ComputeFromTemperature &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_I_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( FV => Storage_C % Value, &
        GV => Storage_G % Value )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( FV, dim = 1 )
    end if

    associate &
      ( M_DD_22 => GV ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => GV ( oV + 1 : oV + nV, G % METRIC_DD_33 ), &
        M_UU_22 => GV ( oV + 1 : oV + nV, G % METRIC_UU_22 ), &
        M_UU_33 => GV ( oV + 1 : oV + nV, G % METRIC_UU_33 ) )
    associate &
      ( FEP_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        M     => FV ( oV + 1 : oV + nV, C % BARYON_MASS ), &
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_BARYON_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_BARYON_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ), &
        E     => FV ( oV + 1 : oV + nV, C % INTERNAL_ENERGY ), &
        G     => FV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P     => FV ( oV + 1 : oV + nV, C % PRESSURE ), &
        T     => FV ( oV + 1 : oV + nV, C % TEMPERATURE ), &
        SB    => FV ( oV + 1 : oV + nV, C % ENTROPY_PER_BARYON ), &
        DS    => FV ( oV + 1 : oV + nV, C % CONSERVED_ENTROPY ), &
        CS    => FV ( oV + 1 : oV + nV, C % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, C % MACH_NUMBER ) )

    call C % Compute_M_Kernel ( M, C % BaryonMassReference )
    call C % Apply_EOS_I_T_Kernel &
           ( P, E, SB, CS, M, N, T, C % AdiabaticIndex, &
             C % SpecificHeatVolume, C % FiducialBaryonDensity, &
             C % FiducialPressure )
    call C % Compute_D_S_G_Kernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % Compute_G_G_Kernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
    call C % Compute_DS_G_Kernel &
           ( DS, N, SB )
    call C % Compute_FE_P_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
             V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.

    if ( associated ( C % Value, Storage_C % Value ) ) &
      call C % Features % Detect ( )

  end subroutine ComputeFromTemperature


  subroutine ComputeFromPrimitiveCommon &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_I_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( FV => Storage_C % Value, &
        GV => Storage_G % Value )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( FV, dim = 1 )
    end if

    associate &
      ( M_DD_22 => GV ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => GV ( oV + 1 : oV + nV, G % METRIC_DD_33 ), &
        M_UU_22 => GV ( oV + 1 : oV + nV, G % METRIC_UU_22 ), &
        M_UU_33 => GV ( oV + 1 : oV + nV, G % METRIC_UU_33 ) )
    associate &
      ( FEP_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        M     => FV ( oV + 1 : oV + nV, C % BARYON_MASS ), &
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_BARYON_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_BARYON_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ), &
        E     => FV ( oV + 1 : oV + nV, C % INTERNAL_ENERGY ), &
        G     => FV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P     => FV ( oV + 1 : oV + nV, C % PRESSURE ), &
        T     => FV ( oV + 1 : oV + nV, C % TEMPERATURE ), &
        SB    => FV ( oV + 1 : oV + nV, C % ENTROPY_PER_BARYON ), &
        DS    => FV ( oV + 1 : oV + nV, C % CONSERVED_ENTROPY ), &
        CS    => FV ( oV + 1 : oV + nV, C % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, C % MACH_NUMBER ) )

    call C % Compute_M_Kernel ( M, C % BaryonMassReference )
    call C % Apply_EOS_I_E_Kernel &
           ( P, T, SB, CS, M, N, E, C % AdiabaticIndex, &
             C % SpecificHeatVolume, C % FiducialBaryonDensity, &
             C % FiducialPressure )
    call C % Compute_D_S_G_Kernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % Compute_G_G_Kernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
    call C % Compute_DS_G_Kernel &
           ( DS, N, SB )
    call C % Compute_FE_P_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
             V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.

    if ( associated ( C % Value, Storage_C % Value ) ) &
      call C % Features % Detect ( )

  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_I_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues

    associate &
      ( FV => Storage_C % Value, &
        GV => Storage_G % Value )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( FV, dim = 1 )
    end if
    
    select type ( FF => C % Features )
    class is ( FluidFeatures_P_Form )

    associate &
      ( M_DD_22 => GV ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => GV ( oV + 1 : oV + nV, G % METRIC_DD_33 ), &
        M_UU_22 => GV ( oV + 1 : oV + nV, G % METRIC_UU_22 ), &
        M_UU_33 => GV ( oV + 1 : oV + nV, G % METRIC_UU_33 ) )
    associate &
      ( FEP_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        M     => FV ( oV + 1 : oV + nV, C % BARYON_MASS ), &
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_BARYON_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_BARYON_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ), &
        E     => FV ( oV + 1 : oV + nV, C % INTERNAL_ENERGY ), &
        G     => FV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P     => FV ( oV + 1 : oV + nV, C % PRESSURE ), &
        T     => FV ( oV + 1 : oV + nV, C % TEMPERATURE ), &
        SB    => FV ( oV + 1 : oV + nV, C % ENTROPY_PER_BARYON ), &
        DS    => FV ( oV + 1 : oV + nV, C % CONSERVED_ENTROPY ), &
        CS    => FV ( oV + 1 : oV + nV, C % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, C % MACH_NUMBER ), &
        Shock => FF % Value ( oV + 1 : oV + nV, FF % SHOCK ) )

    call C % Compute_M_Kernel &
           ( M, C % BaryonMassReference )
    call C % Compute_N_V_E_G_Kernel &
               ( N, V_1, V_2, V_3, E, D, S_1, S_2, S_3, G, M, &
                 M_UU_22, M_UU_33, C % BaryonDensityMin )
    if ( C % UseEntropy ) then
      call C % Compute_SB_G_Kernel &
             ( SB, DS, N )
      call C % Apply_EOS_I_SB_E_Kernel &
             ( P, E, T, SB, CS, M, N, Shock, C % AdiabaticIndex, &
               C % SpecificHeatVolume, C % FiducialBaryonDensity, &
               C % FiducialPressure )
      call C % Compute_G_G_Kernel &
             ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
    else
      call C % Apply_EOS_I_E_Kernel &
             ( P, T, SB, CS, M, N, E, C % AdiabaticIndex, &
               C % SpecificHeatVolume, C % FiducialBaryonDensity, &
               C % FiducialPressure )
    end if
    call C % Compute_DS_G_Kernel &
           ( DS, N, SB )
    call C % Compute_FE_P_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
             V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end select !-- FF
    end associate !-- FV, etc.
    
    if ( associated ( C % Value, Storage_C % Value ) ) &
      call C % Features % Detect ( )

  end subroutine ComputeFromConservedCommon

  
  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Storage_C, Storage_G, iDimension, &
                 nValuesOption, oValueOption )
    
    class ( StorageForm ), intent ( inout ) :: &
      RawFlux
    class ( Fluid_P_I_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_C, &
      Storage_G
    integer ( KDI ), intent ( in ) :: &
      iDimension
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeRawFluxesTemplate_P &
           ( RawFlux, G, Storage_C, Storage_G, iDimension, nValuesOption, &
             oValueOption )

  end subroutine ComputeRawFluxes


  subroutine Apply_EOS_I_T_Kernel &
               ( P, E, SB, CS, M, N, T, Gamma, C_V, N0, P0 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      E, &
      SB, &
      CS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      T
    real ( KDR ), intent ( in ) :: &
      Gamma, &
      C_V, &
      N0, &
      P0

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      SqrtHuge

    SqrtHuge = sqrt ( huge ( 1.0_KDR ) )

    nValues = size ( P )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      E ( iV )  =  C_V  *  N ( iV )  *  T ( iV )

      P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

      if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then

        SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                    *  ( N0 / N ( iV ) ) ** Gamma ) 

        CS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

      else
        SB ( iV )  =  - SqrtHuge
        CS ( iV )  =    0.0_KDR
      end if

    end do !-- iV
    !$OMP end parallel do

  end subroutine Apply_EOS_I_T_Kernel


  subroutine Apply_EOS_I_SB_E_Kernel &
               ( P, E, T, SB, CS, M, N, Shock, Gamma, C_V, N0, P0 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      E, &
      T, &
      SB, &
      CS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      Shock
    real ( KDR ), intent ( in ) :: &
      Gamma, &
      C_V, &
      N0, &
      P0

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      SqrtHuge

    SqrtHuge = sqrt ( huge ( 1.0_KDR ) )

    nValues = size ( P )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      if ( Shock ( iV ) > 0.0_KDR ) then

        P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
          SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                      *  ( N0 / N ( iV ) ) ** Gamma ) 
        else
          SB ( iV )  =  - SqrtHuge
        end if

      else

        P ( iV )  =  P0  *  ( N ( iV ) / N0 ) ** Gamma  &
                     *  exp ( SB ( iV ) / C_V )

        E ( iV )  =  P ( iV )  /  ( Gamma - 1.0_KDR )

      end if

      if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then

         T ( iV )  =  E ( iV )  /  ( C_V  *  N ( iV ) )

        CS ( iV ) = sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

      else
         T ( iV ) = 0.0_KDR
        CS ( iV ) = 0.0_KDR
      end if 

    end do !-- iV
    !$OMP end parallel do

  end subroutine Apply_EOS_I_SB_E_Kernel


  subroutine Apply_EOS_I_E_Kernel ( P, T, SB, CS, M, N, E, Gamma, C_V, N0, P0 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      T, &
      SB, &
      CS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      E
    real ( KDR ), intent ( in ) :: &
      Gamma, &
      C_V, &
      N0, &
      P0

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      SqrtHuge

    SqrtHuge = sqrt ( huge ( 1.0_KDR ) )

    nValues = size ( P )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

      if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then

         T ( iV )  =  E ( iV )  /  ( C_V  *  N ( iV ) )

        SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                    *  ( N0 / N ( iV ) ) ** Gamma ) 

        CS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

      else
         T ( iV )  =    0.0_KDR
        SB ( iV )  =  - SqrtHuge
        CS ( iV )  =    0.0_KDR
      end if

    end do !-- iV
    !$OMP end parallel do

  end subroutine Apply_EOS_I_E_Kernel


  subroutine InitializeBasics &
               ( F, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
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
      oV, &  !-- oVector
      oP, &  !-- oPrimitive
      oC     !-- oConserved

    if ( F % Type == '' ) &
      F % Type = 'a Fluid_P_I'

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE + F % N_FIELDS_DUST &
         + F % N_FIELDS_PERFECT
    if ( F % N_FIELDS == 0 ) &
      F % N_FIELDS = oF + F % N_FIELDS_IDEAL

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    !-- units

    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( F % N_FIELDS ) )
      VariableUnit ( oF + 1 : oF + F % N_FIELDS_IDEAL ) = UNIT % IDENTITY
    end if

    !-- vectors

    oV = F % N_VECTORS_TEMPLATE + F % N_VECTORS_DUST &
         + F % N_VECTORS_PERFECT
    if ( F % N_VECTORS == 0 ) F % N_VECTORS = oV + F % N_VECTORS_IDEAL

  end subroutine InitializeBasics
  
    
end module Fluid_P_I__Form
