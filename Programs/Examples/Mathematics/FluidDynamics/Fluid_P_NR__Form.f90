module Fluid_P_NR__Form

  !-- Fluid_Perfect_NonRelativistic__Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_NON_RELATIVISTIC = 1, &
      N_CONSERVED_NON_RELATIVISTIC = 0, &
      N_FIELDS_NON_RELATIVISTIC    = 1, &
      N_VECTORS_NON_RELATIVISTIC   = 0

  type, public, extends ( Fluid_P_Template ) :: Fluid_P_NR_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_NON_RELATIVISTIC = N_PRIMITIVE_NON_RELATIVISTIC, &
      N_CONSERVED_NON_RELATIVISTIC = N_CONSERVED_NON_RELATIVISTIC, &
      N_FIELDS_NON_RELATIVISTIC    = N_FIELDS_NON_RELATIVISTIC, &
      N_VECTORS_NON_RELATIVISTIC   = N_VECTORS_NON_RELATIVISTIC, &
      MEAN_MOLECULAR_WEIGHT = 0
    real ( KDR ) :: &
      AdiabaticIndex        = 1.4_KDR, &
      MeanMolecularWeight   = 1.0_KDR, &
      FiducialBaryonDensity = 1.0_KDR, &
      FiducialTemperature   = 1.0_KDR
  contains
    procedure, public, pass :: &
      InitializeAllocate_P_NR
    generic, public :: &
      Initialize => InitializeAllocate_P_NR
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetAdiabaticIndex
    procedure, public, pass :: &
      SetMeanMolecularWeight
    procedure, public, pass :: &
      SetFiducialBaryonDensity
    procedure, public, pass :: &
      SetFiducialTemperature
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass :: & 
      ComputeFromPrimitiveTemperature
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, nopass :: &
      Apply_EOS_NR_T_Kernel
    procedure, public, nopass :: &
      Apply_EOS_NR_E_Kernel
  end type Fluid_P_NR_Form

    private :: &
      InitializeBasics

contains


  subroutine InitializeAllocate_P_NR &
               ( F, RiemannSolverType, ReconstructedType, UseLimiter, &
                 VelocityUnit, MassDensityUnit, EnergyDensityUnit, &
                 TemperatureUnit, LimiterParameter, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, PinnedOption, &
                 UnitOption, VectorIndicesOption )

    class ( Fluid_P_NR_Form ), intent ( inout ) :: &
      F
    character ( * ), intent ( in ) :: &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
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
      ClearOption, &
      PinnedOption
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
           ( RiemannSolverType, ReconstructedType, UseLimiter, VelocityUnit, &
             MassDensityUnit, EnergyDensityUnit, TemperatureUnit, &
             LimiterParameter, nValues, VariableOption = Variable, & 
             VectorOption = VectorOption, NameOption = NameOption, & 
             ClearOption = ClearOption, PinnedOption = PinnedOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeAllocate_P_NR


  subroutine SetPrimitiveConserved ( C )

    class ( Fluid_P_NR_Form ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_NON_RELATIVISTIC ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_NON_RELATIVISTIC ) :: &
      ConservedName

    !-- select primitive, conserved

    oP = C % N_PRIMITIVE_TEMPLATE + C % N_PRIMITIVE_DUST &
         + C % N_PRIMITIVE_PERFECT
    oC = C % N_CONSERVED_TEMPLATE + C % N_CONSERVED_DUST &
         + C % N_CONSERVED_PERFECT

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_NON_RELATIVISTIC
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_NON_RELATIVISTIC ) &
      = [ C % INTERNAL_ENERGY ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_NON_RELATIVISTIC
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
!    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_NON_RELATIVISTIC ) &
!      = [ ]
    
    call C % SetPrimitiveConservedTemplate_P ( )

  end subroutine SetPrimitiveConserved


  subroutine SetAdiabaticIndex ( F, AdiabaticIndex )

    class ( Fluid_P_NR_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      AdiabaticIndex

    F % AdiabaticIndex = AdiabaticIndex

  end subroutine SetAdiabaticIndex


  subroutine SetMeanMolecularWeight ( F, MeanMolecularWeight )

    class ( Fluid_P_NR_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      MeanMolecularWeight

    F % MeanMolecularWeight = MeanMolecularWeight

  end subroutine SetMeanMolecularWeight


  subroutine SetFiducialBaryonDensity ( F, FiducialBaryonDensity )

    class ( Fluid_P_NR_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      FiducialBaryonDensity

    F % FiducialBaryonDensity = FiducialBaryonDensity

  end subroutine SetFiducialBaryonDensity


  subroutine SetFiducialTemperature ( F, FiducialTemperature )

    class ( Fluid_P_NR_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      FiducialTemperature

    F % FiducialTemperature = FiducialTemperature

  end subroutine SetFiducialTemperature


  subroutine SetOutput ( F, Output )

    class ( Fluid_P_NR_Form ), intent ( in ) :: &
      F
    type ( StorageForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( F % VELOCITY_U )
    call Output % Initialize &
           ( F, iaSelectedOption &
                  = [ F % COMOVING_DENSITY, F % VELOCITY_U, F % PRESSURE, &
                      F % SOUND_SPEED, F % MACH_NUMBER, F % TEMPERATURE, &
                      F % ENTROPY_PER_BARYON ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  subroutine ComputeFromPrimitiveTemperature &
               ( F, G, nValuesOption, oValueOption ) 

    class ( Fluid_P_NR_Form  ), intent ( inout ) :: &
      F
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption
    
    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( FV => F % Value, &
        GV => G % Value )

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
      ( FEP_1 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        M     => FV ( oV + 1 : oV + nV, F % BARYON_MASS ), &
        N     => FV ( oV + 1 : oV + nV, F % COMOVING_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, F % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, F % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, F % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, F % CONSERVED_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, F % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, F % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, F % MOMENTUM_DENSITY_D ( 3 ) ), &
        E     => FV ( oV + 1 : oV + nV, F % INTERNAL_ENERGY ), &
        G     => FV ( oV + 1 : oV + nV, F % CONSERVED_ENERGY ), &
        P     => FV ( oV + 1 : oV + nV, F % PRESSURE ), &
        Gamma => FV ( oV + 1 : oV + nV, F % ADIABATIC_INDEX ), &
        CS    => FV ( oV + 1 : oV + nV, F % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, F % MACH_NUMBER ), &
        SB    => FV ( oV + 1 : oV + nV, F % ENTROPY_PER_BARYON ), &
        T     => FV ( oV + 1 : oV + nV, F % TEMPERATURE ) )

    call F % ComputeBaryonMassKernel ( M )
    call F % Apply_EOS_NR_T_Kernel &
           ( P, Gamma, SB, E, M, N, T, F % AdiabaticIndex, &
             F % MeanMolecularWeight, F % FiducialBaryonDensity, &
             F % FiducialTemperature, CONSTANT % ATOMIC_MASS_UNIT, &
             CONSTANT % BOLTZMANN )
    !call F % ComputeDensityMomentumKernel &
    !       ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call F % ComputeFromPrimitiveKernel &
           ( M, D, S_1, S_2, S_3, G, N, V_1, V_2, V_3, E, &
             M_DD_22, M_DD_33 )
    call F % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromPrimitiveTemperature


  subroutine ComputeFromPrimitiveCommon &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_NR_Form ), intent ( in ) :: &
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
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ), &
        E     => FV ( oV + 1 : oV + nV, C % INTERNAL_ENERGY ), &
        G     => FV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P     => FV ( oV + 1 : oV + nV, C % PRESSURE ), &
        Gamma => FV ( oV + 1 : oV + nV, C % ADIABATIC_INDEX ), &
        CS    => FV ( oV + 1 : oV + nV, C % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, C % MACH_NUMBER ), &
        SB    => FV ( oV + 1 : oV + nV, C % ENTROPY_PER_BARYON ), &
        T     => FV ( oV + 1 : oV + nV, C % TEMPERATURE ) )

    !call C % ComputeBaryonMassKernel ( M )
    !call C % ComputeDensityMomentumKernel &
    !       ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % ComputeFromPrimitiveKernel &
           ( M, D, S_1, S_2, S_3, G, N, V_1, V_2, V_3, E, &
             M_DD_22, M_DD_33 )
    call C % Apply_EOS_NR_E_Kernel &
           ( P, Gamma, SB, T, M, N, E, C % AdiabaticIndex, &
             C % MeanMolecularWeight, C % FiducialBaryonDensity, &
             C % FiducialTemperature, CONSTANT % ATOMIC_MASS_UNIT, &
             CONSTANT % BOLTZMANN )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Storage_C, C, G, Storage_G, DetectFeaturesOption, &
                 nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_NR_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_G
    logical ( KDL ), intent ( in ), optional :: &
      DetectFeaturesOption
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
    logical ( KDL ) :: &
      DetectFeatures

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
      ( M_UU_22 => GV ( oV + 1 : oV + nV, G % METRIC_UU_22 ), &
        M_UU_33 => GV ( oV + 1 : oV + nV, G % METRIC_UU_33 ) )
    associate &
      ( FEP_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        M     => FV ( oV + 1 : oV + nV, C % BARYON_MASS ), &
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ), &
        E     => FV ( oV + 1 : oV + nV, C % INTERNAL_ENERGY ), &
        G     => FV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P     => FV ( oV + 1 : oV + nV, C % PRESSURE ), &
        Gamma => FV ( oV + 1 : oV + nV, C % ADIABATIC_INDEX ), &
        CS    => FV ( oV + 1 : oV + nV, C % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, C % MACH_NUMBER ), &
        SB    => FV ( oV + 1 : oV + nV, C % ENTROPY_PER_BARYON ), &
        T     => FV ( oV + 1 : oV + nV, C % TEMPERATURE ) )

    !call C % ComputeBaryonMassKernel ( M )
    !call C % ComputeDensityVelocityKernel &
    !       ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M,  )
    call C % ComputeFromConservedKernel &
           ( E, G, M, N, D, V_1, V_2, V_3, S_1, S_2, S_3, &
             M_UU_22, M_UU_33 )
    call C % Apply_EOS_NR_E_Kernel &
           ( P, Gamma, SB, T, M, N, E, C % AdiabaticIndex, &
             C % MeanMolecularWeight, C % FiducialBaryonDensity, &
             C % FiducialTemperature, CONSTANT % ATOMIC_MASS_UNIT, &
             CONSTANT % BOLTZMANN )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.
    
    DetectFeatures = .false.
    if ( present ( DetectFeaturesOption ) ) &
      DetectFeatures = DetectFeaturesOption
    if ( DetectFeatures ) &
      call C % Features % Detect ( )

  end subroutine ComputeFromConservedCommon


  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Storage_C, Storage_G, iDimension, &
                 nValuesOption, oValueOption )
    
    class ( StorageForm ), intent ( inout ) :: &
      RawFlux
    class ( Fluid_P_NR_Form ), intent ( in ) :: &
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


  subroutine Apply_EOS_NR_T_Kernel &
               ( P, Gamma, SB, E, M, N, T, Gamma_0, Mu_0, N_0, T_0, AMU, k_B )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      Gamma, &
      SB, &
      E
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      T
    real ( KDR ), intent ( in ) :: &
      Gamma_0, &
      Mu_0, &
      N_0, &
      T_0, &
      AMU, &
      k_B

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( P )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV )  >  0.0_KDR ) then
        P ( iV )  =  M ( iV )  *  N ( iV )  *  k_B  *  T ( iV )  &
                     /  ( Mu_0 * AMU )  
      else
        P ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      Gamma ( iV )  =  Gamma_0
      E     ( iV )  =  P ( iV )  /  ( Gamma_0 - 1.0_KDR ) 
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV )  >  0.0_KDR ) then
        SB ( iV )  &
          =  ( k_B / Mu_0 ) &
             *  log ( ( N_0  /  N ( iV ) ) &
                      *  ( T ( iV ) / T_0 ) &
                         ** ( 1.0_KDR / ( Gamma_0 - 1.0_KDR ) ) )
      else
        SB ( iV )  =  - 0.1_KDR * huge ( 1.0_KDR )
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine Apply_EOS_NR_T_Kernel


  subroutine Apply_EOS_NR_E_Kernel &
               ( P, Gamma, SB, T, M, N, E, Gamma_0, Mu_0, N_0, T_0, AMU, k_B )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      Gamma, &
      SB, &
      T
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      E
    real ( KDR ), intent ( in ) :: &
      Gamma_0, &
      Mu_0, &
      N_0, &
      T_0, &
      AMU, &
      k_B

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( P )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      Gamma ( iV )  =  Gamma_0
      P     ( iV )  =  E ( iV )  *  ( Gamma_0 - 1.0_KDR ) 
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV )  >  0.0_KDR ) then
        T ( iV )  =  P ( iV )  *  ( Mu_0 * AMU )  &
                     /  ( M ( iV )  *  N ( iV )  * k_B )
      else
        T ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV )  >  0.0_KDR ) then
        SB ( iV )  &
          =  ( k_B / Mu_0 ) &
             *  log ( ( N_0  /  N ( iV ) ) &
                      *  ( T ( iV ) / T_0 ) &
                         ** ( 1.0_KDR / ( Gamma_0 - 1.0_KDR ) ) )
      else
        SB ( iV )  =  - 0.1_KDR * huge ( 1.0_KDR )
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine Apply_EOS_NR_E_Kernel


  subroutine InitializeBasics &
               ( F, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( Fluid_P_NR_Form ), intent ( inout ) :: &
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
      F % Type = 'a Fluid_P_NR'

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE + F % N_FIELDS_DUST &
         + F % N_FIELDS_PERFECT
    if ( F % N_FIELDS == 0 ) &
      F % N_FIELDS = oF + F % N_FIELDS_NON_RELATIVISTIC

    F % MEAN_MOLECULAR_WEIGHT = oF + 1

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_NON_RELATIVISTIC ) &
      = [ 'MeanMolecularWeight' ]
    
    !-- units

    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( F % N_FIELDS ) )
      VariableUnit ( oF + 1 : oF + F % N_FIELDS_NON_RELATIVISTIC ) &
        = UNIT % IDENTITY
    end if

    !-- vectors

    oV = F % N_VECTORS_TEMPLATE + F % N_VECTORS_DUST &
         + F % N_VECTORS_PERFECT
    if ( F % N_VECTORS == 0 ) &
      F % N_VECTORS = oV + F % N_VECTORS_NON_RELATIVISTIC

  end subroutine InitializeBasics
  
    
end module Fluid_P_NR__Form
