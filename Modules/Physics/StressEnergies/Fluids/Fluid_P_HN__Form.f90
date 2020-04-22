module Fluid_P_HN__Form

  !-- Fluid_Perfect_RepresentativeHeavyNucleus__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
!  use FluidFeatures_Template
  use Fluid_P__Template
  use FluidFeatures_P__Form
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_HEAVY_NUCLEUS = 1, &
      N_CONSERVED_HEAVY_NUCLEUS = 1, &
      N_FIELDS_HEAVY_NUCLEUS    = 10, &
      N_VECTORS_HEAVY_NUCLEUS   = 0

  type, public, extends ( Fluid_P_Template ) :: Fluid_P_HN_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_HEAVY_NUCLEUS   = N_PRIMITIVE_HEAVY_NUCLEUS, &
      N_CONSERVED_HEAVY_NUCLEUS   = N_CONSERVED_HEAVY_NUCLEUS, &
      N_FIELDS_HEAVY_NUCLEUS      = N_FIELDS_HEAVY_NUCLEUS, &
      N_VECTORS_HEAVY_NUCLEUS     = N_VECTORS_HEAVY_NUCLEUS, &
      ELECTRON_FRACTION           = 0, &
      CONSERVED_ELECTRON_DENSITY  = 0, &
      MASS_FRACTION_PROTON        = 0, &
      MASS_FRACTION_NEUTRON       = 0, &
      MASS_FRACTION_ALPHA         = 0, &
      MASS_FRACTION_HEAVY         = 0, &
      ATOMIC_NUMBER_HEAVY         = 0, &
      MASS_NUMBER_HEAVY           = 0, &
      CHEMICAL_POTENTIAL_N_P      = 0, &  
        !-- a.k.a. mu_hat. Includes m_n - m_p. (mu_n and mu_p both
        !   measured with respect to m_n.)
      CHEMICAL_POTENTIAL_E      = 0     
        !-- Includes m_e.
  contains
    procedure, public, pass :: &
      Initialize_P_HN
    procedure, public, pass :: &
      SetPrimitiveConserved
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
    procedure, public, pass ( C ) :: &
      ComputeCenterStates
    procedure, public, nopass :: &
      Compute_DE_G_Kernel
    procedure, public, nopass :: &
      Compute_YE_G_Kernel
    procedure, public, nopass :: &
      Apply_EOS_HN_T_Kernel
    procedure, public, nopass :: &
      Apply_EOS_HN_SB_E_Kernel
    procedure, public, nopass :: &
      Apply_EOS_HN_E_Kernel
    procedure, public, nopass :: &
      Apply_EOS_HN_SB_Kernel
  end type Fluid_P_HN_Form

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeRawFluxesKernel, &
      ComputeCenterStatesKernel!, &
!      InterpolateSoundSpeed

    real ( KDR ), private, protected :: &
      OR_Shift, &
      MassDensity_CGS, &
      SpecificEnergy_CGS, &
      Pressure_CGS, &
      Speed_CGS, &
      MeV
    logical ( KDL ), private, protected :: &
      TableInitialized = .false.
      
  interface 
  
    module subroutine Compute_DE_G_Kernel ( DE, N, YE, UseDeviceOption )
     use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
        DE
      real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
        N, &
        YE
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_DE_G_Kernel


    module subroutine Compute_YE_G_Kernel ( YE, DE, N, UseDeviceOption )
      !-- Compute_ProtonFraction_Galilean_Kernel
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        YE, &
        DE
      real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
        N 	 	 
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_YE_G_Kernel
    
    
    module subroutine Apply_EOS_HN_T_Kernel &
             ( P, E, CS, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               M, N, T, YE, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        P, &
        E, &
        CS, &
        SB, &
        X_P, X_N, X_He, X_A, &
        Z, A, &
        Mu_NP, Mu_E
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N, &
        T, &
        YE
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_HN_T_Kernel
    

    module subroutine Apply_EOS_HN_SB_E_Kernel &
             ( P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               M, N, YE, Shock, UseDeviceOption ) 
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        P, &
        T, &
        CS, &
        E, &
        SB, &
        X_P, X_N, X_He, X_A, &
        Z, A, &
        Mu_NP, Mu_E
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N, &
        YE, &
        Shock
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_HN_SB_E_Kernel


    module subroutine Apply_EOS_HN_E_Kernel &
             ( P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               M, N, YE, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        P, &
        T, &
        CS, &
        E, &
        SB, &
        X_P, X_N, X_He, X_A, &
        Z, A, &
        Mu_NP, Mu_E
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N, &
        YE
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_HN_E_Kernel


    module subroutine Apply_EOS_HN_SB_Kernel &
             ( P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               M, N, YE, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        P, &
        T, &
        CS, &
        E, &
        SB, &
        X_P, X_N, X_He, X_A, &
        Z, A, &
        Mu_NP, Mu_E
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N, &
        YE
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_HN_SB_Kernel


    module subroutine ComputeRawFluxesKernel &
             ( F_DE, DE, V_Dim, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        F_DE
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        DE, &
        V_Dim
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeRawFluxesKernel


    module subroutine ComputeCenterStatesKernel &
             ( DE_ICL, DE_ICR, DE_IL, DE_IR, V_Dim_IL, V_Dim_IR, &
               AP_I, AM_I, AC_I, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        DE_ICL, DE_ICR
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        DE_IL, DE_IR, &
        V_Dim_IL, V_Dim_IR, &
        AP_I, &
        AM_I, &
        AC_I
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeCenterStatesKernel

  end interface

contains


  subroutine Initialize_P_HN &
               ( F, FluidType, RiemannSolverType, ReconstructedType, &
                 UseEntropy, UseLimiter, Units, BaryonMassReference, &
                 LimiterParameter, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
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

    call SetUnits ( VariableUnit, F, Units )

    call F % InitializeTemplate_P &
           ( FluidType, RiemannSolverType, ReconstructedType, UseEntropy, &
             UseLimiter, Units, BaryonMassReference, LimiterParameter, &
             nValues, VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    if ( TableInitialized ) &
      return

    call DelayFileAccess ( PROGRAM_HEADER % Communicator % Rank )
    call READTABLE &
           ( '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5' )
!    call READTABLE &
!           ( '../Parameters/HShenEOS_rho220_temp180_ye65_version_1.1' &
!             // '_20120817.h5' )
!    call READTABLE &
!           ( '../Parameters/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1' &
!             // '_20120817.h5' )

    !-- Historical Oak Ridge Shift, accounting for nuclear binding energy
    OR_Shift = 8.9_KDR * UNIT % MEGA_ELECTRON_VOLT &
               / CONSTANT % ATOMIC_MASS_UNIT
    
    MassDensity_CGS     =  UNIT % MASS_DENSITY_CGS
    SpecificEnergy_CGS  =  UNIT % ERG  /  UNIT % GRAM
    Pressure_CGS        =  UNIT % BARYE
    Speed_CGS           =  UNIT % CENTIMETER  /  UNIT % SECOND
    MeV                 =  UNIT % MEGA_ELECTRON_VOLT

    TableInitialized  =  .true.

  end subroutine Initialize_P_HN


  subroutine SetPrimitiveConserved ( C )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_HEAVY_NUCLEUS ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_HEAVY_NUCLEUS ) :: &
      ConservedName

    oP = C % N_PRIMITIVE_TEMPLATE + C % N_PRIMITIVE_DUST &
         + C % N_PRIMITIVE_PERFECT
    oC = C % N_CONSERVED_TEMPLATE + C % N_CONSERVED_DUST &
         + C % N_CONSERVED_PERFECT

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_HEAVY_NUCLEUS
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_HEAVY_NUCLEUS ) &
      = [ C % ELECTRON_FRACTION ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_HEAVY_NUCLEUS
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_HEAVY_NUCLEUS ) &
      = [ C % CONSERVED_ELECTRON_DENSITY ]
    
    do iF = 1, C % N_PRIMITIVE_HEAVY_NUCLEUS
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_HEAVY_NUCLEUS
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )

    call C % SetPrimitiveConservedTemplate_P ( )

  end subroutine SetPrimitiveConserved


  subroutine InitializeBasics &
               ( F, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
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
      F % Type = 'a Fluid_P_HN'

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE + F % N_FIELDS_DUST &
         + F % N_FIELDS_PERFECT
    if ( F % N_FIELDS == 0 ) F % N_FIELDS &
      = oF + F % N_FIELDS_HEAVY_NUCLEUS

    F % ELECTRON_FRACTION           =  oF +  1
    F % CONSERVED_ELECTRON_DENSITY  =  oF +  2
    F % MASS_FRACTION_PROTON        =  oF +  3
    F % MASS_FRACTION_NEUTRON       =  oF +  4
    F % MASS_FRACTION_ALPHA         =  oF +  5
    F % MASS_FRACTION_HEAVY         =  oF +  6
    F % ATOMIC_NUMBER_HEAVY         =  oF +  7
    F % MASS_NUMBER_HEAVY           =  oF +  8
    F % CHEMICAL_POTENTIAL_N_P      =  oF +  9
    F % CHEMICAL_POTENTIAL_E        =  oF + 10

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_HEAVY_NUCLEUS ) &
      = [ 'ElectronFraction        ', &
          'ConservedElectronDensity', &
          'MassFractionProton      ', &
          'MassFractionNeutron     ', &
          'MassFractionAlpha       ', &
          'MassFractionHeavy       ', &
          'AtomicNumberHeavy       ', &
          'MassNumberHeavy         ', &
          'ChemicalPotential_N_P   ', &
          'ChemicalPotential_E     ' ]
    
    !-- units

    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( F % N_FIELDS ) )
      VariableUnit &
        ( oF + 1 : oF + F % N_FIELDS_HEAVY_NUCLEUS ) &
        = UNIT % IDENTITY
    end if

    !-- vectors

    oV = F % N_VECTORS_TEMPLATE + F % N_VECTORS_DUST &
         + F % N_VECTORS_PERFECT
    if ( F % N_VECTORS == 0 ) F % N_VECTORS &
      = oV + F % N_VECTORS_HEAVY_NUCLEUS

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, F, Units )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      F
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units

    VariableUnit ( F % CONSERVED_ELECTRON_DENSITY ) &
      = Units % NumberDensity
    VariableUnit ( F % CHEMICAL_POTENTIAL_N_P ) &
      = Units % Temperature
    VariableUnit ( F % CHEMICAL_POTENTIAL_E ) &
      = Units % Temperature

  end subroutine SetUnits


  subroutine SetOutput ( F, Output )

    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      F
    type ( StorageForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( F % VELOCITY_U )
    call Output % Initialize &
           ( F, iaSelectedOption &
                  = [ F % COMOVING_BARYON_DENSITY, F % VELOCITY_U, &
                      F % PRESSURE, F % MACH_NUMBER, F % TEMPERATURE, &
                      F % ENTROPY_PER_BARYON, F % ELECTRON_FRACTION, &
                      F % MASS_FRACTION_PROTON, F % MASS_FRACTION_NEUTRON, &
                      F % MASS_FRACTION_ALPHA, F % MASS_FRACTION_HEAVY, &
                      F % ATOMIC_NUMBER_HEAVY, F % MASS_NUMBER_HEAVY, &
                      F % CHEMICAL_POTENTIAL_N_P, &
                      F % CHEMICAL_POTENTIAL_E ], &
             VectorOption = [ 'Velocity' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput
  
  
  subroutine ComputeFromTemperature &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
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
    !    Gamma => FV ( oV + 1 : oV + nV, C % ADIABATIC_INDEX ), &
        CS    => FV ( oV + 1 : oV + nV, C % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, C % MACH_NUMBER ), &
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ), &
        DE    => FV ( oV + 1 : oV + nV, C % CONSERVED_ELECTRON_DENSITY ), &
        X_P   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_PROTON ), &
        X_N   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_NEUTRON ), &
        X_He  => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_ALPHA ), &
        X_A   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_HEAVY ), &
        Z     => FV ( oV + 1 : oV + nV, C % ATOMIC_NUMBER_HEAVY ), &
        A     => FV ( oV + 1 : oV + nV, C % MASS_NUMBER_HEAVY ), &
        Mu_NP => FV ( oV + 1 : oV + nV, C % CHEMICAL_POTENTIAL_N_P ), &
        Mu_E  => FV ( oV + 1 : oV + nV, C % CHEMICAL_POTENTIAL_E ) )

    call C % Compute_M_Kernel &
           ( M, C % BaryonMassReference, &
             UseDeviceOption = C % AllocatedDevice )
    call C % Apply_EOS_HN_T_Kernel &
           ( P, E, CS, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
             M, N, T, YE, UseDeviceOption = C % AllocatedDevice )
    call C % Compute_D_S_G_Kernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33, &
             UseDeviceOption = C % AllocatedDevice )
    call C % Compute_G_G_Kernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E, &
             UseDeviceOption = C % AllocatedDevice )
    call C % Compute_DS_G_Kernel &
           ( DS, N, SB, UseDeviceOption = C % AllocatedDevice )
    call C % Compute_DE_G_Kernel &
           ( DE, N, YE, UseDeviceOption = C % AllocatedDevice )
    call C % Compute_FE_P_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
             V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
             UseDeviceOption = C % AllocatedDevice )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromTemperature


  subroutine ComputeFromPrimitiveCommon &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
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
    !    Gamma => FV ( oV + 1 : oV + nV, C % ADIABATIC_INDEX ), &
        CS    => FV ( oV + 1 : oV + nV, C % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, C % MACH_NUMBER ), &
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ), &
        DE    => FV ( oV + 1 : oV + nV, C % CONSERVED_ELECTRON_DENSITY ), &
        X_P   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_PROTON ), &
        X_N   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_NEUTRON ), &
        X_He  => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_ALPHA ), &
        X_A   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_HEAVY ), &
        Z     => FV ( oV + 1 : oV + nV, C % ATOMIC_NUMBER_HEAVY ), &
        A     => FV ( oV + 1 : oV + nV, C % MASS_NUMBER_HEAVY ), &
        Mu_NP => FV ( oV + 1 : oV + nV, C % CHEMICAL_POTENTIAL_N_P ), &
        Mu_E  => FV ( oV + 1 : oV + nV, C % CHEMICAL_POTENTIAL_E ) )

    associate &
      ( T_CFP => PROGRAM_HEADER % Timer ( C % iTimerComputeFromPrimitive ), &
        T_CE  => PROGRAM_HEADER % Timer ( C % iTimerComputeEigenspeed ), &
        T_AE  => PROGRAM_HEADER % Timer ( C % iTimerApply_EOS ) )

    call T_CFP % Start ( )
    call Copy ( C % Value ( :, C % PRESSURE ), P, &
                UseDeviceOption = C % AllocatedDevice )
    call Copy ( C % Value ( :, C % TEMPERATURE ), T, &
                UseDeviceOption = C % AllocatedDevice )

    call C % Compute_M_Kernel &
           ( M, C % BaryonMassReference, &
             UseDeviceOption = C % AllocatedDevice )
    call T_CFP % Stop ( )
    
    call T_AE % Start ( )
!    call C % Apply_EOS_HN_T_Kernel &
!           ( P, E, CS, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
!             M, N, T, YE )
    call C % Apply_EOS_HN_E_Kernel &
           ( P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
             M, N, YE, UseDeviceOption = C % AllocatedDevice )
!    call C % Apply_EOS_HN_SB_Kernel &
!           ( P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
!             M, N, YE )
    call T_AE % Stop ( )
    
    call T_CFP % Start ( )
    call C % Compute_D_S_G_Kernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33, &
             UseDeviceOption = C % AllocatedDevice )
    call C % Compute_G_G_Kernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E, &
             UseDeviceOption = C % AllocatedDevice )
    call C % Compute_DS_G_Kernel &
           ( DS, N, SB, UseDeviceOption = C % AllocatedDevice )
    call C % Compute_DE_G_Kernel &
           ( DE, N, YE, UseDeviceOption = C % AllocatedDevice )
    call T_CFP % Stop ( )
    
    call T_CE % Start ( )
    call C % Compute_FE_P_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
             V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
             UseDeviceOption = C % AllocatedDevice )
    call T_CE % Stop ( )
    
    end associate !-- T_CFP, etc.

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Storage_C, C, G, Storage_G, DetectFeaturesOption, &
                 nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
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
        GE    => FV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P     => FV ( oV + 1 : oV + nV, C % PRESSURE ), &
    !     Gamma => FV ( oV + 1 : oV + nV, C % ADIABATIC_INDEX ), &
        CS    => FV ( oV + 1 : oV + nV, C % SOUND_SPEED ), &
        MN    => FV ( oV + 1 : oV + nV, C % MACH_NUMBER ), &
        T     => FV ( oV + 1 : oV + nV, C % TEMPERATURE ), &
        SB    => FV ( oV + 1 : oV + nV, C % ENTROPY_PER_BARYON ), &
        DS    => FV ( oV + 1 : oV + nV, C % CONSERVED_ENTROPY ), &
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ), &
        DE    => FV ( oV + 1 : oV + nV, C % CONSERVED_ELECTRON_DENSITY ), &
        X_P   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_PROTON ), &
        X_N   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_NEUTRON ), &
        X_He  => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_ALPHA ), &
        X_A   => FV ( oV + 1 : oV + nV, C % MASS_FRACTION_HEAVY ), &
        Z     => FV ( oV + 1 : oV + nV, C % ATOMIC_NUMBER_HEAVY ), &
        A     => FV ( oV + 1 : oV + nV, C % MASS_NUMBER_HEAVY ), &
        Mu_NP => FV ( oV + 1 : oV + nV, C % CHEMICAL_POTENTIAL_N_P ), &
        Mu_E  => FV ( oV + 1 : oV + nV, C % CHEMICAL_POTENTIAL_E ), &
        Shock => FF % Value ( oV + 1 : oV + nV, FF % SHOCK ) )

    associate &
      ( T_CFP => PROGRAM_HEADER % Timer ( C % iTimerComputeFromPrimitive ), &
        T_CE  => PROGRAM_HEADER % Timer ( C % iTimerComputeEigenspeed ), &
        T_AE  => PROGRAM_HEADER % Timer ( C % iTimerApply_EOS ) )
    
    call T_CFP % Start ( )
    call Copy ( C % Value ( :, C % PRESSURE ), P, &
                UseDeviceOption = C % AllocatedDevice )
    call Copy ( C % Value ( :, C % TEMPERATURE ), T, &
                UseDeviceOption = C % AllocatedDevice )

    call C % Compute_M_Kernel &
           ( M, C % BaryonMassReference, &
             UseDeviceOption = C % AllocatedDevice )
    call C % Compute_N_V_E_G_Kernel &
               ( N, V_1, V_2, V_3, E, D, S_1, S_2, S_3, GE, M, &
                 M_UU_22, M_UU_33, C % BaryonDensityMin, &
                 UseDeviceOption = C % AllocatedDevice )
    call C % Compute_YE_G_Kernel &
           ( YE, DE, N, UseDeviceOption = C % AllocatedDevice )
    call T_CFP % Stop ( )
    
    call T_AE % Start ( )
    if ( C % UseEntropy ) then
      call C % Compute_SB_G_Kernel &
             ( SB, DS, N, UseDeviceOption = C % AllocatedDevice )
      call C % Apply_EOS_HN_SB_E_Kernel &
             ( P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               M, N, YE, Shock, UseDeviceOption = C % AllocatedDevice )
      call C % Compute_G_G_Kernel &
             ( GE, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E, &
               UseDeviceOption = C % AllocatedDevice )
    else
      call C % Apply_EOS_HN_E_Kernel &
             ( P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               M, N, YE, UseDeviceOption = C % AllocatedDevice )
    end if
    call C % Compute_DS_G_Kernel &
           ( DS, N, SB, UseDeviceOption = C % AllocatedDevice )
    call T_AE % Stop ( )

!    if ( associated ( C % Value, Value_C ) ) &
!      call InterpolateSoundSpeed &
!             ( CS, P, C % Features, &
!               GV ( oV + 1 : oV + nV, G % CENTER_U ( 1 ) ) )

    call T_CE % Start ( )
    call C % Compute_FE_P_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
             V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
             UseDeviceOption = C % AllocatedDevice )
    call T_CE % Stop ( )
    
    end associate !-- T_CFP, etc.

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end select !-- FF
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
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
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

    integer ( KDI ) :: &
      iElectron, &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( T_CRF => PROGRAM_HEADER % Timer ( C % iTimerComputeRawFluxes ) )

    call C % ComputeRawFluxesTemplate_P &
           ( RawFlux, G, Storage_C, Storage_G, iDimension, nValuesOption, &
             oValueOption )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( Storage_C % Value, dim = 1 )
    end if

    call T_CRF % Start ( )

    call Search ( C % iaConserved, C % CONSERVED_ELECTRON_DENSITY, iElectron )

    associate &
      ( F_DE  => RawFlux % Value ( oV + 1 : oV + nV, iElectron ), &
        DE    => Storage_C % Value ( oV + 1 : oV + nV, &
                                     C % CONSERVED_ELECTRON_DENSITY ), &
        V_Dim => Storage_C % Value ( oV + 1 : oV + nV, &
                                     C % VELOCITY_U ( iDimension ) ) )

    call ComputeRawFluxesKernel &
           ( F_DE, DE, V_Dim, UseDeviceOption = C % AllocatedDevice )

    end associate !-- F_DE, etc.
    
    call T_CRF % Stop ( )
    
    end associate !-- T_CRF, etc. 

  end subroutine ComputeRawFluxes


  subroutine ComputeCenterStates &
               ( C_ICL, C_ICR, C, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iD )

    type ( StorageForm ), intent ( inout ) :: &
      C_ICL, C_ICR
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      C
    type ( StorageForm ), intent ( in ) :: &
      C_IL, C_IR, &
      SS_I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_DD_22, M_DD_33
    integer ( KDI ), intent ( in ) :: &
      iD

    associate &
      ( T_CCS => PROGRAM_HEADER % Timer ( C % iTimerComputeCenterStates ) )

    call C % ComputeCenterStatesTemplate_P &
           ( C_ICL, C_ICR, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iD )

    call T_CCS % Start ( )
    call ComputeCenterStatesKernel &
           ( C_ICL % Value ( :, C % CONSERVED_ELECTRON_DENSITY ), &
             C_ICR % Value ( :, C % CONSERVED_ELECTRON_DENSITY ), &
             C_IL % Value ( :, C % CONSERVED_ELECTRON_DENSITY ), &
             C_IR % Value ( :, C % CONSERVED_ELECTRON_DENSITY ), &
             C_IL % Value ( :, C % VELOCITY_U ( iD ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( iD ) ), &
             SS_I % Value ( :, C % ALPHA_PLUS ), &
             SS_I % Value ( :, C % ALPHA_MINUS ), &
             SS_I % Value ( :, C % ALPHA_CENTER ), &
             UseDeviceOption = C % AllocatedDevice )
    call T_CCS % Stop ( )
    
    end associate

  end subroutine ComputeCenterStates


  ! subroutine InterpolateSoundSpeed ( CS, P, FF, R )

  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     CS, &
  !     P
  !   class ( FluidFeaturesTemplate ), intent ( in ) :: &
  !     FF
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     R

  !   integer ( KDI ) :: &
  !     iC, &
  !     i_1, i_2
  !   real ( KDR ) :: &
  !     R_1, R_2, &
  !     CS_1, CS_2, &
  !     P_1, P_2

  !   select type ( Grid => FF % Grid )
  !   class is ( Chart_SL_Template )

  !   select case ( Grid % CoordinateSystem )
  !   case ( 'SPHERICAL' )

  !   i_1 = size ( R )
  !   do iC = 1, size ( R )
  !     if ( R ( iC )  >  8.0_KDR  *  UNIT % KILOMETER ) then
  !       i_1 = iC
  !       exit
  !     end if
  !   end do !-- iC

  !   i_2 = 1
  !   do iC = size ( R ), 1, -1
  !     if ( R ( iC )  <  13.0_KDR  *  UNIT % KILOMETER ) then
  !       i_2 = iC
  !       exit
  !     end if
  !   end do !-- iC

  !   R_1  =  R ( i_1 )
  !   R_2  =  R ( i_2 )

  !   CS_1  =  CS ( i_1 )
  !   CS_2  =  CS ( i_2 )

  !   !P_1  =  P ( i_1 )
  !   !P_2  =  P ( i_2 )

  !   do iC = i_1 + 1, i_2 - 1
  !     CS ( iC )  =  CS_1  +  ( CS_2 - CS_1 ) / ( R_2 - R_1 )  &
  !                            *  ( R ( iC )  -  R_1 )
  !     !P ( iC )   =  P_1  +  ( P_2 - P_1 ) / ( R_2 - R_1 )  &
  !     !                       *  ( R ( iC )  -  R_1 )
  !   end do !-- iC

  !   end select !-- CoordinateSystem
  !   end select !-- Grid

  ! end subroutine InterpolateSoundSpeed


end module Fluid_P_HN__Form
