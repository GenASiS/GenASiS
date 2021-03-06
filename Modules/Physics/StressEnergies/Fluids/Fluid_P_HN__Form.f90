module Fluid_P_HN__Form

  !-- Fluid_Perfect_RepresentativeHeavyNucleus__Form
  
  use Basics
  use Mathematics
  use StressEnergyBasics
!  use FluidFeatures_Template
  use EOS_P_HN_OConnorOtt__Form
  use Fluid_P__Template
  use FluidFeatures_P__Form
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_HEAVY_NUCLEUS = 1, &
!      N_PRIMITIVE_HEAVY_NUCLEUS = 2, &
      N_CONSERVED_HEAVY_NUCLEUS = 1, &
      N_FIELDS_HEAVY_NUCLEUS    = 11, &
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
      CHEMICAL_POTENTIAL_E        = 0, &
      UNUSED_VARIABLE             = 0
        !-- Includes m_e.
    logical ( KDL ), private :: &
      Allocated_EOS = .false.
    type ( EOS_P_HN_OConnorOtt_Form ), public, pointer :: &
      EOS => null ( )
  contains
    procedure, public, pass :: &
      Initialize_P_HN
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_P_HN
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetOutput
    final :: &
      Finalize
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
    procedure, private, nopass :: &
      Apply_EOS_PrologueKernel
    procedure, private, nopass :: &
      Apply_EOS_EpilogueKernel
!    procedure, public, nopass :: &
!      Apply_EOS_HN_SB_Kernel
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
    
    !-- OConnorOtt NucEOS-specific variables
    real ( KDR ), public, protected :: &
      EOS_RF_Accuracy     !-- EOS_RootFinding_Accuracy
    integer ( KDI ), public, parameter :: &
      EOS_Apply_EOS_HN_T = 1_KDI, &   !-- T input
      EOS_Apply_EOS_HN_E = 0_KDI, &   !-- E input, solve for T
      EOS_Apply_EOS_HN_S = 2_KDI      !-- S input, solve for T
    logical ( KDL ), private, protected :: &
      EOS_Initialized = .false.
    type ( EOS_P_HN_OConnorOtt_Form ), pointer, private, protected :: &
      EOS_Pointer
      
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
             ( N, T, P, E, CS, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               U_V, T_EOS, M, YE, T_L_D, T_L_T, T_YE, Error, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        N, &
        T, &
        P, &
        E, &
        CS, &
        SB, &
        X_P, X_N, X_He, X_A, &
        Z, A, &
        Mu_NP, Mu_E, &
        U_V             !-- Dummy storage for Unused Variable
      real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
        T_EOS
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        YE, &
        T_L_D, &  !-- TableLogDensity
        T_L_T, &  !-- TableLogTemperature
        T_YE      !-- TableElectronFraction
      integer ( KDI ), dimension ( : ), intent ( out ) :: &
        Error
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_HN_T_Kernel
    

    module subroutine Apply_EOS_HN_SB_E_Kernel &
             ( N, P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               U_V, Error_A, T_EOS, M, YE, Shock, T_L_D, T_L_T, T_YE, & 
               Error, UseDeviceOption ) 
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        N, &
        P, &
        T, &
        CS, &
        E, &
        SB, &
        X_P, X_N, X_He, X_A, &
        Z, A, &
        Mu_NP, Mu_E, &
        U_V, &          !-- Dummy storage for Unused Variable
        Error_A 
      real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
        T_EOS
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        YE, &
        Shock, &
        T_L_D, &  !-- TableLogDensity
        T_L_T, &  !-- TableLogTemperature
        T_YE      !-- TableElectronFraction
      integer ( KDI ), dimension ( : ), intent ( out ) :: &
        Error
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_HN_SB_E_Kernel
    
    
    module subroutine Apply_EOS_HN_E_Kernel &
             ( N, P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
               U_V, Error_A, T_EOS, M, YE, T_L_D, T_L_T, T_YE, Error, &
               UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        N, &
        P, &
        T, &
        CS, &
        E, &
        SB, &
        X_P, X_N, X_He, X_A, &
        Z, A, &
        Mu_NP, Mu_E, &
        U_V, &
        Error_A
      real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
        T_EOS
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        YE, &
        T_L_D, &  !-- TableLogDensity
        T_L_T, &  !-- TableLogTemperature
        T_YE      !-- TableElectronFraction
      integer ( KDI ), dimension ( : ), intent ( out ) :: &
        Error
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_HN_E_Kernel


    module subroutine Apply_EOS_PrologueKernel &
             ( N, P, T, E, M, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        N, &
        P, &
        T, &
        E
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_PrologueKernel
    
    
    module subroutine Apply_EOS_EpilogueKernel &
             ( N, P, T, CS, E, Mu_NP, Mu_E, M, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        N, &
        P, &
        T, &
        CS, &
        E, &
        Mu_NP, &
        Mu_E
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_EpilogueKernel


!    module subroutine Apply_EOS_HN_SB_Kernel &
!             ( N, P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
!               U_V, T_EOS, M, N, YE, T_L_D, T_L_T, T_YE, Error, &
!               UseDeviceOption )
!      use Basics
!      real ( KDR ), dimension ( : ), intent ( inout ) :: &
!        N, &
!        P, &
!        T, &
!        CS, &
!        E, &
!        SB, &
!        X_P, X_N, X_He, X_A, &
!        Z, A, &
!        Mu_NP, Mu_E, &
!        U_V
!      real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
!        T_EOS
!      real ( KDR ), dimension ( : ), intent ( in ) :: &
!        M, &
!        YE, &
!        T_L_D, &  !-- TableLogDensity
!        T_L_T, &  !-- TableLogTemperature
!        T_YE      !-- TableElectronFraction
!      integer ( KDI ), dimension ( : ), intent ( out ) :: &
!        Error
!      logical ( KDL ), intent ( in ), optional :: &
!        UseDeviceOption
!    end subroutine Apply_EOS_HN_SB_Kernel


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
                 UseEntropy, UseInitialTemperature, UseLimiter, Units, &
                 BaryonMassReference, LimiterParameter, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F
    character ( * ), intent ( in ) :: &
      FluidType, &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseEntropy, &
      UseInitialTemperature, &
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

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaFluidOutput, &
      iaSelected_EOS
    character ( LDF ) :: &
      EOS_Filename
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit

    call InitializeBasics &
           ( F, Variable, VariableUnit, VariableOption, UnitOption )

    call SetUnits ( VariableUnit, F, Units )

    call F % InitializeTemplate_P &
           ( FluidType, RiemannSolverType, ReconstructedType, UseEntropy, &
             UseInitialTemperature, UseLimiter, Units, BaryonMassReference, &
             LimiterParameter, nValues, VariableOption = Variable, &
             VectorOption = VectorOption, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )
    
    if ( EOS_Initialized ) then
      F % EOS => EOS_Pointer
      return
    end if

    allocate ( F % EOS )
    
    call DelayFileAccess ( PROGRAM_HEADER % Communicator % Rank )
    
    EOS_Filename = '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5'
    call PROGRAM_HEADER % GetParameter ( EOS_Filename, 'EOS_Filename' )
    call F % EOS % Initialize ( EOS_Filename )
    
    allocate ( iaFluidOutput ( 12 ) )
    allocate ( iaSelected_EOS ( 12 ) )
    
    iaFluidOutput &
      = [ F % INTERNAL_ENERGY, &
          F % PRESSURE, &
          F % ENTROPY_PER_BARYON, &
          F % SOUND_SPEED, &
          F % MASS_FRACTION_ALPHA, &
          F % MASS_FRACTION_HEAVY, &
          F % MASS_FRACTION_NEUTRON, &
          F % MASS_FRACTION_PROTON, &
          F % MASS_NUMBER_HEAVY, &
          F % ATOMIC_NUMBER_HEAVY, &
          F % CHEMICAL_POTENTIAL_E, &
          F % CHEMICAL_POTENTIAL_N_P ]

    iaSelected_EOS &
      = [ F % EOS % LOG_ENERGY, &
          F % EOS % LOG_PRESSURE, &
          F % EOS % ENTROPY, &
          F % EOS % SOUND_SPEED_SQUARE, &
          F % EOS % MASS_FRACTION_A, &
          F % EOS % MASS_FRACTION_H, &
          F % EOS % MASS_FRACTION_N, &
          F % EOS % MASS_FRACTION_P, &
          F % EOS % MASS_NUMBER_BAR, &
          F % EOS % ATOMIC_NUMBER_BAR, &
          F % EOS % CHEMICAL_POTENTIAL_E, &
          F % EOS % CHEMICAL_POTENTIAL_HAT ]
    
    call F % EOS % SelectVariables ( iaFluidOutput, iaSelected_EOS )

    F % Allocated_EOS = .true.
    EOS_Pointer => F % EOS

    !-- Historical Oak Ridge Shift, accounting for nuclear binding energy
    OR_Shift = 8.9_KDR * UNIT % MEGA_ELECTRON_VOLT &
               / CONSTANT % ATOMIC_MASS_UNIT
    
    MassDensity_CGS     =  UNIT % MASS_DENSITY_CGS
    SpecificEnergy_CGS  =  UNIT % ERG  /  UNIT % GRAM
    Pressure_CGS        =  UNIT % BARYE
    Speed_CGS           =  UNIT % CENTIMETER  /  UNIT % SECOND
    MeV                 =  UNIT % MEGA_ELECTRON_VOLT
    EOS_RF_Accuracy     =  1.0e-9_KDR

    EOS_Initialized  =  .true.

  end subroutine Initialize_P_HN
  
  
  subroutine AllocateDevice_P_HN ( S, AssociateVariablesOption ) 
    
    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      S
    logical ( KDL ), intent ( in ), optional :: &
      AssociateVariablesOption
      
    type ( c_ptr ) :: &
      D_P  !-- Device Pointer
   
    call S % EOS % AllocateDevice ( )
    
    call S % StorageForm % AllocateDevice ( AssociateVariablesOption )
      
  end subroutine AllocateDevice_P_HN


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
!      = [ C % TEMPERATURE, C % ELECTRON_FRACTION ]

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
    F % UNUSED_VARIABLE             =  of + 11

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
          'ChemicalPotential_E     ', &
          'UnusedVariable          ' ]
    
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
  
  
  impure elemental subroutine Finalize ( F )
    
    type ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F
    
    if ( F % Allocated_EOS ) then
      if ( associated ( EOS_Pointer, F % EOS ) ) &
        nullify ( EOS_Pointer )
      deallocate ( F % EOS )
      !-- FIXME: Need to deallocate device memory for EOS
    else
      nullify ( F % EOS )
    end if
    
    nullify ( F % EOS )
    
  end subroutine Finalize 
  
  
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
        Mu_E  => FV ( oV + 1 : oV + nV, C % CHEMICAL_POTENTIAL_E ), &
        U_V   => FV ( oV + 1 : oV + nV, C % UNUSED_VARIABLE ) )

    call C % Compute_M_Kernel &
           ( M, C % BaryonMassReference, &
             UseDeviceOption = C % AllocatedDevice )
    !call C % Apply_EOS_HN_T_Kernel &
    !       ( N, T, P, E, CS, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, U_V, &
    !         T_EOS, M, YE, T_L_D, T_L_T, T_YE, Error, &
    !         UseDeviceOption = C % AllocatedDevice )
    call C % Apply_EOS_PrologueKernel &
           ( N, P, T, E, M, UseDeviceOption = C % AllocatedDevice )
    
    call Storage_C % ReassociateHost ( AssociateVariablesOption = .false. )
    call C % EOS % ComputeFromTemperature &
           ( Storage_C, &
             iaFluidInput = [ C % COMOVING_BARYON_DENSITY, &
                              C % TEMPERATURE, C % ELECTRON_FRACTION ] )
    call Storage_C % ReassociateHost ( AssociateVariablesOption = .true. )
    
    call C % Apply_EOS_EpilogueKernel &
           ( N, P, T, CS, E, Mu_NP, Mu_E, M, &
             UseDeviceOption = C % AllocatedDevice )
        
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
      
!    call ComputeFromTemperature &
!           ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

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
        Mu_E  => FV ( oV + 1 : oV + nV, C % CHEMICAL_POTENTIAL_E ), &
        U_V   => FV ( oV + 1 : oV + nV, C % UNUSED_VARIABLE ), &
        Error_A => FF % Value ( oV + 1 : oV + nV, FF % EOS_ERROR ) )

    call Copy ( C % Value ( :, C % PRESSURE ), P, &
                UseDeviceOption = C % AllocatedDevice )
    call Copy ( C % Value ( :, C % TEMPERATURE ), T, &
                UseDeviceOption = C % AllocatedDevice )

    call C % Compute_M_Kernel &
           ( M, C % BaryonMassReference, &
             UseDeviceOption = C % AllocatedDevice )

!    call C % Apply_EOS_HN_T_Kernel &
!           ( P, E, CS, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
!             M, N, T, YE )
!    call C % Apply_EOS_HN_E_Kernel &
!           ( N, P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
!             U_V, Error_A, T_EOS, M, YE, T_L_D, T_L_T, T_YE, Error, &
!             UseDeviceOption = C % AllocatedDevice )
!    call C % Apply_EOS_HN_SB_Kernel &
!           ( P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
!             M, N, YE )

!call Show ( '>>> 3.1' )
    call C % Apply_EOS_PrologueKernel &
           ( N, P, T, E, M, UseDeviceOption = C % AllocatedDevice )
    
!call Show ( '>>> 3.2' )
    call Storage_C % ReassociateHost ( AssociateVariablesOption = .false. )
    call C % EOS % ComputeFromEnergy &
           ( Storage_C, &
             iaFluidInput = [ C % COMOVING_BARYON_DENSITY, &
                              C % TEMPERATURE, C % ELECTRON_FRACTION ], &
             iSolve = C % INTERNAL_ENERGY )
    call Storage_C % ReassociateHost ( AssociateVariablesOption = .true. )
    
!call Show ( '>>> 3.3' )
    call C % Apply_EOS_EpilogueKernel &
           ( N, P, T, CS, E, Mu_NP, Mu_E, M, &
             UseDeviceOption = C % AllocatedDevice )
    
!call Show ( '>>> 3.4' )
    call C % Compute_D_S_G_Kernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33, &
             UseDeviceOption = C % AllocatedDevice )
!call Show ( '>>> 3.5' )
    call C % Compute_G_G_Kernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E, &
             UseDeviceOption = C % AllocatedDevice )
!call Show ( '>>> 3.6' )
    call C % Compute_DS_G_Kernel &
           ( DS, N, SB, UseDeviceOption = C % AllocatedDevice )
!call Show ( '>>> 3.7' )
    call C % Compute_DE_G_Kernel &
           ( DE, N, YE, UseDeviceOption = C % AllocatedDevice )

!call Show ( '>>> 3.8' )
    call C % Compute_FE_P_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
             V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
             UseDeviceOption = C % AllocatedDevice )
    
!call Show ( '>>> 3.9' )
    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end select !-- FF
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
        U_V   => FV ( oV + 1 : oV + nV, C % UNUSED_VARIABLE ), &
        Shock => FF % Value ( oV + 1 : oV + nV, FF % SHOCK ), &
        Error_A => FF % Value ( oV + 1 : oV + nV, FF % EOS_ERROR ) )
    
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
    
    if ( C % UseEntropy ) then
      call C % Compute_SB_G_Kernel &
             ( SB, DS, N, UseDeviceOption = C % AllocatedDevice )
      !call C % Apply_EOS_HN_SB_E_Kernel &
      !       ( N, P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
      !         U_V, Error_A, T_EOS, M, YE, Shock, T_L_D, T_L_T, T_YE, Error, &
      !         UseDeviceOption = C % AllocatedDevice )
      
      call C % Apply_EOS_PrologueKernel &
           ( N, P, T, E, M, UseDeviceOption = C % AllocatedDevice )
      
      call Storage_C % ReassociateHost ( AssociateVariablesOption = .false. )
      call C % EOS % ComputeFromEnergyEntropy &
             ( Storage_C, Mask = Shock, Threshold = 0.0_KDR, &
               iaFluidInput = [ C % COMOVING_BARYON_DENSITY, &
                                C % TEMPERATURE, C % ELECTRON_FRACTION ], &
               iEnergy = C % INTERNAL_ENERGY, &
               iEntropy = C % ENTROPY_PER_BARYON )
      call Storage_C % ReassociateHost ( AssociateVariablesOption = .true. )
      
      call C % Apply_EOS_EpilogueKernel &
             ( N, P, T, CS, E, Mu_NP, Mu_E, M, &
               UseDeviceOption = C % AllocatedDevice )
      
      call C % Compute_G_G_Kernel &
             ( GE, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E, &
               UseDeviceOption = C % AllocatedDevice )
    else
!      call C % Apply_EOS_HN_E_Kernel &
!             ( N, P, T, CS, E, SB, X_P, X_N, X_He, X_A, Z, A, Mu_NP, Mu_E, &
!               U_V, Error_A, T_EOS, M, YE, T_L_D, T_L_T, T_YE, Error, &
!               UseDeviceOption = C % AllocatedDevice )
      call C % Apply_EOS_PrologueKernel &
           ( N, P, T, E, M, UseDeviceOption = C % AllocatedDevice )
    
      call Storage_C % ReassociateHost ( AssociateVariablesOption = .false. )
      call C % EOS % ComputeFromEnergy &
             ( Storage_C, &
               iaFluidInput = [ C % COMOVING_BARYON_DENSITY, &
                                C % TEMPERATURE, C % ELECTRON_FRACTION ], &
               iSolve = C % INTERNAL_ENERGY )
      call Storage_C % ReassociateHost ( AssociateVariablesOption = .true. )
      
      call C % Apply_EOS_EpilogueKernel &
             ( N, P, T, CS, E, Mu_NP, Mu_E, M, &
               UseDeviceOption = C % AllocatedDevice )

    end if
    call C % Compute_DS_G_Kernel &
           ( DS, N, SB, UseDeviceOption = C % AllocatedDevice )

!    if ( associated ( C % Value, Value_C ) ) &
!      call InterpolateSoundSpeed &
!             ( CS, P, C % Features, &
!               GV ( oV + 1 : oV + nV, G % CENTER_U ( 1 ) ) )

    call C % Compute_FE_P_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
             V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
             UseDeviceOption = C % AllocatedDevice )
    
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

    call C % ComputeCenterStatesTemplate_P &
           ( C_ICL, C_ICR, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iD )

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
