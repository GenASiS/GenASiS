module Fluid_P_HN__Form

  !-- Fluid_Perfect_RepresentativeHeavyNucleus__Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_HEAVY_NUCLEUS = 2, &
      N_CONSERVED_HEAVY_NUCLEUS = 2, &
      N_FIELDS_HEAVY_NUCLEUS    = 11, &
      N_VECTORS_HEAVY_NUCLEUS   = 0

  type, public, extends ( Fluid_P_Template ) :: Fluid_P_HN_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_HEAVY_NUCLEUS = N_PRIMITIVE_HEAVY_NUCLEUS, &
      N_CONSERVED_HEAVY_NUCLEUS = N_CONSERVED_HEAVY_NUCLEUS, &
      N_FIELDS_HEAVY_NUCLEUS    = N_FIELDS_HEAVY_NUCLEUS, &
      N_VECTORS_HEAVY_NUCLEUS   = N_VECTORS_HEAVY_NUCLEUS, &
      PROTON_FRACTION           = 0, &
      CONSERVED_PROTON_DENSITY  = 0, &
      MASS_FRACTION_PROTON      = 0, &
      MASS_FRACTION_NEUTRON     = 0, &
      MASS_FRACTION_ALPHA       = 0, &
      MASS_FRACTION_HEAVY       = 0, &
      ATOMIC_NUMBER_HEAVY       = 0, &
      MASS_NUMBER_HEAVY         = 0, &
      CHEMICAL_POTENTIAL_N_P    = 0, &  
        !-- a.k.a. mu_hat. Includes m_n - m_p. (mu_n and mu_p both
        !   measured with respect to m_n.)
      CHEMICAL_POTENTIAL_E      = 0     
        !-- Includes m_e.
  contains
    procedure, public, pass :: &
      InitializeAllocate_P_HN
    generic, public :: &
      Initialize => InitializeAllocate_P_HN
  !   procedure, public, pass :: &
  !     SetPrimitiveConserved
  !   procedure, public, pass :: &
  !     SetOutput
  !   procedure, public, pass ( C ) :: &
  !     ComputeFromPrimitiveCommon
  !   procedure, public, pass ( C ) :: &
  !     ComputeFromConservedCommon
  !   procedure, public, pass ( C ) :: &
  !     ComputeRawFluxes
  !   procedure, public, pass ( C ) :: &
  !     ComputeCenterStates
  !   procedure, public, nopass :: &
  !     ComputeConservedProtonKernel
  !   procedure, public, nopass :: &
  !     ComputeConservedEntropyKernel
  !   procedure, public, nopass :: &
  !     ComputeElectronFractionKernel
  !   procedure, public, nopass :: &
  !     ComputeEntropyPerBaryonKernel
  !   procedure, public, nopass :: &
  !     Apply_EOS_HN_T_Kernel
  !   procedure, public, nopass :: &
  !     Apply_EOS_HN_SB_E_Kernel
  end type Fluid_P_HN_Form

    private :: &
      InitializeBasics, &
      SetUnits

    real ( KDR ), private :: &
      OR_Shift, &
      MassDensity_CGS, &
      SpecificEnergy_CGS, &
      Pressure_CGS, &
      MeV

contains


  subroutine InitializeAllocate_P_HN &
               ( F, RiemannSolverType, UseLimiter, Velocity_U_Unit, &
                 MomentumDensity_D_Unit, BaryonMassUnit, NumberDensityUnit, &
                 EnergyDensityUnit, TemperatureUnit, BaryonMassReference, &
                 LimiterParameter, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F
    character ( * ), intent ( in ) :: &
      RiemannSolverType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit, &
      MomentumDensity_D_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
      BaryonMassUnit, &
      NumberDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
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

    call SetUnits ( VariableUnit, F, NumberDensityUnit, TemperatureUnit )

    call F % InitializeTemplate_P &
           ( RiemannSolverType, UseLimiter, Velocity_U_Unit, &
             MomentumDensity_D_Unit, BaryonMassUnit, NumberDensityUnit, &
             EnergyDensityUnit, TemperatureUnit, BaryonMassReference, &
             LimiterParameter, nValues, VariableOption = Variable, &
             VectorOption = VectorOption, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    call READTABLE &
           ( '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5' )
!    call READTABLE &
!           ( '../Parameters/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5' )

    !-- Historical Oak Ridge Shift, accounting for nuclear binding energy
    OR_Shift = 8.9_KDR * UNIT % MEGA_ELECTRON_VOLT &
               / CONSTANT % ATOMIC_MASS_UNIT
    
    MassDensity_CGS     =  UNIT % MASS_DENSITY_CGS
    SpecificEnergy_CGS  =  UNIT % ERG  /  UNIT % GRAM
    Pressure_CGS        =  UNIT % BARYE
    MeV                 =  UNIT % MEGA_ELECTRON_VOLT

  end subroutine InitializeAllocate_P_HN


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

    F % PROTON_FRACTION          =  oF +  1
    F % CONSERVED_PROTON_DENSITY =  oF +  2
    F % MASS_FRACTION_PROTON     =  oF +  3
    F % MASS_FRACTION_NEUTRON    =  oF +  4
    F % MASS_FRACTION_ALPHA      =  oF +  5
    F % MASS_FRACTION_HEAVY      =  oF +  6
    F % ATOMIC_NUMBER_HEAVY      =  oF +  7
    F % MASS_NUMBER_HEAVY        =  oF +  8
    F % CHEMICAL_POTENTIAL_N_P   =  oF +  9
    F % CHEMICAL_POTENTIAL_E     =  oF + 10

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_HEAVY_NUCLEUS ) &
      = [ 'ProtonFraction        ', &
          'ConservedProtonDensity', &
          'MassFractionProton    ', &
          'MassFractionNeutron   ', &
          'MassFractionAlpha     ', &
          'MassFractionHeavy     ', &
          'AtomicNumberHeavy     ', &
          'MassNumberHeavy       ', &
          'ChemicalPotential_N_P ', &
          'ChemicalPotential_E   ']
    
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


  subroutine SetUnits ( VariableUnit, F, NumberDensityUnit, TemperatureUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), intent ( in ) :: &
      NumberDensityUnit, &
      TemperatureUnit

    VariableUnit ( F % CONSERVED_PROTON_DENSITY ) &
      = NumberDensityUnit
    VariableUnit ( F % CHEMICAL_POTENTIAL_N_P ) &
      = TemperatureUnit
    VariableUnit ( F % CHEMICAL_POTENTIAL_E ) &
      = TemperatureUnit

  end subroutine SetUnits


end module Fluid_P_HN__Form
