module Fluid_P_I__Form

  !-- Fluid_Perfect_Ideal__Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_IDEAL = 1, &
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
      AtomicMassUnit, &
      BoltzmannConstant, &
      AdiabaticIndex, &
      MeanMolecularWeight, &
      SpecificHeatVolume, &
      FiducialBaryonDensity, &
      FiducialPressure
  contains
    procedure, public, pass :: &
      InitializeAllocate_P_I
    generic, public :: &
      Initialize => InitializeAllocate_P_I
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
  end type Fluid_P_I_Form

    private :: &
      InitializeBasics


contains


  subroutine InitializeAllocate_P_I &
               ( F, RiemannSolverType, UseLimiter, Velocity_U_Unit, &
                 MomentumDensity_D_Unit, BaryonMassUnit, NumberDensityUnit, &
                 EnergyDensityUnit, TemperatureUnit, BaryonMassReference, &
                 LimiterParameter, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
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

    associate &
      ( amu   => F % AtomicMassUnit, &
        k     => F % BoltzmannConstant, &
        gamma => F % AdiabaticIndex, &
        mu    => F % MeanMolecularWeight, &
        c_v   => F % SpecificHeatVolume, &
        n_0   => F % FiducialBaryonDensity, &
        p_0   => F % FiducialPressure )

    if ( BaryonMassUnit /= UNIT % IDENTITY ) then
      amu = CONSTANT % ATOMIC_MASS_UNIT
    else !-- Dimensionless
      amu = 1.0_KDR
    end if

    if ( TemperatureUnit /= UNIT % IDENTITY ) then
      k = CONSTANT % BOLTZMANN
    else !-- Dimensionless
      k = 1.0_KDR
    end if

    gamma  =  1.4_KDR
    mu     =  1.0_KDR
    c_v    =  k / ( mu * ( gamma - 1.0_KDR ) )

    n_0  =  1.0_KDR
    p_0  =  1.0_KDR
    
    end associate !-- amu, etc.

    call InitializeBasics &
           ( F, Variable, VariableUnit, VariableOption, UnitOption )

    call F % InitializeTemplate_P &
           ( RiemannSolverType, UseLimiter, Velocity_U_Unit, &
             MomentumDensity_D_Unit, BaryonMassUnit, NumberDensityUnit, &
             EnergyDensityUnit, TemperatureUnit, BaryonMassReference, &
             LimiterParameter, nValues, VariableOption = Variable, &
             VectorOption = VectorOption, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeAllocate_P_I


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
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_IDEAL ) &
      = [ C % INTERNAL_ENERGY ]

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

  end subroutine SetFiducialParameters


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
