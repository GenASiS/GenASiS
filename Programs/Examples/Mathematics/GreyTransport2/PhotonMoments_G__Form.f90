module PhotonMoments_G__Form

  !-- PhotonMoments_Grey__Form

  use Basics
  use Mathematics
  use RadiationMoments_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_PM = 0, &
      N_CONSERVED_PM = 0, &
      N_FIELDS_PM    = 2, &
      N_VECTORS_PM   = 0

  type, public, extends ( RadiationMomentsForm ) :: PhotonMoments_G_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_PM           = N_PRIMITIVE_PM, &
      N_CONSERVED_PM           = N_CONSERVED_PM, &
      N_FIELDS_PM              = N_FIELDS_PM, &
      N_VECTORS_PM             = N_VECTORS_PM, &
      TEMPERATURE_PARAMETER    = 0, &
      TEMPERATURE_PARAMETER_EQ = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_PM
    generic, public :: &
      Initialize => InitializeAllocate_PM
    procedure, public, pass :: &
      SetOutput
    final :: &
      Finalize
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, private, pass ( PM ) :: &
      ComputeSpectralParameters_PM
    generic, public :: &
      ComputeSpectralParameters => ComputeSpectralParameters_PM
  end type PhotonMoments_G_Form

    private :: &
      InitializeBasics, &
      SetUnits

contains


  subroutine InitializeAllocate_PM &
               ( PM, RiemannSolverType, UseLimiter, Velocity_U_Unit, &
                 MomentumDensity_U_Unit, MomentumDensity_D_Unit, &
                 EnergyDensityUnit, TemperatureUnit, LimiterParameter, &
                 nValues, VariableOption, VectorOption, NameOption, &
                 ClearOption, UnitOption, VectorIndicesOption )

    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      PM
    character ( * ), intent ( in ) :: &
      RiemannSolverType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
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
           ( PM, Variable, VariableUnit, VariableOption, UnitOption )

    call SetUnits ( VariableUnit, PM, TemperatureUnit )

    call PM % RadiationMomentsForm % Initialize &
           ( RiemannSolverType, UseLimiter, Velocity_U_Unit, &
             MomentumDensity_U_Unit, MomentumDensity_D_Unit, &
             EnergyDensityUnit, LimiterParameter, &
             nValues, VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeAllocate_PM


  subroutine SetOutput ( RM, Output )

    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      RM
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( RM % COMOVING_MOMENTUM_U )

    call Output % Initialize &
           ( RM, iaSelectedOption = [ RM % COMOVING_ENERGY, &
!                                       RM % COMOVING_NUMBER_DENSITY, &
                                      RM % COMOVING_MOMENTUM_U, &
                                      RM % FLUX_FACTOR, &
                                      RM % STRESS_FACTOR, &
                                      RM % TEMPERATURE_PARAMETER, &
                                      RM % TEMPERATURE_PARAMETER_EQ, &
!                                       RM % DEGENERACY_PARAMETER, &
!                                       RM % DEGENERACY_PARAMETER_EQ, &
!                                       RM % ENERGY_AVERAGE, &
!                                       RM % OCCUPANCY_AVERAGE, &
                                      RM % COMOVING_ENERGY_EQ ], &
             VectorOption = [ 'ComovingMomentum_U' ], &
!                               'ComovingNumberFlux             ' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  impure elemental subroutine Finalize ( PM )

    type ( PhotonMoments_G_Form ), intent ( inout ) :: &
      PM

    !-- Trigger finalization in parent

  end subroutine Finalize


  subroutine ComputeFromPrimitiveCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( PhotonMoments_G_Form ), intent ( in ) :: &
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
      ( J    => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY ), &
        J_EQ => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY_EQ ), &
        T    => RMV ( oV + 1 : oV + nV, C % TEMPERATURE_PARAMETER ), &
        T_EQ => RMV ( oV + 1 : oV + nV, C % TEMPERATURE_PARAMETER_EQ ) )

    if ( associated ( C % Interactions ) ) &
      call C % Interactions % ComputeEquilibriumParameters ( T_EQ, C )
    call C % ComputeSpectralParameters ( T, J_EQ, J, T_EQ )

    end associate !-- J, etc.
    nullify ( RMV )

  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( PhotonMoments_G_Form ), intent ( in ) :: &
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
      ( J    => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY ), &
        J_EQ => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY_EQ ), &
        T    => RMV ( oV + 1 : oV + nV, C % TEMPERATURE_PARAMETER ), &
        T_EQ => RMV ( oV + 1 : oV + nV, C % TEMPERATURE_PARAMETER_EQ ) )

    if ( associated ( C % Interactions ) ) &
      call C % Interactions % ComputeEquilibriumParameters ( T_EQ, C )
    call C % ComputeSpectralParameters ( T, J_EQ, J, T_EQ )

    end associate !-- J, etc.
    nullify ( RMV )

  end subroutine ComputeFromConservedCommon


  subroutine ComputeSpectralParameters_PM &
!               ( T, Eta, E_Ave, F_Ave, J_EQ, RM, J, N, T_EQ, Eta_EQ )
               ( T, J_EQ, PM, J, T_EQ )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
!       Eta, &
!       E_Ave, &
!       F_Ave, &
      J_EQ
    class ( PhotonMoments_G_Form ), intent ( in ) :: &
      PM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
!       N, &
      T_EQ!, &
!       Eta_EQ

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      a

    nValues = size ( T )

    a = 4.0_KDR * CONSTANT % STEFAN_BOLTZMANN

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      T    ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
      J_EQ ( iV )  =  a  *  T_EQ ( iV ) ** 4
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeSpectralParameters_PM


  subroutine InitializeBasics &
               ( PM, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      PM
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

    if ( PM % Type == '' ) &
      PM % Type = 'PhotonMoments_G'

    !-- variable indices

    oF = PM % N_FIELDS_TEMPLATE + PM % N_FIELDS_RM
    if ( PM % N_FIELDS == 0 ) &
      PM % N_FIELDS = oF + PM % N_FIELDS_PM

    PM % TEMPERATURE_PARAMETER     =  oF + 1
    PM % TEMPERATURE_PARAMETER_EQ  =  oF + 2

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( PM % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + PM % N_FIELDS_PM ) &
      = [ 'TemperatureParameter   ', &
          'TemperatureParameter_EQ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( PM % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = PM % N_VECTORS_TEMPLATE + PM % N_VECTORS_RM
    if ( PM % N_VECTORS == 0 ) &
      PM % N_VECTORS = oV + PM % N_VECTORS_PM

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, PM, TemperatureUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( PhotonMoments_G_Form ), intent ( in ) :: &
      PM
    type ( MeasuredValueForm ), intent ( in ) :: &
      TemperatureUnit

    VariableUnit ( PM % TEMPERATURE_PARAMETER )     =  TemperatureUnit
    VariableUnit ( PM % TEMPERATURE_PARAMETER_EQ )  =  TemperatureUnit

  end subroutine SetUnits


end module PhotonMoments_G__Form
