module Fluid_P_NR__Form

  !-- Fluid_Perfect_NonRelativistic__Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_NON_RELATIVISTIC = 0, &
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
      SetAdiabaticIndex
    procedure, public, pass :: &
      SetMeanMolecularWeight
    procedure, public, pass :: &
      SetFiducialBaryonDensity
    procedure, public, pass :: &
      SetFiducialTemperature
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, nopass :: &
      Apply_EOS_NR_Kernel
  end type Fluid_P_NR_Form

    private :: &
      InitializeBasics

contains


  subroutine InitializeAllocate_P_NR &
               ( F, VelocityUnit, MassDensityUnit, EnergyDensityUnit, &
                 TemperatureUnit, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( Fluid_P_NR_Form ), intent ( inout ) :: &
      F
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
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
           ( VelocityUnit, MassDensityUnit, EnergyDensityUnit, &
             TemperatureUnit, nValues, VariableOption = Variable, &
             VectorOption = VectorOption, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeAllocate_P_NR


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


  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Value_C, Value_G, iDimension, &
                 nValuesOption, oValueOption )
    
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    class ( Fluid_P_NR_Form ), intent ( in ) :: &
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

    call C % ComputeRawFluxesTemplate_P &
           ( RawFlux, G, Value_C, Value_G, iDimension, nValuesOption, &
             oValueOption )

  end subroutine ComputeRawFluxes


  subroutine SetOutput ( F, Output )

    class ( Fluid_P_NR_Form ), intent ( in ) :: &
      F
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( F % VELOCITY_U )
    call Output % Initialize &
           ( F, iaSelectedOption &
                  = [ F % COMOVING_DENSITY, F % VELOCITY_U, F % PRESSURE, &
                      F % SOUND_SPEED, F % MACH_NUMBER, &
                      F % ENTROPY_PER_BARYON ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  subroutine ComputeFromPrimitiveCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( Fluid_P_NR_Form ), intent ( in ) :: &
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
      
    associate &
      ( FV => Value_C, &
        GV => Value_G )

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

    call C % ComputeBaryonMassKernel ( M )
    call C % ComputeDensityMomentumKernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % ComputeConservedEnergyKernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
    call C % Apply_EOS_NR_Kernel &
           ( P, Gamma, SB, T, M, N, E, C % AdiabaticIndex, &
             C % MeanMolecularWeight, C % FiducialBaryonDensity, &
             C % FiducialTemperature, UNIT % ATOMIC_MASS_UNIT, &
             UNIT % BOLTZMANN )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( Fluid_P_NR_Form ), intent ( in ) :: &
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

    associate &
      ( FV => Value_C, &
        GV => Value_G )

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

    call C % ComputeBaryonMassKernel ( M )
    call C % ComputeDensityVelocityKernel &
           ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33 )
    call C % ComputeInternalEnergyKernel &
           ( E, G, M, N, V_1, V_2, V_3, S_1, S_2, S_3 )
    call C % Apply_EOS_NR_Kernel &
           ( P, Gamma, SB, T, M, N, E, C % AdiabaticIndex, &
             C % MeanMolecularWeight, C % FiducialBaryonDensity, &
             C % FiducialTemperature, UNIT % ATOMIC_MASS_UNIT, &
             UNIT % BOLTZMANN )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromConservedCommon


  subroutine Apply_EOS_NR_Kernel &
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
      T_0
    type ( MeasuredValueForm ), intent ( in ) :: &
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
        SB ( iV )  =  - huge ( 1.0_KDR )
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine Apply_EOS_NR_Kernel


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

    !-- select primitive, conserved

    oP = F % N_PRIMITIVE_TEMPLATE + F % N_PRIMITIVE_DUST &
         + F % N_PRIMITIVE_PERFECT
    oC = F % N_CONSERVED_TEMPLATE + F % N_CONSERVED_DUST &
         + F % N_CONSERVED_PERFECT

    if ( .not. allocated ( F % iaPrimitive ) ) then
      F % N_PRIMITIVE = oP + F % N_PRIMITIVE_NON_RELATIVISTIC
      allocate ( F % iaPrimitive ( F % N_PRIMITIVE ) )
    end if
!    F % iaPrimitive ( oP + 1 : oP + F % N_PRIMITIVE_NON_RELATIVISTIC ) &
!      = [ ]

    if ( .not. allocated ( F % iaConserved ) ) then
      F % N_CONSERVED = oC + F % N_CONSERVED_NON_RELATIVISTIC
      allocate ( F % iaConserved ( F % N_CONSERVED ) )
    end if
!    F % iaConserved ( oC + 1 : oC + F % N_CONSERVED_NON_RELATIVISTIC ) &
!      = [ ]
    
  end subroutine InitializeBasics
  
    
end module Fluid_P_NR__Form
