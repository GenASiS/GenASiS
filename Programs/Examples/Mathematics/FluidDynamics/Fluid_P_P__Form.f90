module Fluid_P_P__Form

  !-- Fluid_Perfect_Polytropic__Form

#include "Preprocessor"

  use iso_c_binding
  use Basics
  use Mathematics
  use Fluid_P__Template
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_POLYTROPIC = 1, &
      N_CONSERVED_POLYTROPIC = 0, &
      N_FIELDS_POLYTROPIC    = 1, &
      N_VECTORS_POLYTROPIC   = 0

  type, public, extends ( Fluid_P_Template ) :: Fluid_P_P_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_POLYTROPIC = N_PRIMITIVE_POLYTROPIC, &
      N_CONSERVED_POLYTROPIC = N_CONSERVED_POLYTROPIC, &
      N_FIELDS_POLYTROPIC    = N_FIELDS_POLYTROPIC, &
      N_VECTORS_POLYTROPIC   = N_VECTORS_POLYTROPIC, &
      POLYTROPIC_PARAMETER   = 0
    real ( KDR ) :: &
      AdiabaticIndex = 1.4_KDR, &
      FiducialPolytropicParameter = 1.0_KDR
  contains
    procedure, public, pass :: &
      InitializeAllocate_P_P
    generic, public :: &
      Initialize => InitializeAllocate_P_P
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetAdiabaticIndex
    procedure, public, pass :: &
      SetFiducialPolytropicParameter
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, nopass :: &
      Apply_EOS_P_KernelHost, &
      Apply_EOS_P_KernelDevice
  end type Fluid_P_P_Form

    private :: &
      InitializeBasics, &
      SetUnits

contains


  subroutine InitializeAllocate_P_P &
               ( F, RiemannSolverType, ReconstructedType, UseLimiter, &
                 VelocityUnit, MassDensityUnit, EnergyDensityUnit, &
                 TemperatureUnit, LimiterParameter, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, PinnedOption, &
                 UnitOption, VectorIndicesOption )

    class ( Fluid_P_P_Form ), intent ( inout ) :: &
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

    call SetUnits ( VariableUnit, F, MassDensityUnit, EnergyDensityUnit )

    call F % InitializeTemplate_P &
           ( RiemannSolverType, ReconstructedType, UseLimiter, VelocityUnit, &
             MassDensityUnit, EnergyDensityUnit, TemperatureUnit, &
             LimiterParameter, nValues, VariableOption = Variable, &
             VectorOption = VectorOption, NameOption = NameOption, &
             ClearOption = ClearOption, PinnedOption = PinnedOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    !-- Non-standard entropy unit
    F % Unit ( F % ENTROPY_PER_BARYON )  =  UNIT % IDENTITY

  end subroutine InitializeAllocate_P_P


  subroutine SetPrimitiveConserved ( C )

    class ( Fluid_P_P_Form ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_POLYTROPIC ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_POLYTROPIC ) :: &
      ConservedName

    oP = C % N_PRIMITIVE_TEMPLATE + C % N_PRIMITIVE_DUST &
         + C % N_PRIMITIVE_PERFECT
    oC = C % N_CONSERVED_TEMPLATE + C % N_CONSERVED_DUST &
         + C % N_CONSERVED_PERFECT

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_POLYTROPIC
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_POLYTROPIC ) &
      = [ C % INTERNAL_ENERGY ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_POLYTROPIC
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
!    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_POLYTROPIC ) &
!      = [ ]
    
    do iF = 1, C % N_PRIMITIVE_POLYTROPIC
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_POLYTROPIC
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )

    call C % SetPrimitiveConservedTemplate_P ( )

  end subroutine SetPrimitiveConserved


  subroutine SetAdiabaticIndex ( F, AdiabaticIndex )

    class ( Fluid_P_P_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      AdiabaticIndex

    F % AdiabaticIndex = AdiabaticIndex

  end subroutine SetAdiabaticIndex


  subroutine SetFiducialPolytropicParameter ( F, FiducialPolytropicParameter )

    class ( Fluid_P_P_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      FiducialPolytropicParameter

    F % FiducialPolytropicParameter = FiducialPolytropicParameter

  end subroutine SetFiducialPolytropicParameter


  subroutine SetOutput ( F, Output )

    class ( Fluid_P_P_Form ), intent ( in ) :: &
      F
    type ( StorageForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( F % VELOCITY_U )
    call Output % Initialize &
           ( F, iaSelectedOption &
                  = [ F % COMOVING_DENSITY, F % VELOCITY_U, F % PRESSURE, &
                      F % MACH_NUMBER, F % ENTROPY_PER_BARYON ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  subroutine ComputeFromPrimitiveCommon &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( Fluid_P_P_Form ), intent ( in ) :: &
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
        GV => Storage_G % Value, &
        D_S_C => Storage_C % D_Selected, &
        D_S_G => Storage_G % D_Selected, &
        I_DD_22 => G % METRIC_DD_22, &
        I_DD_33 => G % METRIC_DD_33, &        
        I_UU_22 => G % METRIC_UU_22, &
        I_UU_33 => G % METRIC_UU_33 )
        
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
        K     => FV ( oV + 1 : oV + nV, C % POLYTROPIC_PARAMETER ) )

    if ( C % AllocatedDevice ) then

      associate &
        ( D_M_DD_22 => D_S_G ( I_DD_22 ), &
          D_M_DD_33 => D_S_G ( I_DD_33 ), &
          D_M_UU_22 => D_S_G ( I_UU_22 ), &
          D_M_UU_33 => D_S_G ( I_UU_33 ) )
      associate &
        ( D_FEP_1 => D_S_C ( C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
          D_FEP_2 => D_S_C ( C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
          D_FEP_3 => D_S_C ( C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
          D_FEM_1 => D_S_C ( C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
          D_FEM_2 => D_S_C ( C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
          D_FEM_3 => D_S_C ( C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
          D_M     => D_S_C ( C % BARYON_MASS ), &
          D_N     => D_S_C ( C % COMOVING_DENSITY ), &
          D_V_1   => D_S_C ( C % VELOCITY_U ( 1 ) ), &
          D_V_2   => D_S_C ( C % VELOCITY_U ( 2 ) ), &
          D_V_3   => D_S_C ( C % VELOCITY_U ( 3 ) ), &
          D_D     => D_S_C ( C % CONSERVED_DENSITY ), &
          D_S_1   => D_S_C ( C % MOMENTUM_DENSITY_D ( 1 ) ), &
          D_S_2   => D_S_C ( C % MOMENTUM_DENSITY_D ( 2 ) ), &
          D_S_3   => D_S_C ( C % MOMENTUM_DENSITY_D ( 3 ) ), &
          D_E     => D_S_C ( C % INTERNAL_ENERGY ), &
          D_G     => D_S_C ( C % CONSERVED_ENERGY ), &
          D_P     => D_S_C ( C % PRESSURE ), &
          D_Gamma => D_S_C ( C % ADIABATIC_INDEX ), &
          D_CS    => D_S_C ( C % SOUND_SPEED ), &
          D_MN    => D_S_C ( C % MACH_NUMBER ), &
          D_SB    => D_S_C ( C % ENTROPY_PER_BARYON ), &
          D_K     => D_S_C ( C % POLYTROPIC_PARAMETER ) )
          
      call C % ComputeBaryonMassKernelDevice ( M, D_M )
      call C % Apply_EOS_P_KernelDevice &
             ( P, Gamma, SB, K, N, E, &
               D_P, D_Gamma, D_SB, D_K, D_N, D_E, &
               C % AdiabaticIndex, C % FiducialPolytropicParameter )
      call C % ComputeDensityMomentumKernelDevice &
             ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 , &
               D_D, D_S_1, D_S_2, D_S_3, D_N, D_M, D_V_1, D_V_2, D_V_3, &
               D_M_DD_22, D_M_DD_33 )
      call C % ComputeConservedEnergyKernelDevice &
             ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E, &
               D_G, D_M, D_N, D_V_1, D_V_2, D_V_3, D_S_1, D_S_2, D_S_3, D_E )
      call C % ComputeEigenspeedsFluidKernelDevice &
             ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
               M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, &
               M_UU_22, M_UU_33, &
               D_FEP_1, D_FEP_2, D_FEP_3, D_FEM_1, D_FEM_2, D_FEM_3, &
               D_CS, D_MN, D_M, D_N, D_V_1, D_V_2, D_V_3, D_S_1, D_S_2, &
               D_S_3, D_P, D_Gamma, D_M_UU_22, D_M_UU_33 )

      end associate !-- D_FEP_1, etc.
      end associate !-- D_MM_DD_22, etc.

    else

      call C % ComputeBaryonMassKernelHost ( M )
      call C % Apply_EOS_P_KernelHost &
             ( P, Gamma, SB, K, N, E, C % AdiabaticIndex, &
               C % FiducialPolytropicParameter )
      call C % ComputeDensityMomentumKernelHost &
             ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
      call C % ComputeConservedEnergyKernelHost &
             ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
      call C % ComputeEigenspeedsFluidKernelHost &
             ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
               M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, &
               M_UU_22, M_UU_33 )
    
    end if

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
    class ( Fluid_P_P_Form ), intent ( in ) :: &
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
        GV => Storage_G % Value, &
        D_S_C => Storage_C % D_Selected, &
        D_S_G => Storage_G % D_Selected, &
        I_DD_22 => G % METRIC_DD_22, &
        I_DD_33 => G % METRIC_DD_33, &        
        I_UU_22 => G % METRIC_UU_22, &
        I_UU_33 => G % METRIC_UU_33 )

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
        K     => FV ( oV + 1 : oV + nV, C % POLYTROPIC_PARAMETER ) )

    if ( C % AllocatedDevice ) then
    
      associate &
        ( D_M_DD_22 => D_S_G ( I_DD_22 ), &
          D_M_DD_33 => D_S_G ( I_DD_33 ), &
          D_M_UU_22 => D_S_G ( I_UU_22 ), &
          D_M_UU_33 => D_S_G ( I_UU_33 ) )
      associate &
        ( D_FEP_1 => D_S_C ( C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
          D_FEP_2 => D_S_C ( C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
          D_FEP_3 => D_S_C ( C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
          D_FEM_1 => D_S_C ( C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
          D_FEM_2 => D_S_C ( C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
          D_FEM_3 => D_S_C ( C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
          D_M     => D_S_C ( C % BARYON_MASS ), &
          D_N     => D_S_C ( C % COMOVING_DENSITY ), &
          D_V_1   => D_S_C ( C % VELOCITY_U ( 1 ) ), &
          D_V_2   => D_S_C ( C % VELOCITY_U ( 2 ) ), &
          D_V_3   => D_S_C ( C % VELOCITY_U ( 3 ) ), &
          D_D     => D_S_C ( C % CONSERVED_DENSITY ), &
          D_S_1   => D_S_C ( C % MOMENTUM_DENSITY_D ( 1 ) ), &
          D_S_2   => D_S_C ( C % MOMENTUM_DENSITY_D ( 2 ) ), &
          D_S_3   => D_S_C ( C % MOMENTUM_DENSITY_D ( 3 ) ), &
          D_E     => D_S_C ( C % INTERNAL_ENERGY ), &
          D_G     => D_S_C ( C % CONSERVED_ENERGY ), &   
          D_P     => D_S_C ( C % PRESSURE ), &
          D_Gamma => D_S_C ( C % ADIABATIC_INDEX ), &
          D_CS    => D_S_C ( C % SOUND_SPEED ), &
          D_MN    => D_S_C ( C % MACH_NUMBER ), &
          D_SB    => D_S_C ( C % ENTROPY_PER_BARYON ), &
          D_K     => D_S_C ( C % POLYTROPIC_PARAMETER ) )
          
      call C % ComputeBaryonMassKernelDevice ( M, D_M )
      call C % ComputeDensityVelocityKernelDevice &
             ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33, &
               D_N, D_V_1, D_V_2, D_V_3, D_D, D_S_1, D_S_2, D_S_3, D_M, &
               D_M_UU_22, D_M_UU_33 )
      call C % ComputeInternalEnergyKernelDevice &
             ( E, G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, &
               D_E, D_G, D_M, D_N, D_V_1, D_V_2, D_V_3, D_S_1, D_S_2, D_S_3 )
      call C % Apply_EOS_P_KernelDevice &
             ( P, Gamma, SB, K, N, E, &
               D_P, D_Gamma, D_SB, D_K, D_N, D_E, &
               C % AdiabaticIndex, C % FiducialPolytropicParameter )
      call C % ComputeEigenspeedsFluidKernelDevice &
             ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
               M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, &
               M_UU_22, M_UU_33, &
               D_FEP_1, D_FEP_2, D_FEP_3, D_FEM_1, D_FEM_2, D_FEM_3, D_CS, &
               D_MN, D_M, D_N, D_V_1, D_V_2, D_V_3, D_S_1, D_S_2, D_S_3, &
               D_P, D_Gamma, D_M_UU_22, D_M_UU_33 )
      
      end associate   !-- D_FEP_1, etc
      end associate   !-- D_MM_DD_22, etc

    else
    
      call C % ComputeBaryonMassKernelHost ( M )
      call C % ComputeDensityVelocityKernelHost &
             ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33 )
      call C % ComputeInternalEnergyKernelHost &
             ( E, G, M, N, V_1, V_2, V_3, S_1, S_2, S_3 )
      call C % Apply_EOS_P_KernelHost &
             ( P, Gamma, SB, K, N, E, C % AdiabaticIndex, &
               C % FiducialPolytropicParameter )
      call C % ComputeEigenspeedsFluidKernelHost &
             ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
               M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )
    
    end if

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.
    
    if ( associated ( C % Value, Storage_C % Value ) ) &
      call C % Features % Detect ( )

  end subroutine ComputeFromConservedCommon


  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Storage_C, Storage_G, iDimension, &
                 nValuesOption, oValueOption )
    
    class ( StorageForm ), intent ( inout ) :: &
      RawFlux
    class ( Fluid_P_P_Form ), intent ( in ) :: &
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


  subroutine Apply_EOS_P_KernelHost ( P, Gamma, SB, K, N, E, Gamma_0, K0 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      Gamma, &
      SB, &
      K
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      E
    real ( KDR ), intent ( in ) :: &
      Gamma_0, &
      K0

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( P )

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      Gamma ( iV )  =  Gamma_0
      P     ( iV )  =  E ( iV )  *  ( Gamma_0 - 1.0_KDR ) 
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR ) then
        K ( iV )  =  P ( iV ) / ( N ( iV ) ** Gamma_0 )
      else
        K ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      if ( K ( iV ) > 0.0_KDR ) then
        SB ( iV )  =    log ( K ( iV ) / K0 ) / ( Gamma_0 - 1.0_KDR )
      else
        SB ( iV )  =  - 0.1 * huge ( 1.0_KDR )
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine Apply_EOS_P_KernelHost


  subroutine Apply_EOS_P_KernelDevice &
               ( P, Gamma, SB, K, N, E, &
                 D_P, D_Gamma, D_SB, D_K, D_N, D_E, &
                 Gamma_0, K0 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      Gamma, &
      SB, &
      K
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      E
    type ( c_ptr ), intent ( in ) :: &
      D_P, &
      D_Gamma, &
      D_SB, & 
      D_K, & 
      D_N, & 
      D_E
    real ( KDR ), intent ( in ) :: &
      Gamma_0, &
      K0

    integer ( KDI ) :: &
      iV, &
      nValues
      
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_Gamma, Gamma )
    call AssociateHost ( D_SB, SB )
    call AssociateHost ( D_K, K )
    call AssociateHost ( D_N, N ) 
    call AssociateHost ( D_E, E )

    nValues = size ( P )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      Gamma ( iV )  =  Gamma_0
      P     ( iV )  =  E ( iV )  *  ( Gamma_0 - 1.0_KDR ) 
    end do !-- iV
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR ) then
        K ( iV )  =  P ( iV ) / ( N ( iV ) ** Gamma_0 )
      else
        K ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      if ( K ( iV ) > 0.0_KDR ) then
        SB ( iV )  =    log ( K ( iV ) / K0 ) / ( Gamma_0 - 1.0_KDR )
      else
        SB ( iV )  =  - 0.1 * huge ( 1.0_KDR )
      end if
    end do !-- iV
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( E )
    call DisassociateHost ( N )
    call DisassociateHost ( K )
    call DisassociateHost ( SB )
    call DisassociateHost ( Gamma )
    call DisassociateHost ( P )
    
  end subroutine Apply_EOS_P_KernelDevice


  subroutine InitializeBasics &
               ( F, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( Fluid_P_P_Form ), intent ( inout ) :: &
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
      F % Type = 'a Fluid_P_P'

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE + F % N_FIELDS_DUST &
         + F % N_FIELDS_PERFECT
    if ( F % N_FIELDS == 0 ) &
      F % N_FIELDS = oF + F % N_FIELDS_POLYTROPIC

    F % POLYTROPIC_PARAMETER = oF + 1

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_POLYTROPIC ) &
      = [ 'PolytropicParameter            ' ]
    
    !-- units

    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( F % N_FIELDS ) )
      VariableUnit ( oF + 1 : oF + F % N_FIELDS_POLYTROPIC ) = UNIT % IDENTITY
    end if

    !-- vectors

    oV = F % N_VECTORS_TEMPLATE + F % N_VECTORS_DUST &
         + F % N_VECTORS_PERFECT
    if ( F % N_VECTORS == 0 ) F % N_VECTORS = oV + F % N_VECTORS_POLYTROPIC

  end subroutine InitializeBasics
  
    
  subroutine SetUnits ( VariableUnit, F, MassDensityUnit, EnergyDensityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_P_P_Form ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit, &
      EnergyDensityUnit

!    VariableUnit ( F % POLYTROPIC_PARAMETER ) &
!      = EnergyDensityUnit  /  MassDensityUnit ** F % AdiabaticIndex

  end subroutine SetUnits


end module Fluid_P_P__Form
