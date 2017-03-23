module Fluid_P_MHN__Form

  !-- Fluid_Perfect_MeanHeavyNucleus__Form

  use Basics
  use Mathematics
  use Fluid_P__Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_MEAN_HEAVY_NUCLEUS = 2, &
      N_CONSERVED_MEAN_HEAVY_NUCLEUS = 1, &
      N_FIELDS_MEAN_HEAVY_NUCLEUS    = 2, &
      N_VECTORS_MEAN_HEAVY_NUCLEUS   = 0

  type, public, extends ( Fluid_P_Template ) :: Fluid_P_MHN_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_MEAN_HEAVY_NUCLEUS = N_PRIMITIVE_MEAN_HEAVY_NUCLEUS, &
      N_CONSERVED_MEAN_HEAVY_NUCLEUS = N_CONSERVED_MEAN_HEAVY_NUCLEUS, &
      N_FIELDS_MEAN_HEAVY_NUCLEUS    = N_FIELDS_MEAN_HEAVY_NUCLEUS, &
      N_VECTORS_MEAN_HEAVY_NUCLEUS   = N_VECTORS_MEAN_HEAVY_NUCLEUS, &
      ELECTRON_FRACTION              = 0, &
      CONSERVED_PROTON_DENSITY       = 0
  contains
    procedure, public, pass :: &
      InitializeAllocate_P_MHN
    generic, public :: &
      Initialize => InitializeAllocate_P_MHN
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, pass ( C ) :: &
      ComputeCenterStates
    procedure, public, nopass :: &
      ComputeConservedProtonKernel
    procedure, public, nopass :: &
      ComputeElectronFractionKernel
    procedure, public, nopass :: &
      Apply_EOS_MHN_T_Kernel
    procedure, public, nopass :: &
      Apply_EOS_MHN_E_Kernel
  end type Fluid_P_MHN_Form

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeRawFluxesKernel, &
      ComputeCenterStatesKernel

    real ( KDR ), private :: &
      OR_Shift, &
      MassDensity_CGS, &
      SpecificEnergy_CGS, &
      Pressure_CGS, &
      MeV

contains


  subroutine InitializeAllocate_P_MHN &
               ( F, RiemannSolverType, UseLimiter, VelocityUnit, &
                 MassDensityUnit, EnergyDensityUnit, TemperatureUnit, &
                 LimiterParameter, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( Fluid_P_MHN_Form ), intent ( inout ) :: &
      F
    character ( * ), intent ( in ) :: &
      RiemannSolverType
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

    call SetUnits ( VariableUnit, F, MassDensityUnit )

    call F % InitializeTemplate_P &
           ( RiemannSolverType, UseLimiter, VelocityUnit, MassDensityUnit, &
             EnergyDensityUnit, TemperatureUnit, LimiterParameter, nValues, &
             VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    call READTABLE &
           ( '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5' )
!    call READTABLE &
!           ( '../Parameters/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5' )

    !-- Historical Oak Ridge Shift, accounting for nuclear binding energy
    OR_Shift = 8.9_KDR * UNIT % MEV / CONSTANT % ATOMIC_MASS_UNIT
    
    MassDensity_CGS     =  UNIT % MASS_DENSITY_CGS
    SpecificEnergy_CGS  =  UNIT % ERG / UNIT % GRAM
    Pressure_CGS        =  UNIT % BARYE
    MeV                 =  UNIT % MEV

  end subroutine InitializeAllocate_P_MHN


  subroutine SetPrimitiveConserved ( C )

    class ( Fluid_P_MHN_Form ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_MEAN_HEAVY_NUCLEUS ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_MEAN_HEAVY_NUCLEUS ) :: &
      ConservedName

    !-- select primitive, conserved

    oP = C % N_PRIMITIVE_TEMPLATE + C % N_PRIMITIVE_DUST &
         + C % N_PRIMITIVE_PERFECT
    oC = C % N_CONSERVED_TEMPLATE + C % N_CONSERVED_DUST &
         + C % N_CONSERVED_PERFECT

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_MEAN_HEAVY_NUCLEUS
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_MEAN_HEAVY_NUCLEUS ) &
      = [ C % TEMPERATURE, C % ELECTRON_FRACTION ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_MEAN_HEAVY_NUCLEUS
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_MEAN_HEAVY_NUCLEUS ) &
      = [ C % CONSERVED_PROTON_DENSITY ]
    
    do iF = 1, C % N_PRIMITIVE_MEAN_HEAVY_NUCLEUS
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_MEAN_HEAVY_NUCLEUS
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )

    call C % SetPrimitiveConservedTemplate_P ( )

  end subroutine SetPrimitiveConserved


  subroutine SetOutput ( F, Output )

    class ( Fluid_P_MHN_Form ), intent ( in ) :: &
      F
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( F % VELOCITY_U )
    call Output % Initialize &
           ( F, iaSelectedOption &
                  = [ F % COMOVING_DENSITY, F % VELOCITY_U, F % PRESSURE, &
                      F % MACH_NUMBER, F % TEMPERATURE, &
                      F % ENTROPY_PER_BARYON, F % ELECTRON_FRACTION ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput
  
  
  subroutine ComputeFromPrimitiveCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( Fluid_P_MHN_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      iV, &
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
        T     => FV ( oV + 1 : oV + nV, C % TEMPERATURE ), &
        SB    => FV ( oV + 1 : oV + nV, C % ENTROPY_PER_BARYON ), &
        DP    => FV ( oV + 1 : oV + nV, C % CONSERVED_PROTON_DENSITY ), &
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ) )

    call C % ComputeBaryonMassKernel &
           ( M )
    call C % Apply_EOS_MHN_T_Kernel &
           ( P, E, Gamma, SB, M, N, T, YE )
    call C % ComputeDensityMomentumKernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % ComputeConservedEnergyKernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
    call C % ComputeConservedProtonKernel &
           ( DP, N, YE )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.
    
    if ( associated ( C % Value, Value_C ) ) &
      call C % Features % Detect ( )

  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( Fluid_P_MHN_Form ), intent ( in ) :: &
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
        T     => FV ( oV + 1 : oV + nV, C % TEMPERATURE ), &
        SB    => FV ( oV + 1 : oV + nV, C % ENTROPY_PER_BARYON ), &
        DP    => FV ( oV + 1 : oV + nV, C % CONSERVED_PROTON_DENSITY ), &
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ) )

    call C % ComputeBaryonMassKernel &
           ( M )
    call C % ComputeDensityVelocityKernel &
           ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33 )
    call C % ComputeInternalEnergyKernel &
           ( E, G, M, N, V_1, V_2, V_3, S_1, S_2, S_3 )
    call C % ComputeElectronFractionKernel &
           ( YE, DP, N )
    call C % Apply_EOS_MHN_E_Kernel &
           ( P, T, Gamma, SB, M, N, E, YE )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.
    
    if ( associated ( C % Value, Value_C ) ) &
      call C % Features % Detect ( )

  end subroutine ComputeFromConservedCommon


  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Value_C, Value_G, iDimension, &
                 nValuesOption, oValueOption )
    
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    class ( Fluid_P_MHN_Form ), intent ( in ) :: &
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
      iProton
    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues

    call C % ComputeRawFluxesTemplate_P &
           ( RawFlux, G, Value_C, Value_G, iDimension, nValuesOption, &
             oValueOption )

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

    call Search ( C % iaConserved, C % CONSERVED_PROTON_DENSITY, iProton )

    associate &
      ( F_DP  => RawFlux ( oV + 1 : oV + nV, iProton ), &
        DP    => Value_C ( oV + 1 : oV + nV, C % CONSERVED_PROTON_DENSITY ), &
        V_Dim => Value_C ( oV + 1 : oV + nV, C % VELOCITY_U ( iDimension ) ) )

    call ComputeRawFluxesKernel ( F_DP, DP, V_Dim )

    end associate !-- F_DE, etc.

  end subroutine ComputeRawFluxes


  subroutine ComputeCenterStates &
               ( C_ICL, C_ICR, C, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iD )

    type ( VariableGroupForm ), intent ( inout ) :: &
      C_ICL, C_ICR
    class ( Fluid_P_MHN_Form ), intent ( in ) :: &
      C
    type ( VariableGroupForm ), intent ( in ) :: &
      C_IL, C_IR, &
      SS_I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_DD_22, M_DD_33
    integer ( KDI ), intent ( in ) :: &
      iD

    call C % ComputeCenterStatesTemplate_P &
           ( C_ICL, C_ICR, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iD )

    call ComputeCenterStatesKernel &
           ( C_ICL % Value ( :, C % CONSERVED_PROTON_DENSITY ), &
             C_ICR % Value ( :, C % CONSERVED_PROTON_DENSITY ), &
             C_IL % Value ( :, C % CONSERVED_PROTON_DENSITY ), &
             C_IR % Value ( :, C % CONSERVED_PROTON_DENSITY ), &
             C_IL % Value ( :, C % VELOCITY_U ( iD ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( iD ) ), &
             SS_I % Value ( :, C % ALPHA_PLUS ), &
             SS_I % Value ( :, C % ALPHA_MINUS ), &
             SS_I % Value ( :, C % ALPHA_CENTER ) )

  end subroutine ComputeCenterStates


  subroutine ComputeConservedProtonKernel ( DP, N, YE )
 	 
    real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
      DP
    real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
      N, &
      YE 	 	 
 	 	 
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( DP )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      DP ( iV )  =  YE ( iV )  *  N ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeConservedProtonKernel


  subroutine ComputeElectronFractionKernel ( YE, DP, N )
 	 
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      YE, &
      DP
    real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
      N 	 	 
 	 	 
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( YE )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( DP ( iV ) < 0.0_KDR ) &
        DP ( iV )  =  0.0_KDR
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR ) then
        YE ( iV ) = DP ( iV ) / N ( iV )
      else
        YE ( iV ) = 0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeElectronFractionKernel
  
  
  subroutine Apply_EOS_MHN_T_Kernel ( P, E, Gamma, SB, M, N, T, YE )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      E, &
      Gamma, &
      SB
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      T, &
      YE

    integer ( KDI ) :: &
      iV, &
      nValues, &
      keytemp, &
      keyerr
    real ( KDR ) :: &
      rfeps
    real ( KDR ) :: &
      N_Temp, &
      T_Temp, &
      cs2, dedt, dpderho, dpdrhoe, munu

    nValues = size ( P )

    !-- Compute P, E, Gamma, SB from N, T, YE

    rfeps = 1.0e-9_KDR
    keytemp = 1_KDI

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      if ( N ( iV ) == 0.0_KDR ) cycle 
      N_Temp   = M ( iV ) * N ( iV ) / MassDensity_CGS
      T_Temp   = T ( iV ) / MeV
      E ( iV ) = ( E ( iV ) / ( M ( iV ) * N ( iV ) )  -  OR_Shift ) &
                 / SpecificEnergy_CGS
      P ( iV ) = P ( iV ) / Pressure_CGS
      call nuc_eos_short &
             ( N_Temp, T_Temp, YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
               cs2, dedt, dpderho, dpdrhoe, munu, &
               keytemp, keyerr, rfeps )
      if ( keyerr /= 0 ) then
        call Show ( 'EOS error', CONSOLE % WARNING )
        call Show ( 'Fluid_P_MHN__Form', 'module', CONSOLE % WARNING )
        call Show ( 'Apply_EOS_MHN_T_Kernel', 'subroutine', CONSOLE % WARNING )
        call Show ( iV, 'iV', CONSOLE % WARNING )
      end if
      call nuc_eos_one ( N_Temp, T_Temp, YE ( iV ), Gamma ( iV ), 19 )      
      P ( iV ) = P ( iV ) * Pressure_CGS
      E ( iV ) = E ( iV ) * SpecificEnergy_CGS  +  OR_Shift
      E ( iV ) = E ( iV ) * M ( iV ) * N ( iV )
    end do
    !$OMP end parallel do
    
  end subroutine Apply_EOS_MHN_T_Kernel
  

  subroutine Apply_EOS_MHN_E_Kernel ( P, T, Gamma, SB, M, N, E, YE )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      T, &
      Gamma, &
      SB
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      E, &
      YE

    integer ( KDI ) :: &
      iV, &
      nValues, &
      keytemp, &
      keyerr
    real ( KDR ) :: &
      rfeps
    real ( KDR ) :: &
      N_Temp, &
      E_Temp, &
      cs2, dedt, dpderho, dpdrhoe, munu

    nValues = size ( P )

    !-- Compute P, T, Gamma, SB from N, E, YE

    rfeps = 1.0e-9_KDR
    keytemp = 0_KDI
    
    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      if ( N ( iV ) == 0.0_KDR ) cycle 
      N_Temp   = M ( iV ) * N ( iV ) / MassDensity_CGS
      E_Temp   = ( E ( iV ) / ( M ( iV ) * N ( iV ) )  -  OR_Shift ) &
                 / SpecificEnergy_CGS
      P ( iV ) = P ( iV ) / Pressure_CGS
      T ( iV ) = T ( iV ) / MeV
      call nuc_eos_short &
             ( N_Temp, T ( iV ), YE ( iV ), E_temp, P ( iV ), SB ( iV ), &
               cs2, dedt, dpderho, dpdrhoe, munu, &
               keytemp, keyerr, rfeps )
      if ( keyerr /= 0 ) then
        call Show ( 'EOS error', CONSOLE % WARNING )
        call Show ( 'Fluid_P_MHN__Form', 'module', CONSOLE % WARNING )
        call Show ( 'Apply_EOS_MHN_E_Kernel', 'subroutine', CONSOLE % WARNING )
        call Show ( iV, 'iV', CONSOLE % WARNING )
      end if
      call nuc_eos_one ( N_Temp, T ( iV ), YE ( iV ), Gamma ( iV ), 19 )
      P ( iV ) = P ( iV ) * Pressure_CGS
      T ( iV ) = T ( iV ) * MeV
    end do
    !$OMP end parallel do
    
  end subroutine Apply_EOS_MHN_E_Kernel


  subroutine InitializeBasics &
               ( F, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( Fluid_P_MHN_Form ), intent ( inout ) :: &
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
      F % Type = 'a Fluid_P_MHN'

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE + F % N_FIELDS_DUST &
         + F % N_FIELDS_PERFECT
    if ( F % N_FIELDS == 0 ) F % N_FIELDS &
      = oF + F % N_FIELDS_MEAN_HEAVY_NUCLEUS

    F % ELECTRON_FRACTION        = oF + 1
    F % CONSERVED_PROTON_DENSITY = oF + 2

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_MEAN_HEAVY_NUCLEUS ) &
      = [ 'ElectronFraction      ', &
          'ConservedProtonDensity' ]
    
    !-- units

    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( F % N_FIELDS ) )
      VariableUnit ( oF + 1 : oF + F % N_FIELDS_MEAN_HEAVY_NUCLEUS ) &
        = UNIT % IDENTITY
    end if

    !-- vectors

    oV = F % N_VECTORS_TEMPLATE + F % N_VECTORS_DUST &
         + F % N_VECTORS_PERFECT
    if ( F % N_VECTORS == 0 ) F % N_VECTORS &
      = oV + F % N_VECTORS_MEAN_HEAVY_NUCLEUS

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, F, MassDensityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_P_MHN_Form ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit

    VariableUnit ( F % CONSERVED_PROTON_DENSITY ) &
      = MassDensityUnit

  end subroutine SetUnits


  subroutine ComputeRawFluxesKernel ( F_DE, DE, V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_DE
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      DE, &
      V_Dim
    
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( F_DE )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      F_DE ( iV )  =  DE ( iV )  *  V_Dim ( iV ) 
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeRawFluxesKernel


  subroutine ComputeCenterStatesKernel &
               ( DP_ICL, DP_ICR, DP_IL, DP_IR, V_Dim_IL, V_Dim_IR, &
                 AP_I, AM_I, AC_I )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      DP_ICL, DP_ICR
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      DP_IL, DP_IR, &
      V_Dim_IL, V_Dim_IR, &
      AP_I, &
      AM_I, &
      AC_I

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      AM_VL, &
      AM_AC, &
      AM_AC_Inv, &
      AP_VR, &
      AP_AC, &
      AP_AC_Inv, &
      SqrtTiny

    nValues = size ( AC_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      AM_VL     =  AM_I ( iV )  +  V_Dim_IL ( iV )
      AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
      AM_AC_Inv =  1.0_KDR &
                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )

      AP_VR     =  AP_I ( iV )  -  V_Dim_IR ( iV )
      AP_AC     =  AP_I ( iV )  -  AC_I ( iV )
      AP_AC_Inv =  1.0_KDR &
                   / sign ( max ( abs ( AP_AC ), SqrtTiny ), AP_AC )

      DP_ICL ( iV )  =  DP_IL ( iV ) * AM_VL * AM_AC_Inv
      DP_ICR ( iV )  =  DP_IR ( iV ) * AP_VR * AP_AC_Inv

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeCenterStatesKernel


end module Fluid_P_MHN__Form
