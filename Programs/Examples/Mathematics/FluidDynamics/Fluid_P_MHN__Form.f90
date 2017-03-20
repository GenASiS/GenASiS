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
      N_FIELDS_MEAN_HEAVY_NUCLEUS    = 4, &
      N_VECTORS_MEAN_HEAVY_NUCLEUS   = 0

  type, public, extends ( Fluid_P_Template ) :: Fluid_P_MHN_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_MEAN_HEAVY_NUCLEUS = N_PRIMITIVE_MEAN_HEAVY_NUCLEUS, &
      N_CONSERVED_MEAN_HEAVY_NUCLEUS = N_CONSERVED_MEAN_HEAVY_NUCLEUS, &
      N_FIELDS_MEAN_HEAVY_NUCLEUS    = N_FIELDS_MEAN_HEAVY_NUCLEUS, &
      N_VECTORS_MEAN_HEAVY_NUCLEUS   = N_VECTORS_MEAN_HEAVY_NUCLEUS, &
      COMOVING_ELECTRON_DENSITY      = 0, &
      CONSERVED_ELECTRON_DENSITY     = 0, &
      ELECTRON_FRACTION              = 0, &
      PHASE_TRANSITION               = 0
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
      ComputeConservedElectronKernel
    procedure, public, nopass :: &
      ComputeComovingElectronKernel
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

contains


  subroutine InitializeAllocate_P_MHN &
               ( F, VelocityUnit, MassDensityUnit, EnergyDensityUnit, &
                 NumberDensityUnit, TemperatureUnit, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( Fluid_P_MHN_Form ), intent ( inout ) :: &
      F
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit, &
      EnergyDensityUnit, &
      NumberDensityUnit, &
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

    call SetUnits ( VariableUnit, F, NumberDensityUnit )

    call F % InitializeTemplate_P &
           ( VelocityUnit, MassDensityUnit, EnergyDensityUnit, &
             TemperatureUnit, nValues, VariableOption = Variable, &
             VectorOption = VectorOption, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    call READTABLE &
           ( '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5' )
!    call READTABLE &
!           ( '../Parameters/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5' )

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
      = [ C % TEMPERATURE, C % COMOVING_ELECTRON_DENSITY ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_MEAN_HEAVY_NUCLEUS
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_MEAN_HEAVY_NUCLEUS ) &
      = [ C % CONSERVED_ELECTRON_DENSITY ]
    
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
                      F % SOUND_SPEED, F % MACH_NUMBER, F % TEMPERATURE, &
                      F % ENTROPY_PER_BARYON, F % ELECTRON_FRACTION, &
                      F % ADIABATIC_INDEX, F % PHASE_TRANSITION ], &
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
        NE    => FV ( oV + 1 : oV + nV, C % COMOVING_ELECTRON_DENSITY ), &
        DE    => FV ( oV + 1 : oV + nV, C % CONSERVED_ELECTRON_DENSITY ), &
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ), &
        PT    => FV ( oV + 1 : oV + nV, C % PHASE_TRANSITION ) )

    call C % ComputeBaryonMassKernel ( M )
    call C % Apply_EOS_MHN_T_Kernel &
           ( P, E, Gamma, SB, YE, PT, M, N, T, NE, &
             CONSTANT % ATOMIC_MASS_UNIT )
    call C % ComputeDensityMomentumKernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % ComputeConservedEnergyKernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
    call C % ComputeConservedElectronKernel ( DE, NE )
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
        NE    => FV ( oV + 1 : oV + nV, C % COMOVING_ELECTRON_DENSITY ), &
        DE    => FV ( oV + 1 : oV + nV, C % CONSERVED_ELECTRON_DENSITY ), &
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ), &
        PT    => FV ( oV + 1 : oV + nV, C % PHASE_TRANSITION ) )

    call C % ComputeBaryonMassKernel ( M )
    call C % ComputeDensityVelocityKernel &
           ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33 )
    call C % ComputeInternalEnergyKernel &
           ( E, G, M, N, V_1, V_2, V_3, S_1, S_2, S_3 )
    call C % ComputeComovingElectronKernel ( NE, DE )
    call C % Apply_EOS_MHN_E_Kernel &
           ( P, T, Gamma, SB, YE, PT, M, N, E, NE, &
             CONSTANT % ATOMIC_MASS_UNIT )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.
    
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
      iElectron
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

    call Search ( C % iaConserved, C % CONSERVED_ELECTRON_DENSITY, iElectron )

    associate &
      ( F_DE  => RawFlux ( oV + 1 : oV + nV, iElectron ), &
        DE    => Value_C ( oV + 1 : oV + nV, C % CONSERVED_ELECTRON_DENSITY ), &
        V_Dim => Value_C ( oV + 1 : oV + nV, C % VELOCITY_U ( iDimension ) ) )

    call ComputeRawFluxesKernel ( F_DE, DE, V_Dim )

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
           ( C_ICL % Value ( :, C % CONSERVED_ELECTRON_DENSITY ), &
             C_ICR % Value ( :, C % CONSERVED_ELECTRON_DENSITY ), &
             C_IL % Value ( :, C % CONSERVED_ELECTRON_DENSITY ), &
             C_IR % Value ( :, C % CONSERVED_ELECTRON_DENSITY ), &
             C_IL % Value ( :, C % VELOCITY_U ( iD ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( iD ) ), &
             SS_I % Value ( :, C % ALPHA_PLUS ), &
             SS_I % Value ( :, C % ALPHA_MINUS ), &
             SS_I % Value ( :, C % ALPHA_CENTER ) )

  end subroutine ComputeCenterStates


  subroutine ComputeConservedElectronKernel ( DE, NE )
 	 
    real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
      DE, & 	 	 
      NE 	 	 
 	 	 
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( DE )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( NE ( iV ) < 0.0_KDR ) &
        NE ( iV )  =  0.0_KDR
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      DE ( iV ) = NE ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeConservedElectronKernel


  subroutine ComputeComovingElectronKernel ( NE, DE )
 	 
    real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
      NE, & 	 	 
      DE 	 	 
 	 	 
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( NE )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( DE ( iV ) < 0.0_KDR ) &
        DE ( iV )  =  0.0_KDR
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      NE ( iV ) = DE ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeComovingElectronKernel
  
  
  subroutine Apply_EOS_MHN_T_Kernel &
        ( P, E, Gamma, SB, YE, PT, M, N, T, NE, AMU )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      E, &
      Gamma, &
      SB, &
      YE, &
      PT
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      T, &
      NE
    real ( KDR ), intent ( in ) :: &
      AMU

    integer ( KDI ) :: &
      iV, &
      nValues, &
      keytemp, &
      keyerr
    real ( KDR ) :: &
      rfeps, &
      OR_Shift
    real ( KDR ) :: &
      N_Temp, &
      T_Temp, &
      cs2, dedt, dpderho, dpdrhoe, munu

    nValues = size ( P )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR ) then
        YE ( iV )  =  NE ( iV )  *  M ( iV ) * AMU  /  N ( iV )
      else
        YE ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

    !-- FIXME: Compute P, E, Gamma, SB from N, T, YE

    rfeps = 1.0e-9_KDR
    keytemp = 1_KDI

    !-- Historical Oak Ridge Shift, accounting for nuclear binding energy
    OR_Shift = 8.9_KDR * UNIT % MEV / AMU
    
    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      if ( N ( iV ) == 0.0_KDR ) cycle 
      N_Temp   = N ( iV ) / UNIT % MASS_DENSITY_CGS
      T_Temp   = T ( iV ) / UNIT % MEV 
      E ( iV ) = ( ( E ( iV ) / N ( iV ) ) - OR_Shift ) &
                   / ( UNIT % ERG / UNIT % GRAM )
      P ( iV ) = P ( iV ) / UNIT % BARYE
      call nuc_eos_short &
             ( N_Temp, T_Temp, YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
               cs2, dedt, dpderho, dpdrhoe, munu, &
               keytemp, keyerr, rfeps )
      call nuc_eos_one ( N_Temp, T_Temp, YE ( iV ), Gamma ( iV ), 19 )      
      P ( iV ) = P ( iV ) * UNIT % BARYE
      E ( iV ) = E ( iV ) * UNIT % ERG / UNIT % GRAM + OR_Shift
      E ( iV ) = E ( iV ) * N ( iV )
      if ( Gamma ( iV ) < 1.0_KDR ) then
        PT ( iV ) = 1.0_KDR
      else
        PT ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP end parallel do
    
  end subroutine Apply_EOS_MHN_T_Kernel
  

  subroutine Apply_EOS_MHN_E_Kernel &
               ( P, T, Gamma, SB, YE, PT, M, N, E, NE, AMU )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      T, &
      Gamma, &
      SB, &
      YE, &
      PT
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      E, &
      NE
    real ( KDR ), intent ( in ) :: &
      AMU

    integer ( KDI ) :: &
      iV, &
      nValues, &
      keytemp, &
      keyerr
    real ( KDR ) :: &
      rfeps, &
      OR_Shift
    real ( KDR ) :: &
      N_Temp, &
      E_Temp, &
      cs2, dedt, dpderho, dpdrhoe, munu

    nValues = size ( P )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR ) then
        YE ( iV )  =  NE ( iV )  *  M ( iV ) * AMU  /  N ( iV )
      else
        YE ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

    !-- FIXME: Compute P, T, Gamma, SB from N, T, YE

    rfeps = 1.0e-9_KDR
    keytemp = 0_KDI
    
    !-- Historical Oak Ridge Shift, accounting for nuclear binding energy
    OR_Shift = 8.9_KDR * UNIT % MEV / AMU

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      if ( N ( iV ) == 0.0_KDR ) cycle 
      N_Temp   = N ( iV ) / UNIT % MASS_DENSITY_CGS
      E_Temp   = ( ( E ( iV ) / N ( iV ) ) - OR_Shift ) &
                 / ( UNIT % ERG / UNIT % GRAM )
      P ( iV ) = P ( iV ) / UNIT % BARYE
      T ( iV ) = T ( iV ) / UNIT % MEV
      call nuc_eos_short &
             ( N_Temp, T ( iV ), YE ( iV ), E_temp, P ( iV ), SB ( iV ), &
               cs2, dedt, dpderho, dpdrhoe, munu, &
               keytemp, keyerr, rfeps )
      call nuc_eos_one ( N_Temp, T ( iV ), YE ( iV ), Gamma ( iV ), 19 )
      P ( iV ) = P ( iV ) * UNIT % BARYE
      T ( iV ) = T ( iV ) * UNIT % MEV
      if ( Gamma ( iV ) < 1.0_KDR ) then
        PT ( iV ) = 1.0_KDR
      else
        PT ( iV ) = 0.0_KDR
      end if
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

    F % COMOVING_ELECTRON_DENSITY  = oF + 1
    F % CONSERVED_ELECTRON_DENSITY = oF + 2
    F % ELECTRON_FRACTION          = oF + 3
    F % PHASE_TRANSITION           = oF + 4

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_MEAN_HEAVY_NUCLEUS ) &
      = [ 'ComovingElectronDensity        ', &
          'ConservedElectronDensity       ', &
          'ElectronFraction               ', &
          'PhaseTransition                ' ]
    
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


  subroutine SetUnits ( VariableUnit, F, NumberDensityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_P_MHN_Form ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), intent ( in ) :: &
      NumberDensityUnit

    VariableUnit ( F % COMOVING_ELECTRON_DENSITY ) &
      = NumberDensityUnit
    VariableUnit ( F % CONSERVED_ELECTRON_DENSITY ) &
      = NumberDensityUnit

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
               ( DE_ICL, DE_ICR, DE_IL, DE_IR, V_Dim_IL, V_Dim_IR, &
                 AP_I, AM_I, AC_I )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      DE_ICL, DE_ICR
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      DE_IL, DE_IR, &
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
      AP_VR, &
      AP_AC

    nValues = size ( AC_I )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      AM_VL  =  AM_I ( iV )  +  V_Dim_IL ( iV )
      AM_AC  =  AM_I ( iV )  +  AC_I ( iV )

      AP_VR  =  AP_I ( iV )  -  V_Dim_IR ( iV )
      AP_AC  =  AP_I ( iV )  -  AC_I ( iV )

      DE_ICL ( iV )  =  DE_IL ( iV ) * AM_VL / max ( AM_AC, tiny ( 0.0_KDR ) )
      DE_ICR ( iV )  =  DE_IR ( iV ) * AP_VR / max ( AP_AC, tiny ( 0.0_KDR ) )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeCenterStatesKernel


end module Fluid_P_MHN__Form
