module Fluid_P_MHN__Form

  !-- Fluid_Perfect_MeanHeavyNucleus__Form

  use Basics
  use Mathematics
  use Fluid_P__Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_MEAN_HEAVY_NUCLEUS = 1, &
      N_CONSERVED_MEAN_HEAVY_NUCLEUS = 1, &
      N_FIELDS_MEAN_HEAVY_NUCLEUS    = 3, &
      N_VECTORS_MEAN_HEAVY_NUCLEUS   = 0

  type, public, extends ( Fluid_P_Template ) :: Fluid_P_MHN_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_MEAN_HEAVY_NUCLEUS = N_PRIMITIVE_MEAN_HEAVY_NUCLEUS, &
      N_CONSERVED_MEAN_HEAVY_NUCLEUS = N_CONSERVED_MEAN_HEAVY_NUCLEUS, &
      N_FIELDS_MEAN_HEAVY_NUCLEUS    = N_FIELDS_MEAN_HEAVY_NUCLEUS, &
      N_VECTORS_MEAN_HEAVY_NUCLEUS   = N_VECTORS_MEAN_HEAVY_NUCLEUS, &
      COMOVING_ELECTRON_DENSITY      = 0, &
      CONSERVED_ELECTRON_DENSITY     = 0, &
      ELECTRON_FRACTION              = 0
  contains
    procedure, public, pass :: &
      InitializeAllocate_P_MHN
    generic, public :: &
      Initialize => InitializeAllocate_P_MHN
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass :: & 
      ComputeFromProfile
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, nopass :: &
      ComputeConservedElectronKernel
    procedure, public, nopass :: &
      ComputeComovingElectronKernel
    procedure, public, nopass :: &
      Apply_EOS_MHN_T_Kernel
    procedure, public, nopass :: &
      Apply_EOS_MHN_E_Kernel
    procedure, public, nopass :: &
      ResetTemperature
  end type Fluid_P_MHN_Form

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeRawFluxesKernel

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

  end subroutine InitializeAllocate_P_MHN


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
                      F % ENTROPY_PER_BARYON, F % ELECTRON_FRACTION ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput
  
  
  subroutine ComputeFromProfile ( F, G, nValuesOption, oValueOption ) 

    class ( Fluid_P_MHN_Form  ), intent ( inout ) :: &
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
        T     => FV ( oV + 1 : oV + nV, F % TEMPERATURE ), &
        SB    => FV ( oV + 1 : oV + nV, F % ENTROPY_PER_BARYON ), &
        NE    => FV ( oV + 1 : oV + nV, F % COMOVING_ELECTRON_DENSITY ), &
        DE    => FV ( oV + 1 : oV + nV, F % CONSERVED_ELECTRON_DENSITY ), &
        YE    => FV ( oV + 1 : oV + nV, F % ELECTRON_FRACTION ) )

    call F % ComputeBaryonMassKernel ( M )
    call F % Apply_EOS_MHN_T_Kernel &
           ( P, E, Gamma, SB, YE, M, N, T, NE, UNIT % ATOMIC_MASS_UNIT )
    call F % ComputeDensityMomentumKernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call F % ComputeConservedEnergyKernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
    call F % ComputeConservedElectronKernel ( DE, NE )
    call F % Apply_EOS_MHN_E_Kernel &
           ( P, T, Gamma, SB, YE, M, N, E, NE, UNIT % ATOMIC_MASS_UNIT )
    call F % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromProfile


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
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ) )

    call C % ResetTemperature ( C % Value ( :, C % TEMPERATURE ), T )

    call C % ComputeBaryonMassKernel ( M )
    call C % ComputeDensityMomentumKernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % ComputeConservedEnergyKernel &
           ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )
    call C % ComputeConservedElectronKernel ( DE, NE )
    call C % Apply_EOS_MHN_E_Kernel &
           ( P, T, Gamma, SB, YE, M, N, E, NE, UNIT % ATOMIC_MASS_UNIT )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( C, G, Value, nValuesOption, oValueOption )

    class ( Fluid_P_MHN_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues

    associate &
      ( FV => Value, &
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
        YE    => FV ( oV + 1 : oV + nV, C % ELECTRON_FRACTION ) )

    call C % ComputeBaryonMassKernel ( M )
    call C % ComputeDensityVelocityKernel &
           ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33 )
    call C % ComputeInternalEnergyKernel &
           ( E, G, M, N, V_1, V_2, V_3, S_1, S_2, S_3 )
    call C % ComputeComovingElectronKernel ( NE, DE )
    call C % Apply_EOS_MHN_E_Kernel &
           ( P, T, Gamma, SB, YE, M, N, E, NE, UNIT % ATOMIC_MASS_UNIT )
    call C % ComputeEigenspeedsFluidKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
             M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, M_UU_22, M_UU_33 )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.
    
  end subroutine ComputeFromConservedCommon


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
        ( P, E, Gamma, SB, YE, M, N, T, NE, AMU )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      E, &
      Gamma, &
      SB, &
      YE
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      T, &
      NE
    type ( MeasuredValueForm ), intent ( in ) :: &
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
    end do
    !$OMP end parallel do
    
  end subroutine Apply_EOS_MHN_T_Kernel
  

  subroutine Apply_EOS_MHN_E_Kernel ( P, T, Gamma, SB, YE, M, N, E, NE, AMU )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      T, &
      Gamma, &
      SB, &
      YE
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      E, &
      NE
    type ( MeasuredValueForm ), intent ( in ) :: &
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
    end do
    !$OMP end parallel do
    
  end subroutine Apply_EOS_MHN_E_Kernel


  subroutine ResetTemperature ( T_Value, T )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      T_Value
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      T

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( T )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      T ( iV )  =  T_Value ( iV )
    end do
    !$OMP end parallel do

  end subroutine ResetTemperature


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
          'ElectronFraction               ' ]
    
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

    !-- select primitive, conserved

    oP = F % N_PRIMITIVE_TEMPLATE + F % N_PRIMITIVE_DUST &
         + F % N_PRIMITIVE_PERFECT
    oC = F % N_CONSERVED_TEMPLATE + F % N_CONSERVED_DUST &
         + F % N_CONSERVED_PERFECT

    if ( .not. allocated ( F % iaPrimitive ) ) then
      F % N_PRIMITIVE = oP + F % N_PRIMITIVE_MEAN_HEAVY_NUCLEUS
      allocate ( F % iaPrimitive ( F % N_PRIMITIVE ) )
    end if
    F % iaPrimitive ( oP + 1 : oP + F % N_PRIMITIVE_MEAN_HEAVY_NUCLEUS ) &
      = [ F % COMOVING_ELECTRON_DENSITY ]

    if ( .not. allocated ( F % iaConserved ) ) then
      F % N_CONSERVED = oC + F % N_CONSERVED_MEAN_HEAVY_NUCLEUS
      allocate ( F % iaConserved ( F % N_CONSERVED ) )
    end if
    F % iaConserved ( oC + 1 : oC + F % N_CONSERVED_MEAN_HEAVY_NUCLEUS ) &
      = [ F % CONSERVED_ELECTRON_DENSITY ]
    
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


end module Fluid_P_MHN__Form
