module Fluid_P__Template

  !-- Fluid_Perfect__Template

  use Basics
  use Mathematics
  use Fluid_D__Form
  use FluidFeatures_Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_PERFECT = 1, &
      N_CONSERVED_PERFECT = 2, &
      N_FIELDS_PERFECT    = 8, &
      N_VECTORS_PERFECT   = 0

  type, public, extends ( Fluid_D_Form ), abstract :: Fluid_P_Template
    integer ( KDI ) :: &
      N_PRIMITIVE_PERFECT = N_PRIMITIVE_PERFECT, &
      N_CONSERVED_PERFECT = N_CONSERVED_PERFECT, &
      N_FIELDS_PERFECT    = N_FIELDS_PERFECT, &
      N_VECTORS_PERFECT   = N_VECTORS_PERFECT, &
      INTERNAL_ENERGY     = 0, &
      CONSERVED_ENERGY    = 0, &
      PRESSURE            = 0, &
      TEMPERATURE         = 0, &
      ENTROPY_PER_BARYON  = 0, &
      CONSERVED_ENTROPY   = 0, &
      SOUND_SPEED         = 0, &
      MACH_NUMBER         = 0
    logical ( KDL ) :: &
      UseEntropy
  contains
    procedure, public, pass :: &
      InitializeTemplate_P
    procedure, public, pass :: &
      SetPrimitiveConservedTemplate_P
    procedure, public, pass ( C ) :: &
      ComputeFluxes
    procedure, public, pass ( C ) :: &
      ComputeCenterStates
    procedure, public, pass ( C ) :: &
      ComputeRawFluxesTemplate_P
    procedure, public, pass ( C ) :: &
      ComputeCenterStatesTemplate_P
    procedure, public, nopass :: &
      ComputeConservedEnergy_G_Kernel
    procedure, public, nopass :: &
      ComputeConservedEntropy_G_Kernel
    procedure, public, nopass :: &
      ComputeInternalEnergy_G_Kernel
    procedure, public, nopass :: &
      ComputeEntropyPerBaryon_G_Kernel
    procedure, public, nopass :: &
      ComputeEigenspeeds_P_G_Kernel
  end type Fluid_P_Template

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeCenterSpeedKernel, &
      ComputeCenterStatesKernel, &
      ComputeFluxes_HLLC_Kernel, &
      ComputeRawFluxesTemplate_P_Kernel

contains


  subroutine InitializeTemplate_P &
               ( F, RiemannSolverType, UseEntropy, UseLimiter, &
                 Velocity_U_Unit, MomentumDensity_D_Unit, BaryonMassUnit, &
                 NumberDensityUnit, EnergyDensityUnit, TemperatureUnit, &
                 BaryonMassReference, LimiterParameter, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

    class ( Fluid_P_Template ), intent ( inout ) :: &
      F
    character ( * ), intent ( in ) :: &
      RiemannSolverType
    logical ( KDL ), intent ( in ) :: &
      UseEntropy, &
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

    call SetUnits &
           ( VariableUnit, F, Velocity_U_Unit, EnergyDensityUnit, &
             TemperatureUnit, NumberDensityUnit )

    call F % Fluid_D_Form % Initialize &
           ( RiemannSolverType, UseLimiter, Velocity_U_Unit, &
             MomentumDensity_D_Unit, BaryonMassUnit, NumberDensityUnit, &
             BaryonMassReference, LimiterParameter, nValues, &
             VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    F % UseEntropy  =  UseEntropy

  end subroutine InitializeTemplate_P


  subroutine SetPrimitiveConservedTemplate_P ( C )

    class ( Fluid_P_Template ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_PERFECT ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_PERFECT ) :: &
      ConservedName

    oP = C % N_PRIMITIVE_TEMPLATE + C % N_PRIMITIVE_DUST
    oC = C % N_CONSERVED_TEMPLATE + C % N_CONSERVED_DUST

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_PERFECT
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_PERFECT ) &
      = [ C % INTERNAL_ENERGY ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_PERFECT
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_PERFECT ) &
      = [ C % CONSERVED_ENERGY, C % CONSERVED_ENTROPY ]
    
    do iF = 1, C % N_PRIMITIVE_PERFECT
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_PERFECT
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )
    
    call C % Fluid_D_Form % SetPrimitiveConserved ( )

  end subroutine SetPrimitiveConservedTemplate_P


  subroutine ComputeFluxes &
               ( Increment, F_I, F_IL, F_IR, SS_I, DF_I, C, Grid, G, &
                 C_IL, C_IR, G_I, iDimension )

    class ( * ), intent ( inout ) :: &
      Increment
    type ( StorageForm ), intent ( inout ) :: &
      F_I, &
      F_IL, F_IR, &
      SS_I, &
      DF_I
    class ( Fluid_P_Template ), intent ( inout ) :: &
      C
    class ( * ), intent ( in ), target :: &
      Grid
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( StorageForm ), intent ( in ) :: &
      C_IL, C_IR, &
      G_I
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iF, &  !-- iField
      iDensity, &
      iMomentum
    real ( KDR ), dimension ( : ), allocatable, target :: &
      M_UU_11
    real ( KDR ), dimension ( : ), pointer :: &
      M_UU, &
      M_DD_22, &
      M_DD_33

    select type ( I => Increment )
    class is ( IncrementDivergence_FV_Form )

    select case ( trim ( C % RiemannSolverType ) )
    case ( 'HLL' )

      call C % ComputeFluxes_HLL &
             ( Increment, F_I, F_IL, F_IR, SS_I, DF_I, Grid, G, &
               C_IL, C_IR, G_I, iDimension )

    case ( 'HLLC' )

      call C % ComputeFluxes_HLL &
             ( Increment, F_I, F_IL, F_IR, SS_I, DF_I, Grid, G, &
               C_IL, C_IR, G_I, iDimension )

      select case ( iDimension )
      case ( 1 )
        allocate ( M_UU_11 ( C % nValues ) )
        !$OMP parallel do private ( iV )
        do iV = 1, C % nValues
          M_UU_11 ( iV )  =  1.0_KDR
        end do !-- iV
        !$OMP end parallel do
        M_UU => M_UU_11
      case ( 2 )
        M_UU => G_I % Value ( :, G % METRIC_UU_22 )
      case ( 3 )
        M_UU => G_I % Value ( :, G % METRIC_UU_33 )
      end select !-- iDimension

      M_DD_22 => G_I % Value ( :, G % METRIC_DD_22 )
      M_DD_33 => G_I % Value ( :, G % METRIC_DD_33 )

      call Search ( C % iaConserved, C % CONSERVED_BARYON_DENSITY, &
                    iDensity )
      call Search ( C % iaConserved, C % MOMENTUM_DENSITY_D ( iDimension ), &
                    iMomentum )
      call ComputeCenterSpeedKernel &
             ( SS_I % Value ( :, C % ALPHA_CENTER ), &
               F_IL % Value ( :, iDensity ), &
               F_IR % Value ( :, iDensity ), &
               F_IL % Value ( :, iMomentum ), &
               F_IR % Value ( :, iMomentum ), &
               C_IL % Value ( :, C % BARYON_MASS ), &
               C_IR % Value ( :, C % BARYON_MASS ), &
               C_IL % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
               C_IR % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
               C_IL % Value ( :, C % MOMENTUM_DENSITY_D ( iDimension ) ), &
               C_IR % Value ( :, C % MOMENTUM_DENSITY_D ( iDimension ) ), &
               SS_I % Value ( :, C % ALPHA_PLUS ), &
               SS_I % Value ( :, C % ALPHA_MINUS ), &
               M_UU )

      associate &
        ( C_ICL => I % Storage % Current_ICL, &
          C_ICR => I % Storage % Current_ICR )

      call C % ComputeCenterStates &
             ( C_ICL, C_ICR, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iDimension )

      !-- Overwrite F_IL and F_IR with F_ICL and F_ICR
      call C % ComputeRawFluxes &
             ( F_IL % Value, G, C_ICL % Value, G_I % Value, iDimension )
      call C % ComputeRawFluxes &
             ( F_IR % Value, G, C_ICR % Value, G_I % Value, iDimension )

      associate ( FF => C % Features )
      do iF = 1, C % N_CONSERVED
        call ComputeFluxes_HLLC_Kernel &
               ( F_I  % Value ( :, iF ), &
                 F_IL % Value ( :, iF ), &
                 F_IR % Value ( :, iF ), &
                 SS_I % Value ( :, C % ALPHA_PLUS ), &
                 SS_I % Value ( :, C % ALPHA_MINUS ), &
                 SS_I % Value ( :, C % ALPHA_CENTER ), &
                 FF   % Value ( :, FF % DIFFUSIVE_FLUX_I ( iDimension ) ) )
      end do !-- iF
      end associate !-- FF

      end associate !-- C_ICL, etc.

    end select !-- RiemannSolverType

    class default
      call Show ( 'Increment type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_P__Template', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeRiemannSolverInput', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Increment

    nullify ( M_UU, M_DD_22, M_DD_33 )

  end subroutine ComputeFluxes


  subroutine ComputeCenterStates &
               ( C_ICL, C_ICR, C, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iD )

    type ( StorageForm ), intent ( inout ) :: &
      C_ICL, C_ICR
    class ( Fluid_P_Template ), intent ( in ) :: &
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

  end subroutine ComputeCenterStates


  subroutine ComputeRawFluxesTemplate_P &
               ( RawFlux, C, G, Value_C, Value_G, iDimension, &
                 nValuesOption, oValueOption )
    
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    class ( Fluid_P_Template ), intent ( in ) :: &
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
      iMomentum, &
      iEnergy, &
      iEntropy, &
      oV, &  !-- oValue
      nV     !-- nValues

    call C % Fluid_D_Form % ComputeRawFluxes &
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
    
    call Search &
           ( C % iaConserved, C % MOMENTUM_DENSITY_D ( iDimension ), &
             iMomentum )
    call Search &
           ( C % iaConserved, C % CONSERVED_ENERGY, iEnergy )
    call Search &
           ( C % iaConserved, C % CONSERVED_ENTROPY, iEntropy )
    
    associate &
      ( F_S_Dim => RawFlux ( oV + 1 : oV + nV, iMomentum ), &
        F_G     => RawFlux ( oV + 1 : oV + nV, iEnergy ), &
        F_DS    => RawFlux ( oV + 1 : oV + nV, iEntropy ), &
        G       => Value_C ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P       => Value_C ( oV + 1 : oV + nV, C % PRESSURE ), &
        DS      => Value_C ( oV + 1 : oV + nV, C % CONSERVED_ENTROPY ), &
        V_Dim   => Value_C ( oV + 1 : oV + nV, C % VELOCITY_U ( iDimension ) ))

    call ComputeRawFluxesTemplate_P_Kernel &
           ( F_S_Dim, F_G, F_DS, G, P, DS, V_Dim )

    end associate !-- F_S_Dim, etc.

  end subroutine ComputeRawFluxesTemplate_P
  

  subroutine ComputeCenterStatesTemplate_P &
               ( C_ICL, C_ICR, C, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iD )

    type ( StorageForm ), intent ( inout ) :: &
      C_ICL, C_ICR
    class ( Fluid_P_Template ), intent ( in ) :: &
      C
    type ( StorageForm ), intent ( in ) :: &
      C_IL, C_IR, &
      SS_I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_DD_22, M_DD_33
    integer ( KDI ), intent ( in ) :: &
      iD

    call ComputeCenterStatesKernel &
           ( C_ICL % Value ( :, C % VELOCITY_U ( 1 ) ), &
             C_ICR % Value ( :, C % VELOCITY_U ( 1 ) ), &
             C_ICL % Value ( :, C % VELOCITY_U ( 2 ) ), &
             C_ICR % Value ( :, C % VELOCITY_U ( 2 ) ), &
             C_ICL % Value ( :, C % VELOCITY_U ( 3 ) ), &
             C_ICR % Value ( :, C % VELOCITY_U ( 3 ) ), &
             C_ICL % Value ( :, C % VELOCITY_U ( iD ) ), &
             C_ICR % Value ( :, C % VELOCITY_U ( iD ) ), &
             C_ICL % Value ( :, C % BARYON_MASS ), &
             C_ICR % Value ( :, C % BARYON_MASS ), &
             C_ICL % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
             C_ICR % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
             C_ICL % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) ), &
             C_ICR % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) ), &
             C_ICL % Value ( :, C % MOMENTUM_DENSITY_D ( 2 ) ), &
             C_ICR % Value ( :, C % MOMENTUM_DENSITY_D ( 2 ) ), &
             C_ICL % Value ( :, C % MOMENTUM_DENSITY_D ( 3 ) ), &
             C_ICR % Value ( :, C % MOMENTUM_DENSITY_D ( 3 ) ), &
             C_ICL % Value ( :, C % MOMENTUM_DENSITY_D ( iD ) ), &
             C_ICR % Value ( :, C % MOMENTUM_DENSITY_D ( iD ) ), &
             C_ICL % Value ( :, C % CONSERVED_ENERGY ), &
             C_ICR % Value ( :, C % CONSERVED_ENERGY ), &
             C_ICL % Value ( :, C % PRESSURE ), &
             C_ICR % Value ( :, C % PRESSURE ), &
             C_ICL % Value ( :, C % CONSERVED_ENTROPY ), &
             C_ICR % Value ( :, C % CONSERVED_ENTROPY ), &
             C_IL % Value ( :, C % VELOCITY_U ( 1 ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( 1 ) ), &
             C_IL % Value ( :, C % VELOCITY_U ( 2 ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( 2 ) ), &
             C_IL % Value ( :, C % VELOCITY_U ( 3 ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( 3 ) ), &
             C_IL % Value ( :, C % VELOCITY_U ( iD ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( iD ) ), &
             C_IL % Value ( :, C % BARYON_MASS ), &
             C_IR % Value ( :, C % BARYON_MASS ), &
             C_IL % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
             C_IR % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
             C_IL % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) ), &
             C_IR % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) ), &
             C_IL % Value ( :, C % MOMENTUM_DENSITY_D ( 2 ) ), &
             C_IR % Value ( :, C % MOMENTUM_DENSITY_D ( 2 ) ), &
             C_IL % Value ( :, C % MOMENTUM_DENSITY_D ( 3 ) ), &
             C_IR % Value ( :, C % MOMENTUM_DENSITY_D ( 3 ) ), &
             C_IL % Value ( :, C % MOMENTUM_DENSITY_D ( iD ) ), &
             C_IR % Value ( :, C % MOMENTUM_DENSITY_D ( iD ) ), &
             C_IL % Value ( :, C % CONSERVED_ENERGY ), &
             C_IR % Value ( :, C % CONSERVED_ENERGY ), &
             C_IL % Value ( :, C % PRESSURE ), &
             C_IR % Value ( :, C % PRESSURE ), &
             C_IL % Value ( :, C % CONSERVED_ENTROPY ), &
             C_IR % Value ( :, C % CONSERVED_ENTROPY ), &
             SS_I % Value ( :, C % ALPHA_PLUS ), &
             SS_I % Value ( :, C % ALPHA_MINUS ), &
             SS_I % Value ( :, C % ALPHA_CENTER ), &
             M_DD_22, M_DD_33 )

  end subroutine ComputeCenterStatesTemplate_P


  subroutine ComputeConservedEnergy_G_Kernel &
               ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )

    !-- Galilean

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      V_1, V_2, V_3, &
      S_1, S_2, S_3, &
      E

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( G )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                            +  S_2 ( iV ) * V_2 ( iV )  &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeConservedEnergy_G_Kernel


  subroutine ComputeConservedEntropy_G_Kernel ( DS, N, SB )
 	 
    real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
      DS
    real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
      N, &
      SB 	 	 
 	 	 
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( DS )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      DS ( iV )  =  SB ( iV )  *  N ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeConservedEntropy_G_Kernel


  subroutine ComputeInternalEnergy_G_Kernel &
               ( E, G, M, N, V_1, V_2, V_3, S_1, S_2, S_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      G
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      M, &
      N, &
      V_1, V_2, V_3, &
      S_1, S_2, S_3

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( E )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      E ( iV )  =  G ( iV )  -  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                            +  S_2 ( iV ) * V_2 ( iV ) &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( E ( iV ) < 0.0_KDR ) then
        E ( iV ) = 0.0_KDR
        G ( iV ) = E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                            +  S_2 ( iV ) * V_2 ( iV ) &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
      end if
    end do
    !$OMP end parallel do

  end subroutine ComputeInternalEnergy_G_Kernel


  subroutine ComputeEntropyPerBaryon_G_Kernel ( SB, DS, N )
 	 
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      SB, &
      DS
    real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
      N 	 	 
 	 	 
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( SB )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR ) then
        SB ( iV ) = DS ( iV ) / N ( iV )
      else
        SB ( iV ) = 0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeEntropyPerBaryon_G_Kernel
  
  
  subroutine ComputeEigenspeeds_P_G_Kernel &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
                 V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33 )

    !-- Galilean

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      MN
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V_1, V_2, V_3, & 
      CS, &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( FEP_1 )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( CS ( iV ) > 0.0_KDR ) then
        MN ( iV ) = sqrt (                      V_1 ( iV ) * V_1 ( iV )  &
                            +  M_DD_22 ( iV ) * V_2 ( iV ) * V_2 ( iV )  &
                            +  M_DD_33 ( iV ) * V_3 ( iV ) * V_3 ( iV ) ) &
                    /  CS ( iV )
      else
        MN ( iV ) = 0.0_KDR
      end if 
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      FEP_1 ( iV )  =  V_1 ( iV )  +  CS ( iV ) 
      FEP_2 ( iV )  =  V_2 ( iV )  +  sqrt ( M_UU_22 ( iV ) ) * CS ( iV ) 
      FEP_3 ( iV )  =  V_3 ( iV )  +  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
      FEM_1 ( iV )  =  V_1 ( iV )  -  CS ( iV )
      FEM_2 ( iV )  =  V_2 ( iV )  -  sqrt ( M_UU_22 ( iV ) ) * CS ( iV )
      FEM_3 ( iV )  =  V_3 ( iV )  -  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeEigenspeeds_P_G_Kernel


  subroutine InitializeBasics &
               ( F, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( Fluid_P_Template ), intent ( inout ) :: &
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
      oV     !-- oVector

    if ( F % Type == '' ) &
      F % Type = 'a Fluid_P'

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE + F % N_FIELDS_DUST
    if ( F % N_FIELDS == 0 ) &
      F % N_FIELDS = oF + F % N_FIELDS_PERFECT

    F % INTERNAL_ENERGY    = oF + 1
    F % CONSERVED_ENERGY   = oF + 2
    F % PRESSURE           = oF + 3
    F % TEMPERATURE        = oF + 4
    F % ENTROPY_PER_BARYON = oF + 5
    F % CONSERVED_ENTROPY  = oF + 6
    F % SOUND_SPEED        = oF + 7
    F % MACH_NUMBER        = oF + 8

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_PERFECT ) &
      = [ 'InternalEnergy  ', &
          'ConservedEnergy ', &
          'Pressure        ', &
          'Temperature     ', &
          'EntropyPerBaryon', &
          'ConservedEntropy', &
          'SoundSpeed      ', &
          'MachNumber      ' ]
    
    !-- units

    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( F % N_FIELDS ) )
      VariableUnit ( oF + 1 : oF + F % N_FIELDS_PERFECT ) = UNIT % IDENTITY
    end if

    !-- vectors

    oV = F % N_VECTORS_TEMPLATE + F % N_VECTORS_DUST
    if ( F % N_VECTORS == 0 ) &
      F % N_VECTORS = oV + F % N_VECTORS_PERFECT

  end subroutine InitializeBasics
  
  
  subroutine SetUnits &
               ( VariableUnit, F, Velocity_U_Unit, EnergyDensityUnit, &
                 TemperatureUnit, NumberDensityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_P_Template ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit, &
      TemperatureUnit, &
      NumberDensityUnit

    VariableUnit ( F % INTERNAL_ENERGY )    = EnergyDensityUnit
    VariableUnit ( F % CONSERVED_ENERGY )   = EnergyDensityUnit
    VariableUnit ( F % PRESSURE )           = EnergyDensityUnit
    VariableUnit ( F % TEMPERATURE )        = TemperatureUnit
    if ( TemperatureUnit /= UNIT % IDENTITY ) then
      VariableUnit ( F % ENTROPY_PER_BARYON ) = UNIT % BOLTZMANN
      VariableUnit ( F % CONSERVED_ENTROPY )  = UNIT % BOLTZMANN &
                                                * NumberDensityUnit
    end if
    VariableUnit ( F % SOUND_SPEED )        = Velocity_U_Unit ( 1 )
    VariableUnit ( F % MACH_NUMBER )        = UNIT % IDENTITY

  end subroutine SetUnits


  subroutine ComputeCenterSpeedKernel &
               ( AC_I, F_D_IL, F_D_IR, F_S_IL, F_S_IR, M_IL, M_IR, &
                 D_IL, D_IR, S_IL, S_IR, AP_I, AM_I, M_UU )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      AC_I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      F_D_IL, F_D_IR, &
      F_S_IL, F_S_IR, &
      M_IL, M_IR, &
      D_IL, D_IR, &
      S_IL, S_IR, &
      AP_I, &
      AM_I, &
      M_UU

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      MD_Numerator, &
      S_Numerator, &
      MD_Numerator_Inv

    nValues = size ( AC_I )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      MD_Numerator &
        =  AP_I ( iV ) * M_IR ( iV ) * D_IR ( iV )  &
           +  AM_I ( iV ) * M_IL ( iV ) * D_IL ( iV ) &
           -  F_D_IR ( iV )  +  F_D_IL ( iV )

      S_Numerator &
        =  AP_I ( iV ) * S_IR ( iV )  +  AM_I ( iV ) * S_IL ( iV ) &
           -  F_S_IR ( iV )  +  F_S_IL ( iV )

      MD_Numerator_Inv  &
        =  max ( MD_Numerator, 0.0_KDR )  &
           /  max ( MD_Numerator ** 2, tiny ( 0.0_KDR ) )

      AC_I ( iV )  &
        =  M_UU ( iV )  *  S_Numerator  *  MD_Numerator_Inv

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeCenterSpeedKernel


  subroutine ComputeCenterStatesKernel &
               ( V_1_ICL, V_1_ICR, V_2_ICL, V_2_ICR, V_3_ICL, V_3_ICR, &
                 V_Dim_ICL, V_Dim_ICR, M_ICL, M_ICR, D_ICL, D_ICR, &
                 S_1_ICL, S_1_ICR, S_2_ICL, S_2_ICR, S_3_ICL, S_3_ICR, &
                 S_Dim_ICL, S_Dim_ICR, G_ICL, G_ICR, P_ICL, P_ICR, &
                 DS_ICL, DS_ICR, V_1_IL, V_1_IR, V_2_IL, V_2_IR, &
                 V_3_IL, V_3_IR, V_Dim_IL, V_Dim_IR, M_IL, M_IR, D_IL, D_IR, &
                 S_1_IL, S_1_IR, S_2_IL, S_2_IR, S_3_IL, S_3_IR, &
                 S_Dim_IL, S_Dim_IR, G_IL, G_IR, P_IL, P_IR, DS_IL, DS_IR, &
                 AP_I, AM_I, AC_I, M_DD_22, M_DD_33 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      V_1_ICL, V_1_ICR, &
      V_2_ICL, V_2_ICR, &
      V_3_ICL, V_3_ICR, &
      V_Dim_ICL, V_Dim_ICR, &
      M_ICL, M_ICR, &
      D_ICL, D_ICR, &
      S_1_ICL, S_1_ICR, &
      S_2_ICL, S_2_ICR, &
      S_3_ICL, S_3_ICR, &
      S_Dim_ICL, S_Dim_ICR, &
      G_ICL, G_ICR, &
      P_ICL, P_ICR, &
      DS_ICL, DS_ICR
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V_1_IL, V_1_IR, &
      V_2_IL, V_2_IR, &
      V_3_IL, V_3_IR, &
      V_Dim_IL, V_Dim_IR, &
      M_IL, M_IR, &
      D_IL, D_IR, &
      S_1_IL, S_1_IR, &
      S_2_IL, S_2_IR, &
      S_3_IL, S_3_IR, &
      S_Dim_IL, S_Dim_IR, &
      G_IL, G_IR, &
      P_IL, P_IR, &
      DS_IL, DS_IR, &
      AP_I, &
      AM_I, &
      AC_I, &
      M_DD_22, M_DD_33

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

      V_1_ICL ( iV )  =  V_1_IL ( iV )
      V_1_ICR ( iV )  =  V_1_IR ( iV )

      V_2_ICL ( iV )  =  V_2_IL ( iV )
      V_2_ICR ( iV )  =  V_2_IR ( iV )

      V_3_ICL ( iV )  =  V_3_IL ( iV )
      V_3_ICR ( iV )  =  V_3_IR ( iV )

      V_Dim_ICL ( iV )  =  AC_I ( iV )
      V_Dim_ICR ( iV )  =  AC_I ( iV )

      AM_VL     =  AM_I ( iV )  +  V_Dim_IL ( iV )
      AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
!      AM_AC_Inv =  1.0_KDR &
!                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
      AM_AC_Inv =  1.0_KDR &
                   / max ( abs ( AM_AC ), SqrtTiny )

      AP_VR     =  AP_I ( iV )  -  V_Dim_IR ( iV )
      AP_AC     =  AP_I ( iV )  -  AC_I ( iV )
!      AP_AC_Inv =  1.0_KDR &
!                   / sign ( max ( abs ( AP_AC ), SqrtTiny ), AP_AC )
      AP_AC_Inv =  1.0_KDR &
                   / max ( abs ( AP_AC ), SqrtTiny )

      M_ICL ( iV )  =  M_IL ( iV )
      M_ICR ( iV )  =  M_IR ( iV )

      D_ICL ( iV )  =  D_IL ( iV ) * AM_VL * AM_AC_Inv
      D_ICR ( iV )  =  D_IR ( iV ) * AP_VR * AP_AC_Inv

      S_1_ICL ( iV )  =  M_ICL ( iV )  *  D_ICL ( iV )  *  V_1_ICL ( iV )
      S_1_ICR ( iV )  =  M_ICR ( iV )  *  D_ICR ( iV )  *  V_1_ICR ( iV )

      S_2_ICL ( iV )  =  M_DD_22 ( iV )  &
                         *  M_ICL ( iV )  *  D_ICL ( iV )  *  V_2_ICL ( iV )
      S_2_ICR ( iV )  =  M_DD_22 ( iV )  &
                         *  M_ICR ( iV )  *  D_ICR ( iV )  *  V_2_ICR ( iV )

      S_3_ICL ( iV )  =  M_DD_33 ( iV )  &
                         *  M_ICL ( iV )  *  D_ICL ( iV )  *  V_3_ICL ( iV )
      S_3_ICR ( iV )  =  M_DD_33 ( iV )  &
                         *  M_ICR ( iV )  *  D_ICR ( iV )  *  V_3_ICR ( iV )

      P_ICL ( iV )  =  P_IL ( iV )  +  S_Dim_IL  ( iV ) * AM_VL &
                                    -  S_Dim_ICL ( iV ) * AM_AC
      P_ICR ( iV )  =  P_IR ( iV )  -  S_Dim_IR  ( iV ) * AP_VR &
                                    +  S_Dim_ICR ( iV ) * AP_AC

      G_ICL ( iV )  =  ( G_IL ( iV ) * AM_VL &
                         +  V_Dim_IL ( iV ) * P_IL ( iV ) &
                         -  AC_I ( iV ) * P_ICL ( iV ) ) &
                       * AM_AC_Inv
      G_ICR ( iV )  =  ( G_IR ( iV ) * AP_VR &
                         -  V_Dim_IR ( iV ) * P_IR ( iV ) &
                         +  AC_I ( iV ) * P_ICR ( iV ) ) &
                       * AP_AC_Inv

      DS_ICL ( iV )  =  DS_IL ( iV ) * AM_VL * AM_AC_Inv
      DS_ICR ( iV )  =  DS_IR ( iV ) * AP_VR * AP_AC_Inv

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeCenterStatesKernel


  subroutine ComputeFluxes_HLLC_Kernel &
               ( F_I, F_ICL, F_ICR, AP_I, AM_I, AC_I, DF_I )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      F_ICL, F_ICR, &
      AP_I, &
      AM_I, &
      AC_I, &
      DF_I

    integer ( KDI ) :: &
      iV, &
      nValues
    
    nValues = size ( F_I )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      !-- If flagged for diffusive flux, leave HLL flux in place
      if ( DF_I ( iV ) > 0.0_KDR ) &
        cycle

      if ( AP_I ( iV ) /= 0.0_KDR .and. AM_I ( iV ) /= 0.0_KDR ) then
        !-- Use the appropriate center state flux
        if ( AC_I ( iV ) >= 0.0_KDR ) then
          F_I ( iV )  =  F_ICL ( iV )
        else !-- AC < 0
          F_I ( iV )  =  F_ICR ( iV )
        end if !-- AC >= 0
      else  !-- AP or AM == 0
        !-- Leave HLL flux in place (which is upwind in this case)
      end if !-- AP and AM /= 0

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeFluxes_HLLC_Kernel


  subroutine ComputeRawFluxesTemplate_P_Kernel &
               ( F_S_Dim, F_G, F_DS, G, P, DS, V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_S_Dim, &
      F_G, &
      F_DS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      G, &
      P, &
      DS, &
      V_Dim

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( F_S_Dim )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      F_S_Dim ( iV )  =  F_S_Dim ( iV )  +  P ( iV )
      F_G     ( iV )  =  ( G ( iV )  +  P ( iV ) ) * V_Dim ( iV )
      F_DS    ( iV )  =  DS ( iV )  *  V_Dim ( iV ) 
    end do
    !$OMP end parallel do

  end subroutine ComputeRawFluxesTemplate_P_Kernel


end module Fluid_P__Template
