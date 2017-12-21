module Fluid_P__Template

  !-- Fluid_Perfect__Template

  use Basics
  use Mathematics
  use Fluid_D__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_PERFECT = 0, &
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
  contains
    procedure, public, pass :: &
      InitializeTemplate_P
    procedure, public, pass :: &
      SetPrimitiveConservedTemplate_P
    procedure, public, pass ( C ) :: &
      ComputeFluxes
    procedure, public, pass ( C ) :: &
      ComputeCenterStates
  !   procedure, public, pass ( C ) :: &
  !     ComputeRawFluxesTemplate_P
    procedure, public, pass ( C ) :: &
      ComputeCenterStatesTemplate_P
  !   procedure, public, nopass :: &
  !     ComputeConservedEnergyKernel
  !   procedure, public, nopass :: &
  !     ComputeInternalEnergyKernel
  !   procedure, public, nopass :: &
  !     ComputeEigenspeedsFluidKernel
  end type Fluid_P_Template

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeCenterSpeedKernel, &
      ComputeCenterStatesKernel!, &
      ! ComputeFluxes_HLLC_Kernel, &
      ! ComputeRawFluxesTemplate_P_Kernel

contains


  subroutine InitializeTemplate_P &
               ( F, RiemannSolverType, UseLimiter, Velocity_U_Unit, &
                 MomentumDensity_D_Unit, BaryonMassUnit, NumberDensityUnit, &
                 EnergyDensityUnit, TemperatureUnit, BaryonMassReference, &
                 LimiterParameter, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( Fluid_P_Template ), intent ( inout ) :: &
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
!    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_PERFECT ) &
!      = [ ]

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
    type ( VariableGroupForm ), intent ( inout ) :: &
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
    type ( VariableGroupForm ), intent ( in ) :: &
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

    !   select type ( FF => C % Features )
    !   class is ( FluidFeaturesTemplate )
    !   do iF = 1, C % N_CONSERVED
    !     call ComputeFluxes_HLLC_Kernel &
    !            ( F_I  % Value ( :, iF ), &
    !              F_IL % Value ( :, iF ), &
    !              F_IR % Value ( :, iF ), &
    !              SS_I % Value ( :, C % ALPHA_PLUS ), &
    !              SS_I % Value ( :, C % ALPHA_MINUS ), &
    !              SS_I % Value ( :, C % ALPHA_CENTER ), &
    !              FF   % Value ( :, FF % DIFFUSIVE_FLUX_I ( iDimension ) ) )
    !   end do !-- iF
    !   end select !-- FF

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

    type ( VariableGroupForm ), intent ( inout ) :: &
      C_ICL, C_ICR
    class ( Fluid_P_Template ), intent ( in ) :: &
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

  end subroutine ComputeCenterStates


  subroutine ComputeCenterStatesTemplate_P &
               ( C_ICL, C_ICR, C, C_IL, C_IR, SS_I, M_DD_22, M_DD_33, iD )

    type ( VariableGroupForm ), intent ( inout ) :: &
      C_ICL, C_ICR
    class ( Fluid_P_Template ), intent ( in ) :: &
      C
    type ( VariableGroupForm ), intent ( in ) :: &
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
             C_IL % Value ( :, C % VELOCITY_U ( 1 ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( 1 ) ), &
             C_IL % Value ( :, C % VELOCITY_U ( 2 ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( 2 ) ), &
             C_IL % Value ( :, C % VELOCITY_U ( 3 ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( 3 ) ), &
             C_IL % Value ( :, C % VELOCITY_U ( iD ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( iD ) ), &
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
             SS_I % Value ( :, C % ALPHA_PLUS ), &
             SS_I % Value ( :, C % ALPHA_MINUS ), &
             SS_I % Value ( :, C % ALPHA_CENTER ), &
             M_DD_22, M_DD_33 )

  end subroutine ComputeCenterStatesTemplate_P


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
               ( AC_I, F_D_IL, F_D_IR, F_S_IL, F_S_IR, D_IL, D_IR, &
                 S_IL, S_IR, AP_I, AM_I, M_UU )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      AC_I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      F_D_IL, F_D_IR, &
      F_S_IL, F_S_IR, &
      D_IL, D_IR, &
      S_IL, S_IR, &
      AP_I, &
      AM_I, &
      M_UU

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      D_Numerator, &
      S_Numerator, &
      D_Numerator_Inv

    nValues = size ( AC_I )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      D_Numerator &
        =  AP_I ( iV ) * D_IR ( iV )  +  AM_I ( iV ) * D_IL ( iV ) &
           -  F_D_IR ( iV )  +  F_D_IL ( iV )

      S_Numerator &
        =  AP_I ( iV ) * S_IR ( iV )  +  AM_I ( iV ) * S_IL ( iV ) &
           -  F_S_IR ( iV )  +  F_S_IL ( iV )

      D_Numerator_Inv  &
        =  max ( D_Numerator, 0.0_KDR )  &
           /  max ( D_Numerator ** 2, tiny ( 0.0_KDR ) )

      AC_I ( iV )  &
        =  M_UU ( iV )  *  S_Numerator  *  D_Numerator_Inv

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeCenterSpeedKernel


  subroutine ComputeCenterStatesKernel &
               ( V_1_ICL, V_1_ICR, V_2_ICL, V_2_ICR, V_3_ICL, V_3_ICR, &
                 V_Dim_ICL, V_Dim_ICR, D_ICL, D_ICR, S_1_ICL, S_1_ICR, &
                 S_2_ICL, S_2_ICR, S_3_ICL, S_3_ICR, S_Dim_ICL, S_Dim_ICR, &
                 G_ICL, G_ICR, P_ICL, P_ICR, &
                 V_1_IL, V_1_IR, V_2_IL, V_2_IR, V_3_IL, V_3_IR, &
                 V_Dim_IL, V_Dim_IR, D_IL, D_IR, S_1_IL, S_1_IR, &
                 S_2_IL, S_2_IR, S_3_IL, S_3_IR, S_Dim_IL, S_Dim_IR, &
                 G_IL, G_IR, P_IL, P_IR, &
                 AP_I, AM_I, AC_I, M_DD_22, M_DD_33 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      V_1_ICL, V_1_ICR, &
      V_2_ICL, V_2_ICR, &
      V_3_ICL, V_3_ICR, &
      V_Dim_ICL, V_Dim_ICR, &
      D_ICL, D_ICR, &
      S_1_ICL, S_1_ICR, &
      S_2_ICL, S_2_ICR, &
      S_3_ICL, S_3_ICR, &
      S_Dim_ICL, S_Dim_ICR, &
      G_ICL, G_ICR, &
      P_ICL, P_ICR
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V_1_IL, V_1_IR, &
      V_2_IL, V_2_IR, &
      V_3_IL, V_3_IR, &
      V_Dim_IL, V_Dim_IR, &
      D_IL, D_IR, &
      S_1_IL, S_1_IR, &
      S_2_IL, S_2_IR, &
      S_3_IL, S_3_IR, &
      S_Dim_IL, S_Dim_IR, &
      G_IL, G_IR, &
      P_IL, P_IR, &
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

      D_ICL ( iV )  =  D_IL ( iV ) * AM_VL * AM_AC_Inv
      D_ICR ( iV )  =  D_IR ( iV ) * AP_VR * AP_AC_Inv

      S_1_ICL ( iV )  =  D_ICL ( iV )  *  V_1_ICL ( iV )
      S_1_ICR ( iV )  =  D_ICR ( iV )  *  V_1_ICR ( iV )

      S_2_ICL ( iV )  =  M_DD_22 ( iV )  *  D_ICL ( iV )  *  V_2_ICL ( iV )
      S_2_ICR ( iV )  =  M_DD_22 ( iV )  *  D_ICR ( iV )  *  V_2_ICR ( iV )

      S_3_ICL ( iV )  =  M_DD_33 ( iV )  *  D_ICL ( iV )  *  V_3_ICL ( iV )
      S_3_ICR ( iV )  =  M_DD_33 ( iV )  *  D_ICR ( iV )  *  V_3_ICR ( iV )

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

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeCenterStatesKernel


end module Fluid_P__Template
