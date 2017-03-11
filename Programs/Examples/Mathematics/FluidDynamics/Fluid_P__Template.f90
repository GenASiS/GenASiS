module Fluid_P__Template

  !-- Fluid_Perfect__Template

  use Basics
  use Mathematics
  use Fluid_D__Form
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_PERFECT = 1, &
      N_CONSERVED_PERFECT = 1, &
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
      ADIABATIC_INDEX     = 0, &
      SOUND_SPEED         = 0, &
      MACH_NUMBER         = 0, &
      TEMPERATURE         = 0, &
      ENTROPY_PER_BARYON  = 0
  contains
    procedure, public, pass :: &
      InitializeTemplate_P
    procedure, public, pass ( C ) :: &
      ComputeRawFluxesTemplate_P
    procedure, public, pass ( C ) :: &
      ComputeRiemannSolverInput
    procedure, public, nopass :: &
      ComputeConservedEnergyKernel
    procedure, public, nopass :: &
      ComputeInternalEnergyKernel
    procedure, public, nopass :: &
      ComputeEigenspeedsFluidKernel
  end type Fluid_P_Template

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeCenterSpeed, &
      ComputeRawFluxesTemplate_P_Kernel

  public :: &
    ApplySourcesCurvilinear_Fluid_P

    private :: &
      ApplySourcesCurvilinearKernel
  
contains


  subroutine InitializeTemplate_P &
               ( F, VelocityUnit, MassDensityUnit, EnergyDensityUnit, &
                 TemperatureUnit, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( Fluid_P_Template ), intent ( inout ) :: &
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

    call SetUnits &
           ( VariableUnit, F, VelocityUnit, MassDensityUnit, &
             EnergyDensityUnit, TemperatureUnit )

    call F % Fluid_D_Form % Initialize &
           ( VelocityUnit, MassDensityUnit, nValues, &
             VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeTemplate_P


  subroutine ComputeFromConservedSelf ( C, G, nValuesOption, oValueOption )

    !-- FIXME: Intel compiler does not recognize inheritance from Fluid_D in
    !          extensions of this template

    class ( Fluid_P_Template ), intent ( inout ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromConservedCommon &
           ( C % Value, G, G % Value, nValuesOption, oValueOption )
    
  end subroutine ComputeFromConservedSelf


  subroutine ComputeFromConservedOther &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    !-- FIXME: Intel compiler does not recognize inheritance from Fluid_D in
    !          extensions of this template

    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value_C
    class ( Fluid_P_Template ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromConservedCommon &
           ( Value_C, G, Value_G, nValuesOption, oValueOption )
    
  end subroutine ComputeFromConservedOther


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
      iEnergy
    integer ( KDI ) :: &
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
    
    associate &
      ( F_S_Dim => RawFlux ( oV + 1 : oV + nV, iMomentum ), &
        F_G     => RawFlux ( oV + 1 : oV + nV, iEnergy ), &
        G       => Value_C ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P       => Value_C ( oV + 1 : oV + nV, C % PRESSURE ), &
        V_Dim   => Value_C ( oV + 1 : oV + nV, C % VELOCITY_U ( iDimension ) ))

    call ComputeRawFluxesTemplate_P_Kernel ( F_S_Dim, F_G, G, P, V_Dim )

    end associate !-- F_S_Dim, etc.

  end subroutine ComputeRawFluxesTemplate_P
  

  subroutine ComputeRiemannSolverInput &
               ( Increment, SolverSpeeds_I, DiffusionFactor_I, C, Grid, &
                 C_IL, C_IR, iDimension )
    
    class ( * ), intent ( inout ) :: &
      Increment
    type ( VariableGroupForm ), intent ( inout ) :: &
      SolverSpeeds_I, &
      DiffusionFactor_I
    class ( Fluid_P_Template ), intent ( in ) :: &
      C
    class ( * ), intent ( in ) :: &
      Grid
    type ( VariableGroupForm ), intent ( in ) :: &
      C_IL, &
      C_IR
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iV
    real ( KDR ), dimension ( : ), allocatable, target :: &
      M_UU_11
    real ( KDR ), dimension ( : ), pointer :: &
      M_UU
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( I => Increment )
    class is ( IncrementDivergence_FV_Form )

      select case ( trim ( I % RiemannSolverType ) )
      case ( 'HLL' )

        call C % ComputeRiemannSolverInput_HLL &
               ( SolverSpeeds_I, DiffusionFactor_I, Grid, C_IL, C_IR, &
                 iDimension )

      case ( 'HLLC' )

        call C % ComputeRiemannSolverInput_HLL &
               ( SolverSpeeds_I, DiffusionFactor_I, Grid, C_IL, C_IR, &
                 iDimension )

        select type ( Grid )
        class is ( Chart_SL_Template )
          G => Grid % Geometry ( )
        class default
          call Show ( 'Grid type not recognized', CONSOLE % ERROR )
          call Show ( 'Fluid_P__Template', 'module', CONSOLE % ERROR )
          call Show ( 'ComputeRiemannSolverInput', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- Grid

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
          M_UU => G % Value ( :, G % METRIC_UU_22 )
        case ( 3 )
          M_UU => G % Value ( :, G % METRIC_UU_33 )
        end select !-- iDimension

        call ComputeCenterSpeed &
               ( SolverSpeeds_I % Value ( :, C % ALPHA_CENTER ), &
                 C_IL % Value ( :, C % BARYON_MASS ), &
                 C_IR % Value ( :, C % BARYON_MASS ), &
                 C_IL % Value ( :, C % CONSERVED_DENSITY ), &
                 C_IR % Value ( :, C % CONSERVED_DENSITY ), &
                 C_IL % Value ( :, C % MOMENTUM_DENSITY_D ( iDimension ) ), &
                 C_IR % Value ( :, C % MOMENTUM_DENSITY_D ( iDimension ) ), &
                 C_IL % Value ( :, C % VELOCITY_U ( iDimension ) ), &
                 C_IR % Value ( :, C % VELOCITY_U ( iDimension ) ), &
                 C_IL % Value ( :, C % PRESSURE ), &
                 C_IR % Value ( :, C % PRESSURE ), &
                 SolverSpeeds_I % Value ( :, C % ALPHA_PLUS ), &
                 SolverSpeeds_I % Value ( :, C % ALPHA_MINUS ), &
                 M_UU )

      end select !-- RiemannSolverType

    class default
      call Show ( 'Increment type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_P__Template', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeRiemannSolverInput', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Increment

    nullify ( M_UU, G )

  end subroutine ComputeRiemannSolverInput


  subroutine ComputeConservedEnergyKernel &
               ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E )

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

  end subroutine ComputeConservedEnergyKernel


  subroutine ComputeInternalEnergyKernel &
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

  end subroutine ComputeInternalEnergyKernel


  subroutine ComputeEigenspeedsFluidKernel &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
                 M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, &
                 M_UU_22, M_UU_33 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      CS, &
      MN
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      V_1, V_2, V_3, & 
      S_1, S_2, S_3, &
      P, &
      Gamma, &
      M_UU_22, M_UU_33

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( FEP_1 )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
        CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / ( M ( iV ) * N ( iV ) ) )
      else
        CS ( iV ) = 0.0_KDR
      end if 
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( P ( iV ) > 0.0_KDR ) then
        MN ( iV ) = sqrt ( (    S_1 ( iV ) * V_1 ( iV )  &
                             +  S_2 ( iV ) * V_2 ( iV )  &
                             +  S_3 ( iV ) * V_3 ( iV ) ) &
                           / ( Gamma ( iV ) * P ( iV ) ) )
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

  end subroutine ComputeEigenspeedsFluidKernel


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
      oV, &  !-- oVector
      oP, &  !-- oPrimitive
      oC     !-- oConserved

    if ( F % Type == '' ) &
      F % Type = 'Fluid_P'

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE + F % N_FIELDS_DUST
    if ( F % N_FIELDS == 0 ) F % N_FIELDS = oF + F % N_FIELDS_PERFECT

    F % INTERNAL_ENERGY    = oF + 1
    F % CONSERVED_ENERGY   = oF + 2
    F % PRESSURE           = oF + 3
    F % ADIABATIC_INDEX    = oF + 4
    F % SOUND_SPEED        = oF + 5
    F % MACH_NUMBER        = oF + 6
    F % TEMPERATURE        = oF + 7
    F % ENTROPY_PER_BARYON = oF + 8

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_PERFECT ) &
      = [ 'InternalEnergy                 ', &
          'ConservedEnergy                ', &
          'Pressure                       ', &
          'AdiabaticIndex                 ', &
          'SoundSpeed                     ', &
          'MachNumber                     ', &
          'Temperature                    ', &
          'EntropyPerBaryon               ' ]
    
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
    if ( F % N_VECTORS == 0 ) F % N_VECTORS = oV + F % N_VECTORS_PERFECT

    !-- select primitive, conserved

    oP = F % N_PRIMITIVE_TEMPLATE + F % N_PRIMITIVE_DUST
    oC = F % N_CONSERVED_TEMPLATE + F % N_CONSERVED_DUST

    if ( .not. allocated ( F % iaPrimitive ) ) then
      F % N_PRIMITIVE = oP + F % N_PRIMITIVE_PERFECT
      allocate ( F % iaPrimitive ( F % N_PRIMITIVE ) )
    end if
    F % iaPrimitive ( oP + 1 : oP + F % N_PRIMITIVE_PERFECT ) &
      = [ F % INTERNAL_ENERGY ]

    if ( .not. allocated ( F % iaConserved ) ) then
      F % N_CONSERVED = oC + F % N_CONSERVED_PERFECT
      allocate ( F % iaConserved ( F % N_CONSERVED ) )
    end if
    F % iaConserved ( oC + 1 : oC + F % N_CONSERVED_PERFECT ) &
      = [ F % CONSERVED_ENERGY ]
    
  end subroutine InitializeBasics
  
  
  subroutine SetUnits &
               ( VariableUnit, F, VelocityUnit, MassDensityUnit, &
                 EnergyDensityUnit, TemperatureUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_P_Template ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit

    VariableUnit ( F % INTERNAL_ENERGY )    = EnergyDensityUnit
    VariableUnit ( F % CONSERVED_ENERGY )   = EnergyDensityUnit
    VariableUnit ( F % PRESSURE )           = EnergyDensityUnit
    VariableUnit ( F % ADIABATIC_INDEX )    = UNIT % IDENTITY
    VariableUnit ( F % SOUND_SPEED )        = VelocityUnit ( 1 )
    VariableUnit ( F % MACH_NUMBER )        = UNIT % IDENTITY
    VariableUnit ( F % TEMPERATURE )        = TemperatureUnit
    VariableUnit ( F % ENTROPY_PER_BARYON ) = UNIT % BOLTZMANN

  end subroutine SetUnits


  subroutine ComputeCenterSpeed &
               ( AC_I, M_IL, M_IR, D_IL, D_IR, S_IL, S_IR, V_IL, V_IR, &
                 P_IL, P_IR, AP_I, AM_I, M_UU )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      AC_I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_IL, M_IR, &
      D_IL, D_IR, &
      S_IL, S_IR, &
      V_IL, V_IR, &
      P_IL, P_IR, &
      AP_I, &
      AM_I, &
      M_UU
      
    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      AV_IL, AV_IR, &
      Denominator

    nValues = size ( AC_I )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      AV_IL  =  AM_I ( iV )  +  V_IL ( iV )
      AV_IR  =  AP_I ( iV )  -  V_IR ( iV )
 
      Denominator  =     M_IL ( iV ) * D_IL ( iV ) * AV_IL &
                      +  M_IR ( iV ) * D_IR ( iV ) * AV_IR

      if ( Denominator /= 0.0_KDR ) then
        AC_I ( iV )  =  M_UU ( iV ) &
                        * ( S_IL ( iV ) * AV_IL  +  S_IR ( iV ) * AV_IR  &
                            +  P_IR ( iV ) - P_IL ( iV ) ) &
                        / Denominator
      else
        AC_I ( iV )  =  0.0_KDR
      end if

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeCenterSpeed


  subroutine ComputeRawFluxesTemplate_P_Kernel ( F_S_Dim, F_G, G, P, V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_S_Dim, &
      F_G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      G, &
      P, &
      V_Dim

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( F_S_Dim )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      F_S_Dim ( iV ) = F_S_Dim ( iV )  +  P ( iV )
      F_G     ( iV ) =     ( G ( iV )  +  P ( iV ) ) * V_Dim ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeRawFluxesTemplate_P_Kernel


  subroutine ApplySourcesCurvilinear_Fluid_P ( S, Increment, Fluid, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iMomentum_1, &
      iMomentum_2
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( F => Fluid )
    class is ( Fluid_P_Template )

    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    if ( trim ( Chart % CoordinateSystem ) == 'CARTESIAN' ) &
      return

    call ApplySourcesCurvilinearKernel &
           ( Increment % Value ( :, iMomentum_1 ), &
             Increment % Value ( :, iMomentum_2 ), &
             Chart % CoordinateSystem, Chart % IsProperCell, &
             F % Value ( :, F % PRESSURE ), &
             F % Value ( :, F % MOMENTUM_DENSITY_D ( 2 ) ), &
             F % Value ( :, F % MOMENTUM_DENSITY_D ( 3 ) ), &
             F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             F % Value ( :, F % VELOCITY_U ( 3 ) ), &
             S % dLogVolumeJacobian_dX ( 1 ) % Value, &
             S % dLogVolumeJacobian_dX ( 2 ) % Value, &
             TimeStep, Chart % nDimensions )

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Fluid_P__Template', 'module', CONSOLE % ERROR )
      call Show ( 'ApplySourcesCurvilinear', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    end select !-- F

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplySourcesCurvilinear_Fluid_P


  subroutine ApplySourcesCurvilinearKernel &
               ( KVM_1, KVM_2, CoordinateSystem, IsProperCell, &
                 P, S_2, S_3, V_2, V_3, dLVJ_dX1, dLVJ_dX2, dT, nDimensions )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KVM_1, KVM_2
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      P, &
      S_2, S_3, &
      V_2, V_3, &
      dLVJ_dX1, dLVJ_dX2
    real ( KDR ), intent ( in ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( KVM_1 )

    select case ( trim ( CoordinateSystem ) )
    case ( 'CYLINDRICAL' )

      !$OMP parallel do private ( iV )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        KVM_1 ( iV ) &
          =  KVM_1 ( iV )  +  ( V_3 ( iV ) * S_3 ( iV )  +  P ( iV ) ) &
                              * dLVJ_dX1 ( iV ) * dT
      end do
      !$OMP end parallel do

    case ( 'SPHERICAL' )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        KVM_1 ( iV ) &
          =  KVM_1 ( iV )  +  ( V_2 ( iV ) * S_2 ( iV )  &
                                +  V_3 ( iV ) * S_3 ( iV )  &
                                +  2.0 * P ( iV ) ) &
                              * 0.5_KDR * dLVJ_dX1 ( iV ) * dT
      end do
      !$OMP end parallel do

      if ( nDimensions > 1 ) then
        !$OMP parallel do private ( iV )
        do iV = 1, nV
          if ( .not. IsProperCell ( iV ) ) &
            cycle
          KVM_2 ( iV )  &
            =  KVM_2 ( iV )  +  ( V_3 ( iV ) * S_3 ( iV )  +  P ( iV ) ) &
                                * dLVJ_dX2 ( iV ) * dT
        end do
        !$OMP end parallel do
      end if

    end select !-- CoordinateSystem

  end subroutine ApplySourcesCurvilinearKernel


end module Fluid_P__Template
