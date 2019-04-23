module Fluid_P__Template

#include "Preprocessor"

  !-- Fluid_Perfect__Template

  use iso_c_binding
  use Basics
  use Mathematics
  use Fluid_D__Form
  use FluidFeatures_Template
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_PERFECT = 0, &
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
      ComputeConservedEnergyKernelHost, &
      ComputeConservedEnergyKernelDevice
    procedure, public, nopass :: &
      ComputeInternalEnergyKernelHost, &
      ComputeInternalEnergyKernelDevice
    procedure, public, nopass :: &
      ComputeEigenspeedsFluidKernelHost, &
      ComputeEigenspeedsFluidKernelDevice
  end type Fluid_P_Template

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeCenterSpeedKernel, &
      ComputeCenterStatesKernel, &
      ComputeFluxes_HLLC_Kernel, &
      ComputeRawFluxesTemplate_P_KernelHost, &
      ComputeRawFluxesTemplate_P_KernelDevice

contains


  subroutine InitializeTemplate_P &
               ( F, RiemannSolverType, ReconstructedType, UseLimiter, &
                 VelocityUnit, MassDensityUnit, EnergyDensityUnit, &
                 TemperatureUnit, LimiterParameter, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, PinnedOption, &
                 UnitOption, VectorIndicesOption )

    class ( Fluid_P_Template ), intent ( inout ) :: &
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

    call SetUnits &
           ( VariableUnit, F, VelocityUnit, MassDensityUnit, &
             EnergyDensityUnit, TemperatureUnit )

    call F % Fluid_D_Form % Initialize &
           ( RiemannSolverType, ReconstructedType, UseLimiter, VelocityUnit, &
             MassDensityUnit, LimiterParameter, nValues, &
             VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             PinnedOption = PinnedOption, UnitOption = VariableUnit, &
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
      = [ C % CONSERVED_ENERGY ]
    
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


  ! subroutine ComputeFromConservedSelf ( C, G, nValuesOption, oValueOption )

  !   !-- FIXME: Intel compiler does not recognize inheritance from Fluid_D in
  !   !          extensions of this template

  !   class ( Fluid_P_Template ), intent ( inout ) :: &
  !     C
  !   class ( GeometryFlatForm ), intent ( in ) :: &
  !     G
  !   integer ( KDI ), intent ( in ), optional :: &
  !     nValuesOption, &
  !     oValueOption

  !   call C % ComputeFromConservedCommon &
  !          ( C % Value, G, G % Value, nValuesOption, oValueOption )
    
  ! end subroutine ComputeFromConservedSelf


  ! subroutine ComputeFromConservedOther &
  !              ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

  !   !-- FIXME: Intel compiler does not recognize inheritance from Fluid_D in
  !   !          extensions of this template

  !   real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
  !     Value_C
  !   class ( Fluid_P_Template ), intent ( in ) :: &
  !     C
  !   class ( GeometryFlatForm ), intent ( in ) :: &
  !     G
  !   real ( KDR ), dimension ( :, : ), intent ( in ) :: &
  !     Value_G
  !   integer ( KDI ), intent ( in ), optional :: &
  !     nValuesOption, &
  !     oValueOption

  !   call C % ComputeFromConservedCommon &
  !          ( Value_C, G, Value_G, nValuesOption, oValueOption )
    
  ! end subroutine ComputeFromConservedOther


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

      call Search &
             ( C % iaConserved, C % CONSERVED_DENSITY, &
               iDensity )
      call Search &
             ( C % iaConserved, C % MOMENTUM_DENSITY_D ( iDimension ), &
               iMomentum )
      call ComputeCenterSpeedKernel &
             ( SS_I % Value ( :, C % ALPHA_CENTER ), &
               F_IL % Value ( :, iDensity ), &
               F_IR % Value ( :, iDensity ), &
               F_IL % Value ( :, iMomentum ), &
               F_IR % Value ( :, iMomentum ), &
               C_IL % Value ( :, C % CONSERVED_DENSITY ), &
               C_IR % Value ( :, C % CONSERVED_DENSITY ), &
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
             ( F_IL, G, C_ICL, G_I, iDimension )
      call C % ComputeRawFluxes &
             ( F_IR, G, C_ICR, G_I, iDimension )

      select type ( FF => C % Features )
      class is ( FluidFeaturesTemplate )
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
      end select !-- FF

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
               ( RawFlux, C, G, Storage_C, Storage_G, iDimension, &
                 nValuesOption, oValueOption )
    
    class ( StorageForm ), intent ( inout ) :: &
      RawFlux
    class ( Fluid_P_Template ), intent ( in ) :: &
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

    integer ( KDI ) :: &
      iMomentum, &
      iEnergy
    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues

    call C % Fluid_D_Form % ComputeRawFluxes &
           ( RawFlux, G, Storage_C, Storage_G, iDimension, nValuesOption, &
             oValueOption )
             
    associate &
      ( Value_RF => RawFlux % Value, &
        Value_C  => Storage_C % Value, &
        Value_G  => Storage_G % Value, &
        D_S_RF   => RawFlux % D_Selected, &
        D_S_C    => Storage_C % D_Selected, &
        D_S_G    => Storage_G % D_Selected )

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
      ( F_S_Dim => Value_RF ( oV + 1 : oV + nV, iMomentum ), &
        F_G     => Value_RF ( oV + 1 : oV + nV, iEnergy ), &
        G       => Value_C ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P       => Value_C ( oV + 1 : oV + nV, C % PRESSURE ), &
        V_Dim   => Value_C ( oV + 1 : oV + nV, C % VELOCITY_U ( iDimension ) ))
    
    if ( C % AllocatedDevice ) then
      associate &
        ( D_F_S_Dim => D_S_RF ( iMomentum ), &
          D_F_G     => D_S_RF ( iEnergy ), &
          D_G       => D_S_C ( C % CONSERVED_ENERGY ), &
          D_P       => D_S_C ( C % PRESSURE ), &
          D_V_Dim   => D_S_C ( C % VELOCITY_U ( iDimension ) ))
      call ComputeRawFluxesTemplate_P_KernelDevice &
             ( F_S_Dim, F_G, G, P, V_Dim, &
               D_F_S_Dim, D_F_G, D_G, D_P, D_V_Dim )
      end associate
    else
      call ComputeRawFluxesTemplate_P_KernelHost &
             ( F_S_Dim, F_G, G, P, V_Dim )
    end if

    end associate !-- F_S_Dim, etc.
    
    end associate !-- Value_RF

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
             C_ICL % Value ( :, C % CONSERVED_DENSITY ), &
             C_ICR % Value ( :, C % CONSERVED_DENSITY ), &
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
             C_IL % Value ( :, C % CONSERVED_DENSITY ), &
             C_IR % Value ( :, C % CONSERVED_DENSITY ), &
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


  subroutine ComputeConservedEnergyKernelHost &
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

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues 
      G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                            +  S_2 ( iV ) * V_2 ( iV )  &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeConservedEnergyKernelHost


  subroutine ComputeConservedEnergyKernelDevice &
               ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E , &
                 D_G, D_M, D_N, D_V_1, D_V_2, D_V_3, &
                 D_S_1, D_S_2, D_S_3, D_E )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      V_1, V_2, V_3, &
      S_1, S_2, S_3, &
      E
    type ( c_ptr ), intent ( in ) :: &
      D_G, &
      D_M, &
      D_N, &
      D_V_1, D_V_2, D_V_3, &
      D_S_1, D_S_2, D_S_3, &
      D_E

    integer ( KDI ) :: &
      iV, &
      nValues

    call AssociateHost ( D_G, G )
    call AssociateHost ( D_M, M )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )
    call AssociateHost ( D_E, E )

    nValues = size ( G )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                            +  S_2 ( iV ) * V_2 ( iV )  &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
    end do !-- iV
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    call DisassociateHost ( E )
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( M )
    call DisassociateHost ( G )
    
  end subroutine ComputeConservedEnergyKernelDevice


  subroutine ComputeInternalEnergyKernelHost &
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

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      E ( iV )  =  G ( iV )  -  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                            +  S_2 ( iV ) * V_2 ( iV ) &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
    end do
    !$OMP end parallel do

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      if ( E ( iV ) < 0.0_KDR ) then
        E ( iV ) = 0.0_KDR
        G ( iV ) = E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                            +  S_2 ( iV ) * V_2 ( iV ) &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
      end if
    end do
    !$OMP end parallel do

  end subroutine ComputeInternalEnergyKernelHost


  subroutine ComputeInternalEnergyKernelDevice &
               ( E, G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, &
                 D_E, D_G, D_M, D_N, D_V_1, D_V_2, D_V_3, &
                 D_S_1, D_S_2, D_S_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      G
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      M, &
      N, &
      V_1, V_2, V_3, &
      S_1, S_2, S_3
    type ( c_ptr ), intent ( inout ) :: &
      D_E, &
      D_G, &
      D_M, &
      D_N, &
      D_V_1, D_V_2, D_V_3, &
      D_S_1, D_S_2, D_S_3

    integer ( KDI ) :: &
      iV, &
      nValues
      
    call AssociateHost ( D_E, E )
    call AssociateHost ( D_G, G )
    call AssociateHost ( D_M, M )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )

    nValues = size ( E )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      E ( iV )  =  G ( iV )  -  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                            +  S_2 ( iV ) * V_2 ( iV ) &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
    end do
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      if ( E ( iV ) < 0.0_KDR ) then
        E ( iV ) = 0.0_KDR
        G ( iV ) = E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                            +  S_2 ( iV ) * V_2 ( iV ) &
                                            +  S_3 ( iV ) * V_3 ( iV ) )
      end if
    end do
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( M )
    call DisassociateHost ( G )
    call DisassociateHost ( E )
    
  end subroutine ComputeInternalEnergyKernelDevice


  subroutine ComputeEigenspeedsFluidKernelHost &
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

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
        CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / ( M ( iV ) * N ( iV ) ) )
      else
        CS ( iV ) = 0.0_KDR
      end if 
    end do
    !$OMP end parallel do

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
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

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      FEP_1 ( iV )  =  V_1 ( iV )  +  CS ( iV ) 
      FEP_2 ( iV )  =  V_2 ( iV )  +  sqrt ( M_UU_22 ( iV ) ) * CS ( iV ) 
      FEP_3 ( iV )  =  V_3 ( iV )  +  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
      FEM_1 ( iV )  =  V_1 ( iV )  -  CS ( iV )
      FEM_2 ( iV )  =  V_2 ( iV )  -  sqrt ( M_UU_22 ( iV ) ) * CS ( iV )
      FEM_3 ( iV )  =  V_3 ( iV )  -  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeEigenspeedsFluidKernelHost


  subroutine ComputeEigenspeedsFluidKernelDevice &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, MN, &
                 M, N, V_1, V_2, V_3, S_1, S_2, S_3, P, Gamma, &
                 M_UU_22, M_UU_33, &
                 D_FEP_1, D_FEP_2, D_FEP_3, D_FEM_1, D_FEM_2, D_FEM_3, &
                 D_CS, D_MN, D_M, D_N, D_V_1, D_V_2, D_V_3, D_S_1, D_S_2, &
                 D_S_3, D_P, D_Gamma, D_M_UU_22, D_M_UU_33 )

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
    type ( c_ptr ), intent ( in ) :: &
      D_FEP_1, D_FEP_2, D_FEP_3, &
      D_FEM_1, D_FEM_2, D_FEM_3, &
      D_CS, &
      D_MN, &
      D_M, &
      D_N, &
      D_V_1, D_V_2, D_V_3, &
      D_S_1, D_S_2, D_S_3, &
      D_P, &
      D_Gamma, &
      D_M_UU_22, D_M_UU_33

    integer ( KDI ) :: &
      iV, &
      nValues
      
    call AssociateHost ( D_FEP_1, FEP_1 )
    call AssociateHost ( D_FEP_2, FEP_2 )
    call AssociateHost ( D_FEP_3, FEP_3 )
    call AssociateHost ( D_FEM_1, FEM_1 )
    call AssociateHost ( D_FEM_2, FEM_2 )
    call AssociateHost ( D_FEM_3, FEM_3 )
    call AssociateHost ( D_CS, CS )
    call AssociateHost ( D_MN, MN )
    call AssociateHost ( D_M, M )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_Gamma, Gamma )
    call AssociateHost ( D_M_UU_22, M_UU_22 )
    call AssociateHost ( D_M_UU_33, M_UU_33 )

    nValues = size ( FEP_1 )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
        CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / ( M ( iV ) * N ( iV ) ) )
      else
        CS ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
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

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      FEP_1 ( iV )  =  V_1 ( iV )  +  CS ( iV ) 
      FEP_2 ( iV )  =  V_2 ( iV )  +  sqrt ( M_UU_22 ( iV ) ) * CS ( iV ) 
      FEP_3 ( iV )  =  V_3 ( iV )  +  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
      FEM_1 ( iV )  =  V_1 ( iV )  -  CS ( iV )
      FEM_2 ( iV )  =  V_2 ( iV )  -  sqrt ( M_UU_22 ( iV ) ) * CS ( iV )
      FEM_3 ( iV )  =  V_3 ( iV )  -  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
    end do
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( M_UU_33 )
    call DisassociateHost ( M_UU_22 )
    call DisassociateHost ( Gamma )
    call DisassociateHost ( P )
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( M )
    call DisassociateHost ( MN )
    call DisassociateHost ( CS )
    call DisassociateHost ( FEM_3 )
    call DisassociateHost ( FEM_2 )
    call DisassociateHost ( FEM_1 )
    call DisassociateHost ( FEP_3 )
    call DisassociateHost ( FEP_2 )
    call DisassociateHost ( FEP_1 )

  end subroutine ComputeEigenspeedsFluidKernelDevice


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
      F % Type = 'Fluid_P'

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE + F % N_FIELDS_DUST
    if ( F % N_FIELDS == 0 ) &
      F % N_FIELDS = oF + F % N_FIELDS_PERFECT

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
    if ( F % N_VECTORS == 0 ) &
      F % N_VECTORS = oV + F % N_VECTORS_PERFECT

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


  subroutine ComputeRawFluxesTemplate_P_KernelHost &
               ( F_S_Dim, F_G, G, P, V_Dim )

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

    !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
    do iV = 1, nValues
      F_S_Dim ( iV ) = F_S_Dim ( iV )  +  P ( iV )
      F_G     ( iV ) =     ( G ( iV )  +  P ( iV ) ) * V_Dim ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeRawFluxesTemplate_P_KernelHost

  
  subroutine ComputeRawFluxesTemplate_P_KernelDevice &
               ( F_S_Dim, F_G, G, P, V_Dim, D_F_S_Dim, D_F_G, D_G, D_P, &
                 D_V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_S_Dim, &
      F_G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      G, &
      P, &
      V_Dim
    type ( c_ptr ), intent ( in ) :: &
      D_F_S_Dim, &
      D_F_G, &
      D_G, &
      D_P, &
      D_V_Dim

    integer ( KDI ) :: &
      iV, &
      nValues
      
    call AssociateHost ( D_F_S_Dim, F_S_Dim )
    call AssociateHost ( D_F_G, F_G )
    call AssociateHost ( D_G, G )
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_V_Dim, V_Dim )

    nValues = size ( F_S_Dim )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, nValues
      F_S_Dim ( iV ) = F_S_Dim ( iV )  +  P ( iV )
      F_G     ( iV ) =     ( G ( iV )  +  P ( iV ) ) * V_Dim ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( V_Dim )
    call DisassociateHost ( P )
    call DisassociateHost ( G )
    call DisassociateHost ( F_G )
    call DisassociateHost ( F_S_Dim )
    
  end subroutine ComputeRawFluxesTemplate_P_KernelDevice


end module Fluid_P__Template
