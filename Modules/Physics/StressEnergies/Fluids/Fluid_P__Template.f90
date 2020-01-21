module Fluid_P__Template

  !-- Fluid_Perfect__Template

  use Basics
  use Mathematics
  use StressEnergyBasics
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
    procedure ( CFT ), public, pass ( C ), deferred :: &
      ComputeFromTemperature
    procedure, public, pass ( C ) :: &
      ComputeFluxes
    procedure, public, pass ( C ) :: &
      ComputeCenterStates
    procedure, public, pass ( C ) :: &
      ComputeRawFluxesTemplate_P
    procedure, public, pass ( C ) :: &
      ComputeCenterStatesTemplate_P
    procedure, public, nopass :: &
      Compute_G_G_Kernel
    procedure, public, nopass :: &
      Compute_DS_G_Kernel
    procedure, public, nopass :: &
      Compute_N_V_E_G_Kernel
    procedure, public, nopass :: &
      Compute_SB_G_Kernel
    procedure, public, nopass :: &
      Compute_FE_P_G_Kernel
  end type Fluid_P_Template

  abstract interface
    subroutine CFT ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )
      use Basics
      use Mathematics
      import Fluid_P_Template
      class ( StorageForm ), intent ( inout ), target :: &
        Storage_C
      class ( Fluid_P_Template ), intent ( in ) :: &
        C
      class ( GeometryFlatForm ), intent ( in ) :: &
        G
      class ( StorageForm ), intent ( in ) :: &
        Storage_G
      integer ( KDI ), intent ( in ), optional :: &
        nValuesOption, &
        oValueOption
    end subroutine CFT
  end interface

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeCenterSpeedKernel, &
      ComputeCenterStatesKernel, &
      ComputeFluxes_HLLC_Kernel, &
      ComputeRawFluxesTemplate_P_Kernel
      
  interface
  
    module subroutine Compute_G_G_Kernel &
                 ( G, M, N, V_1, V_2, V_3, S_1, S_2, S_3, E, UseDeviceOption )
      !-- Compute_ConservedEnergy_Galilean_Kernel
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        G
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N, &
        V_1, V_2, V_3, &
        S_1, S_2, S_3, &
        E
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_G_G_Kernel


    module subroutine Compute_DS_G_Kernel ( DS, N, SB, UseDeviceOption )
      !-- Compute_ConservedEntropy_Galilean_Kernel
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
        DS
      real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
        N, &
        SB 	 	 
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_DS_G_Kernel


    module subroutine Compute_N_V_E_G_Kernel &
                 ( N, V_1, V_2, V_3, E, D, S_1, S_2, S_3, G, M, &
                   M_UU_22, M_UU_33, N_Min, UseDeviceOption )
      !-- Compute_ComovingBaryonDensity_Velocity_InternalEnergyDensity_Galilean
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        N, &
        V_1, V_2, V_3, &
        E, &
        D, &
        S_1, S_2, S_3, &
        G
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        M_UU_22, M_UU_33
      real ( KDR ), intent ( in ) :: &
        N_Min
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_N_V_E_G_Kernel


    module subroutine Compute_SB_G_Kernel ( SB, DS, N, UseDeviceOption )
      !-- Compute_EntropyPerBaryon_Galilean_Kernel
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        SB, &
        DS
      real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
        N 	 	 
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_SB_G_Kernel
    
    
    module subroutine Compute_FE_P_G_Kernel &
                 ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, MN, &
                   V_1, V_2, V_3, CS, M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
                   UseDeviceOption )
      !-- Compute_FastEigenspeeds_Perfect_Galilean_Kernel
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        FEP_1, FEP_2, FEP_3, &
        FEM_1, FEM_2, FEM_3, &
        MN
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        V_1, V_2, V_3, & 
        CS, &
        M_DD_22, M_DD_33, &
        M_UU_22, M_UU_33
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_FE_P_G_Kernel
    
    
    module subroutine ComputeCenterSpeedKernel &
                 ( AC_I, F_D_IL, F_D_IR, F_S_IL, F_S_IR, M_IL, M_IR, &
                   D_IL, D_IR, S_IL, S_IR, AP_I, AM_I, M_UU, UseDeviceOption )
      use Basics
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
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeCenterSpeedKernel


    module subroutine ComputeCenterStatesKernel &
                 ( V_ICL, V_ICR, S_ICL, S_ICR, M_ICL, M_ICR, D_ICL, D_ICR, &
                   G_ICL, G_ICR, P_ICL, P_ICR, DS_ICL, DS_ICR, &
                   V_IL, V_IR, S_IL, S_IR, M_IL, M_IR, D_IL, D_IR, &
                   G_IL, G_IR, P_IL, P_IR, DS_IL, DS_IR, &
                   AP_I, AM_I, AC_I, M_DD_22, M_DD_33, iD, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
        V_ICL, V_ICR, &
        S_ICL, S_ICR
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        M_ICL, M_ICR, &
        D_ICL, D_ICR, &
        G_ICL, G_ICR, &
        P_ICL, P_ICR, &
        DS_ICL, DS_ICR
      real ( KDR ), dimension ( :, : ), intent ( in ) :: &
        V_IL, V_IR, &
        S_IL, S_IR
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M_IL, M_IR, &
        D_IL, D_IR, &
        G_IL, G_IR, &
        P_IL, P_IR, &
        DS_IL, DS_IR, &
        AP_I, &
        AM_I, &
        AC_I, &
        M_DD_22, M_DD_33
      integer ( KDI ), intent ( in ) :: &
        iD
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeCenterStatesKernel


    module subroutine ComputeFluxes_HLLC_Kernel &
                 ( F_I, F_ICL, F_ICR, AP_I, AM_I, AC_I, DF_I, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        F_I
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        F_ICL, F_ICR, &
        AP_I, &
        AM_I, &
        AC_I, &
        DF_I
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeFluxes_HLLC_Kernel


    module subroutine ComputeRawFluxesTemplate_P_Kernel &
                 ( F_S_Dim, F_G, F_DS, G, P, DS, V_Dim, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        F_S_Dim, &
        F_G, &
        F_DS
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        G, &
        P, &
        DS, &
        V_Dim
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeRawFluxesTemplate_P_Kernel

  end interface

contains


  subroutine InitializeTemplate_P &
               ( F, FluidType, RiemannSolverType, ReconstructedType, &
                 UseEntropy, UseLimiter, Units, BaryonMassReference, &
                 LimiterParameter, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( Fluid_P_Template ), intent ( inout ) :: &
      F
    character ( * ), intent ( in ) :: &
      FluidType, &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseEntropy, &
      UseLimiter
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
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

    call SetUnits ( VariableUnit, F, Units )

    call F % Initialize_D &
           ( FluidType, RiemannSolverType, ReconstructedType, UseLimiter, &
             Units, BaryonMassReference, LimiterParameter, nValues, &
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
             ( F_IL, G, C_ICL, G_I, iDimension )
      call C % ComputeRawFluxes &
             ( F_IR, G, C_ICR, G_I, iDimension )

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
      iEnergy, &
      iEntropy, &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( Value_RF => RawFlux % Value, &
        Value_C  => Storage_C % Value )

    call C % Fluid_D_Form % ComputeRawFluxes &
           ( RawFlux, G, Storage_C, Storage_G, iDimension, nValuesOption, &
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
      ( F_S_Dim => Value_RF ( oV + 1 : oV + nV, iMomentum ), &
        F_G     => Value_RF ( oV + 1 : oV + nV, iEnergy ), &
        F_DS    => Value_RF ( oV + 1 : oV + nV, iEntropy ), &
        G       => Value_C ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        P       => Value_C ( oV + 1 : oV + nV, C % PRESSURE ), &
        DS      => Value_C ( oV + 1 : oV + nV, C % CONSERVED_ENTROPY ), &
        V_Dim   => Value_C ( oV + 1 : oV + nV, C % VELOCITY_U ( iDimension ) ))

    call ComputeRawFluxesTemplate_P_Kernel &
           ( F_S_Dim, F_G, F_DS, G, P, DS, V_Dim )

    end associate !-- F_S_Dim, etc.
    
    end associate !-- Value_RF, Value_C

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
           ( C_ICL % Value ( :, C % VELOCITY_U ( 1 ) &
                                  : C % VELOCITY_U ( 3 ) ), &
             C_ICR % Value ( :, C % VELOCITY_U ( 1 ) &
                                  : C % VELOCITY_U ( 3 ) ), &
             C_ICL % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) &
                                  : C % MOMENTUM_DENSITY_D ( 3 ) ), &
             C_ICR % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) &
                                  : C % MOMENTUM_DENSITY_D ( 3 ) ), &
             C_ICL % Value ( :, C % BARYON_MASS ), &
             C_ICR % Value ( :, C % BARYON_MASS ), &
             C_ICL % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
             C_ICR % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
             C_ICL % Value ( :, C % CONSERVED_ENERGY ), &
             C_ICR % Value ( :, C % CONSERVED_ENERGY ), &
             C_ICL % Value ( :, C % PRESSURE ), &
             C_ICR % Value ( :, C % PRESSURE ), &
             C_ICL % Value ( :, C % CONSERVED_ENTROPY ), &
             C_ICR % Value ( :, C % CONSERVED_ENTROPY ), &
             
             C_IL % Value ( :, C % VELOCITY_U ( 1 ) &
                                 : C % VELOCITY_U ( 3 ) ), &
             C_IR % Value ( :, C % VELOCITY_U ( 1 ) &
                                 : C % VELOCITY_U ( 3 ) ), &
             C_IL % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) &
                                 : C % MOMENTUM_DENSITY_D ( 3 ) ), &
             C_IR % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) &
                                 : C % MOMENTUM_DENSITY_D ( 3 ) ), &
             C_IL % Value ( :, C % BARYON_MASS ), &
             C_IR % Value ( :, C % BARYON_MASS ), &
             C_IL % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
             C_IR % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
             C_IL % Value ( :, C % CONSERVED_ENERGY ), &
             C_IR % Value ( :, C % CONSERVED_ENERGY ), &
             C_IL % Value ( :, C % PRESSURE ), &
             C_IR % Value ( :, C % PRESSURE ), &
             C_IL % Value ( :, C % CONSERVED_ENTROPY ), &
             C_IR % Value ( :, C % CONSERVED_ENTROPY ), &
             SS_I % Value ( :, C % ALPHA_PLUS ), &
             SS_I % Value ( :, C % ALPHA_MINUS ), &
             SS_I % Value ( :, C % ALPHA_CENTER ), &
             M_DD_22, M_DD_33, iD )

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
  
  
  subroutine SetUnits ( VariableUnit, F, Units )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_P_Template ), intent ( in ) :: &
      F
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units

    VariableUnit ( F % INTERNAL_ENERGY )    = Units % EnergyDensity
    VariableUnit ( F % CONSERVED_ENERGY )   = Units % EnergyDensity
    VariableUnit ( F % PRESSURE )           = Units % EnergyDensity
    VariableUnit ( F % TEMPERATURE )        = Units % Temperature
    if ( Units % Temperature % Label /= '' ) then
      VariableUnit ( F % ENTROPY_PER_BARYON ) = UNIT % BOLTZMANN
      VariableUnit ( F % CONSERVED_ENTROPY )  = UNIT % BOLTZMANN &
                                                * Units % NumberDensity
    end if
    VariableUnit ( F % SOUND_SPEED ) = Units % Velocity_U ( 1 )
    VariableUnit ( F % MACH_NUMBER ) = UNIT % IDENTITY

  end subroutine SetUnits


end module Fluid_P__Template
