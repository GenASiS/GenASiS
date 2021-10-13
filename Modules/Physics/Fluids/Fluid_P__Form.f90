module Fluid_P__Form

  !-- Fluid_Perfect__Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_D__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_P = 1, &
      N_BALANCED_P = 1, &
      N_FIELDS_P    = 7, &
      N_VECTORS_P   = 0

  type, public, extends ( Fluid_D_Form ) :: Fluid_P_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_P      = N_PRIMITIVE_P, &
      N_BALANCED_P       = N_BALANCED_P, &
      N_FIELDS_P         = N_FIELDS_P, &
      N_VECTORS_P        = N_VECTORS_P, &
      ENERGY_DENSITY_C   = 0, &
      ENERGY_DENSITY_B   = 0, &
      PRESSURE           = 0, &
      TEMPERATURE        = 0, &
      ENTROPY_PER_BARYON = 0, &
      SOUND_SPEED        = 0, &
      MACH_NUMBER        = 0
   logical ( KDL ) :: &
     UseInitialTemperature
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    procedure, public, pass :: &
      SetUseInitialTemperature
    procedure, public, pass :: &
      ComputeFromInitial
    procedure, public, pass :: &
      ComputeFromTemperature
    procedure, public, pass ( CS ) :: &
      ComputeFluxes
  !   procedure, public, pass ( C ) :: &
  !     ComputeCenterStates
  !   procedure, public, pass ( C ) :: &
  !     ComputeCenterStatesTemplate_P
  !   procedure, public, nopass :: &
  !     Compute_SB_G_Kernel
  !   procedure, public, nopass :: &
  !     Compute_FE_P_G_Kernel
    final :: &
      Finalize
    procedure, public, nopass :: &
      Compute_D_S_G_G_Kernel
    procedure, public, nopass :: &
      Compute_N_V_E_G_Kernel
    procedure, public, nopass :: &
      ComputeFluxes_G_Kernel
  end type Fluid_P_Form


  interface

    module subroutine Compute_D_S_G_G_Kernel & 	 	 
             ( N, V_1, V_2, V_3, M, E, CS, M_DD_11, M_DD_22, M_DD_33, N_Min, &
               D, S_1, S_2, S_3, G, MN, UseDeviceOption )
      !-- Compute_DensityB_Momentum_EnergyB_Galileo_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
        N, & 	 	 
        V_1, V_2, V_3
      real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
        M, &
        E, &
        CS, &
        M_DD_11, M_DD_22, M_DD_33
      real ( KDR ), intent ( in ) :: &
        N_Min
      real ( KDR ), dimension ( : ), intent ( out ) :: & 	 	 
        D, & 	 	 
        S_1, S_2, S_3, &
        G, &
        MN
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_D_S_G_G_Kernel 	 	 

    module subroutine Compute_N_V_E_G_Kernel &
             ( D, S_1, S_2, S_3, G, M, M_UU_11, M_UU_22, M_UU_33, N_Min, &
               N, V_1, V_2, V_3, E, UseDeviceOption )
      !-- Compute_DensityC_Velocity_EnergyC_Galileo_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        D, &
        S_1, S_2, S_3, &
        G
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        M_UU_11, M_UU_22, M_UU_33
      real ( KDR ), intent ( in ) :: &
        N_Min
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        N, &
        V_1, V_2, V_3, &
        E
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_N_V_E_G_Kernel

    module subroutine ComputeFluxes_G_Kernel &
             ( D, S_1, S_2, S_3, G, P, V_Dim, iDim, &
               F_D, F_S_1, F_S_2, F_S_3, F_G, UseDeviceOption )
      !-- ComputeFluxes_Galileo_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        D, &
        S_1, S_2, S_3, &
        G, &
        P, &
        V_Dim
      integer ( KDI ), intent ( in ) :: &
        iDim
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        F_D, &
        F_S_1, F_S_2, F_S_3, &
        F_G
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeFluxes_G_Kernel

  end interface


contains


  subroutine InitializeAllocate_F &
               ( F, G, Units_F, FieldOption, VectorOption, NameOption, &
                 UnitOption, VectorIndicesOption, iaPrimitiveOption, &
                 iaBalancedOption, nFieldsOption, IgnorabilityOption )

    class ( Fluid_P_Form ), intent ( inout ) :: &
      F
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( Units_F_Form ), dimension ( : ), intent ( in ) :: &
      Units_F
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaPrimitiveOption, &
      iaBalancedOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iP, &  !-- iPrimitive
      iB, &  !-- iBalanced
      iC, &  !-- iChart
      oF, &  !-- oField
      oP, &  !-- oPrimitive
      oB, &  !-- oBalanced
      nFields, &
      nPrimitive, &
      nBalanced
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaBalanced
    type ( MeasuredValueForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( F % Type  ==  '' ) &
      F % Type  =  'a Fluid_P' 
    
    !-- Field indices

    oF  =  F % N_FIELDS_CS  +  F % N_FIELDS_D

    F % ENERGY_DENSITY_C    =  oF + 1
    F % ENERGY_DENSITY_B    =  oF + 2
    F % PRESSURE            =  oF + 3
    F % TEMPERATURE         =  oF + 4
    F % ENTROPY_PER_BARYON  =  oF + 5
    F % SOUND_SPEED         =  oF + 6
    F % MACH_NUMBER         =  oF + 7

    nFields  =  oF  +  F % N_FIELDS_P
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + F % N_FIELDS_P ) &
      = [ 'EnergyDensity_C ', &
          'EnergyDensity_B ', &
          'Pressure        ', &
          'Temperature     ', &
          'EntropyPerBaryon', &
          'SoundSpeed      ', &
          'MachNumber      ' ]

    !-- Units

    associate ( nC  =>  G % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC
      FieldUnit ( F % ENERGY_DENSITY_C, iC ) &
        =  Units_F ( iC ) % EnergyDensity
      FieldUnit ( F % ENERGY_DENSITY_B, iC ) &
        =  Units_F ( iC ) % SqrtDet_M  *  Units_F ( iC ) % EnergyDensity
      FieldUnit ( F % PRESSURE, iC ) &
        =  Units_F ( iC ) % EnergyDensity
      FieldUnit ( F % TEMPERATURE, iC ) &
        =  Units_F ( iC ) % Temperature
      FieldUnit ( F % ENTROPY_PER_BARYON, iC ) &
        =  Units_F ( iC ) % Energy  /  Units_F ( iC ) % Temperature
      FieldUnit ( F % SOUND_SPEED, iC ) &
        =  Units_F ( iC ) % Velocity_U ( 1 )
      FieldUnit ( F % MACH_NUMBER, iC ) &
        =  UNIT % IDENTITY
    end do !-- iC

    end associate !-- nC

    !-- Vector indices: no additional vectors

    !-- Vector names: no additional vectors

    !-- Primitive fields

    oP  =  F % N_PRIMITIVE_CS  +  F % N_PRIMITIVE_D

    if ( present ( iaPrimitiveOption ) ) then
      nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( iaPrimitive, source = iaPrimitiveOption )
    else
      nPrimitive  =  oP  +  F % N_PRIMITIVE_P
      allocate ( iaPrimitive ( nPrimitive ) )
    end if !-- iaPrimitiveOption

    iaPrimitive ( oP  +  1 : oP  +  F % N_PRIMITIVE_P )  &
      =  [ F % ENERGY_DENSITY_C ]

    !-- Balanced fields

    oB  =  F % N_BALANCED_CS  +  F % N_BALANCED_D

    if ( present ( iaBalancedOption ) ) then
      nBalanced  =  size ( iaBalancedOption )
      allocate ( iaBalanced, source = iaBalancedOption )
    else
      nBalanced  =  oB  +  F % N_BALANCED_P
      allocate ( iaBalanced ( nBalanced ) )
    end if !-- iaPrimitiveOption

    iaBalanced ( oB  +  1 : oB  +  F % N_BALANCED_P )  &
      =  [ F % ENERGY_DENSITY_B ]

    !-- Fluid_D

    call F % Fluid_D_Form % Initialize &
           ( G, Units_F, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = NameOption, &
             UnitOption = FieldUnit, &
             VectorIndicesOption = VectorIndicesOption, &
             iaPrimitiveOption = iaPrimitive, &
             iaBalancedOption = iaBalanced, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    !-- Parameters

    F % UseInitialTemperature  =  .false.

  end subroutine InitializeAllocate_F


  subroutine SetUseInitialTemperature ( F, UseInitialTemperature )

    class ( Fluid_P_Form ), intent ( inout ) :: &
      F
    logical ( KDL ), intent ( in ) :: &
      UseInitialTemperature

    F % UseInitialTemperature  =  UseInitialTemperature

    call Show ( 'Setting UseInitialTemperatre of a Fluid_P', &
                F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( F % UseInitialTemperature, 'UseInitialTemperature', &
                F % IGNORABILITY + 1 )

  end subroutine SetUseInitialTemperature


  subroutine ComputeFromInitial ( CS )

    class ( Fluid_P_Form ), intent ( inout ) :: &
      CS

    call Show ( 'ComputeFromInitial', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    if ( CS % UseInitialTemperature ) then
      call CS % ComputeFromTemperature ( )
    else
      call CS % ComputeFromPrimitive ( CS )
    end if

    select type ( G  =>  CS % Geometry )
      class is ( Gravitation_N_H_Form )
    !-- FIXME Constant_G
    call G % Solve &
           ( CS, Constant_G = 1.0_KDR, &
             iBaryonMass = CS % BARYON_MASS, &
             iBaryonDensity = CS % BARYON_DENSITY_B )
    end select !-- G

  end subroutine ComputeFromInitial


  subroutine ComputeFromTemperature ( F )

    class ( Fluid_P_Form ), intent ( inout ) :: &
      F

    call Show ( 'ComputeFromTemperature should be overridden', &
                CONSOLE % WARNING)
    call Show ( F % Name, 'Fluid', CONSOLE % WARNING )

  end subroutine ComputeFromTemperature


  subroutine ComputeFluxes ( FS, CS, FS_CS, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    class ( Fluid_P_Form ), intent ( in ) :: &
      CS
    class ( FieldSetForm ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
    
    integer ( KDI ) :: &
      iDensity, &
      iEnergy
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum

    call Search &
           ( CS % iaBalanced, CS % BARYON_DENSITY_B, iDensity )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D_1, iMomentum ( 1 ) )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D_2, iMomentum ( 2 ) )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D_3, iMomentum ( 3 ) )
    call Search &
           ( CS % iaBalanced, CS % ENERGY_DENSITY_B, iEnergy )

    associate &
      ( FSV  =>  FS    % Storage ( iC ) % Value, &
        CSV  =>  FS_CS % Storage ( iC ) % Value )
    associate &
      ( F_D      =>  FSV ( :, iDensity ), &
        F_S_1    =>  FSV ( :, iMomentum ( 1 ) ), &
        F_S_2    =>  FSV ( :, iMomentum ( 2 ) ), &
        F_S_3    =>  FSV ( :, iMomentum ( 3 ) ), &
        F_G      =>  FSV ( :, iEnergy ), &
          D      =>  CSV ( :, CS % BARYON_DENSITY_B ), &
          S_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_D_1 ), &
          S_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
          S_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ), &
          G      =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
          P      =>  CSV ( :, CS % PRESSURE ), &
          V_Dim  =>  CSV ( :, CS % VELOCITY_U ( iD ) ) )
 
    call ComputeFluxes_G_Kernel &
           ( D, S_1, S_2, S_3, G, P, V_Dim, iD, F_D, F_S_1, F_S_2, F_S_3, F_G, &
             UseDeviceOption = CS % DeviceMemory )
  
    end associate !-- F_D, etc.
    end associate !-- FSV, etc.

  end subroutine ComputeFluxes


  impure elemental subroutine Finalize ( F )

    type ( Fluid_P_Form ), intent ( inout ) :: &
      F

  end subroutine Finalize


end module Fluid_P__Form
