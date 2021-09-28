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
!    logical ( KDL ) :: &
!      UseInitialTemperature, &
!      UseEntropy
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
  !   procedure, public, pass :: &
  !     SetPrimitiveConservedTemplate_P
  !   procedure, public, pass :: &
  !     ComputeFromInitial
  !   procedure ( CFT ), public, pass ( C ), deferred :: &
  !     ComputeFromTemperature
  !   procedure, public, pass ( C ) :: &
  !     ComputeFluxes
  !   procedure, public, pass ( C ) :: &
  !     ComputeCenterStates
  !   procedure, public, pass ( C ) :: &
  !     ComputeRawFluxesTemplate_P
  !   procedure, public, pass ( C ) :: &
  !     ComputeCenterStatesTemplate_P
  !   procedure, public, nopass :: &
  !     Compute_G_G_Kernel
  !   procedure, public, nopass :: &
  !     Compute_DS_G_Kernel
  !   procedure, public, nopass :: &
  !     Compute_N_V_E_G_Kernel
  !   procedure, public, nopass :: &
  !     Compute_SB_G_Kernel
  !   procedure, public, nopass :: &
  !     Compute_FE_P_G_Kernel
    final :: &
      Finalize
  end type Fluid_P_Form


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
      if ( Units_F ( iC ) % Temperature % Label /= '' ) then
        FieldUnit ( F % ENTROPY_PER_BARYON, iC ) &
          =  UNIT % BOLTZMANN
      end if
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

  end subroutine InitializeAllocate_F

  
  impure elemental subroutine Finalize ( F )

    type ( Fluid_P_Form ), intent ( inout ) :: &
      F

  end subroutine Finalize


end module Fluid_P__Form
