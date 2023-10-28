module Interactions_BM__Form

  use Basics
  use Mathematics
  use Fluids
  use Units_R__Form

 implicit none
 private

  integer ( KDI ), private, parameter :: &
    N_FIELDS_I = 8

  type, public, extends ( FieldSet_BM_Form ) :: Interactions_BM_Form
    integer ( KDI ) :: &
      N_FIELDS_I = N_FIELDS_I
    integer ( KDI ) :: &
      EMISSIVITY_J  = 0, &
      EMISSIVITY_H  = 0, &
      EMISSIVITY_N  = 0, &
      OPACITY_J     = 0, &
      OPACITY_H     = 0, &
      OPACITY_N     = 0, &
      EQUILIBRIUM_J = 0, &
      EQUILIBRIUM_N = 0
    class ( Fluid_P_Form ), pointer :: &
      Fluid => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    generic, public :: &
      Initialize => InitializeAllocate_I
    procedure, public, pass ( I ) :: &
      SetStream
    procedure, public, pass :: &
      Compute
    ! procedure, private, pass ( I ) :: &
    !   ComputeEquilibrium_T
    ! procedure, private, pass ( I ) :: &
    !   ComputeEquilibrium_T_Eta
    final :: &
      Finalize
    procedure, public, nopass :: &
      Compute_J_Eq_Ph_G_Kernel
  end type Interactions_BM_Form

  interface

    module subroutine Compute_J_Eq_Ph_G_Kernel ( J_Eq, T, UseDeviceOption )
      !-- Compute_J_Eq_Photons_Grey_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        J_Eq
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        T
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_J_Eq_Ph_G_Kernel

  end interface


contains


  subroutine InitializeAllocate_I &
               ( I, F, Units_R, FieldOption, NameOption, UnitOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Interactions_BM_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F
    class ( Units_R_Form ), dimension ( : ), intent ( in ) :: &
      Units_R
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      nFields
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( I % Type  ==  '' ) &
      I % Type  =  'an Interactions' 
    
    Name  =  'Interactions'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    I % Fluid  =>  F

    !-- Field indices

    I % EMISSIVITY_J   =  1
    I % EMISSIVITY_H   =  2
    I % EMISSIVITY_N   =  3
    I % OPACITY_J      =  4
    I % OPACITY_H      =  5
    I % OPACITY_N      =  6
    I % EQUILIBRIUM_J  =  7
    I % EQUILIBRIUM_N  =  8

    nFields  =  I % N_FIELDS_I
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( 1 : I % N_FIELDS_I ) &
      = [ 'Emissivity_J ', &
          'Emissivity_H ', &
          'Emissivity_N ', &
          'Opacity_J    ', &
          'Opacity_H    ', &
          'Opacity_N    ', &
          'Equilibrium_J', &
          'Equilibrium_N' ]
          
    !-- Units

    associate ( nC  =>  F % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    ! do iC  =  1, nC
      ! FieldUnit ( I % EMISSIVITY_J ) &
      !   =  Units % EnergyDensity  *  Units % Length ** (-1)
      ! FieldUnit ( I % EMISSIVITY_H ) &
      !   =  Units % EnergyDensity  *  Units % Length ** (-1)
      ! FieldUnit ( I % EMISSIVITY_N ) &
      !   =  Units % NumberDensity  *  Units % Length ** (-1)
      ! FieldUnit ( I % OPACITY_J ) &
      !   =  Units % Length ** (-1)
      ! FieldUnit ( I % OPACITY_H ) &
      !   =  Units % Length ** (-1)
      ! FieldUnit ( I % OPACITY_N ) &
      !   =  Units % Length ** (-1)
      ! FieldUnit ( I % EQUILIBRIUM_J ) &
      !   =  Units % EnergyDensity
      ! FieldUnit ( I % EQUILIBRIUM_N ) &
      !   =  Units % NumberDensity
    ! end do !-- iC

    end associate !-- nC

    !-- FieldSet

    call I % FieldSet_BM_Form % Initialize &
           ( F % Atlas, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = F % DeviceMemory, &
             PinnedMemoryOption = F % PinnedMemory, &
             DevicesCommunicateOption = F % DevicesCommunicate, &
             UnitOption = FieldUnit, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_I


  subroutine SetStream ( S, I )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( Interactions_BM_Form ), intent ( in ) :: &
      I

    call S % AddFieldSet &
           ( I, &
             iaSelectedOption &
               = [ I % EMISSIVITY_J, &
                   I % EMISSIVITY_H, &
                   I % EMISSIVITY_N, &
                   I % OPACITY_J, &
                   I % OPACITY_H, &
                   I % OPACITY_N, &
                   I % EQUILIBRIUM_J, &
                   I % EQUILIBRIUM_N ] )

  end subroutine SetStream


  subroutine Compute ( I )

    class ( Interactions_BM_Form ), intent ( inout ) :: &
      I

    !-- To be filled in by extension

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_BM_Form ), intent ( inout ) :: &
      I

    nullify ( I % Fluid )

  end subroutine Finalize


end module Interactions_BM__Form
