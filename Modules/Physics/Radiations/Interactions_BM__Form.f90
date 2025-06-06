module Interactions_BM__Form

  use Basics
  use Mathematics
  use Fluids
  use Units_R__Form
  use RadiationMoments_BM__Form

 implicit none
 private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_I = 4

  type, public, extends ( FieldSet_BM_Form ) :: Interactions_BM_Form
    integer ( KDI ) :: &
      N_FIELDS_I = N_FIELDS_I
    integer ( KDI ) :: &
      EMISSIVITY_J = 0, &
      EMISSIVITY_H = 0, &
      OPACITY_J    = 0, &
      OPACITY_H    = 0
    class ( Fluid_P_Form ), pointer :: &
      Fluid => null ( )
    class ( RadiationMoments_BM_Form ), pointer :: &
      Radiation    => null ( ), &
      RadiationBar => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    generic, public :: &
      Initialize => InitializeAllocate_I
    procedure, public, pass ( I ) :: &
      SetStream
    procedure, public, pass :: &
      ComputeAll
    procedure, public, pass :: &
      ComputeSingle
    generic, public :: &
      Compute  =>  ComputeAll, &
                   ComputeSingle
    final :: &
      Finalize
  end type Interactions_BM_Form


contains


  subroutine InitializeAllocate_I &
               ( I, R, Units_R, F, FieldOption, NameOption, UnitOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Interactions_BM_Form ), intent ( inout ) :: &
      I
    class ( RadiationMoments_BM_Form ), intent ( inout ), target :: &
      R
    class ( Units_R_Form ), dimension ( : ), intent ( in ) :: &
      Units_R
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F
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
      iC, &
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

    I % Fluid      =>  F
    I % Radiation  =>  R

    !-- Field indices

    I % EMISSIVITY_J  =  1
    I % EMISSIVITY_H  =  2
    I % OPACITY_J     =  3
    I % OPACITY_H     =  4

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
      = [ 'Emissivity_J', &
          'Emissivity_H', &
          'Opacity_J   ', &
          'Opacity_H   ' ]
          
    !-- Units

    associate ( nC  =>  F % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC
      FieldUnit ( I % EMISSIVITY_J, iC ) &
        =  Units_R ( iC ) % EnergyDensity  &
           *  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
      FieldUnit ( I % EMISSIVITY_H, iC ) &
        =  Units_R ( iC ) % EnergyDensity  &
           *  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
      FieldUnit ( I % OPACITY_J, iC ) &
        =  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
      FieldUnit ( I % OPACITY_H, iC ) &
        =  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
    end do !-- iC

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

    ! !-- Fluid SplitSource

    ! allocate ( I % Fluid % SplitSource )
    ! call I % Fluid % SplitSource % Initialize &
    !        ( F % Atlas, &
    !          FieldOption = F % Balanced, &
    !          NameOption = 'Fluid_SplitSource', &
    !          DeviceMemoryOption = F % DeviceMemory, &
    !          PinnedMemoryOption = F % PinnedMemory, &
    !          DevicesCommunicateOption = F % DevicesCommunicate, &
    !          nFieldsOption = F % nBalanced, &
    !          IgnorabilityOption = F % IGNORABILITY )

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
                   I % OPACITY_J, &
                   I % OPACITY_H ] )

  end subroutine SetStream


  subroutine ComputeAll ( I )

    class ( Interactions_BM_Form ), intent ( inout ) :: &
      I

    call Show ( 'Must be replaced by extension', CONSOLE % ERROR )
    call Show ( 'Interactions_BM__Form', 'module', CONSOLE % ERROR )
    call Show ( 'ComputeAll', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ComputeAll


  subroutine ComputeSingle ( I, iC, iV )

    class ( Interactions_BM_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      iC, &
      iV

    call Show ( 'Must be replaced by extension', CONSOLE % ERROR )
    call Show ( 'Interactions_BM__Form', 'module', CONSOLE % ERROR )
    call Show ( 'ComputeSingle', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ComputeSingle


  impure elemental subroutine Finalize ( I )

    type ( Interactions_BM_Form ), intent ( inout ) :: &
      I

    nullify ( I % RadiationBar )
    nullify ( I % Radiation )
    nullify ( I % Fluid )

  end subroutine Finalize


end module Interactions_BM__Form
