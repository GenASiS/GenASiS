module Gravitation_NH_C__Form

  !-- Gravitation_NewtonHeader_Chart_Form

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_N  = 1, &
      N_VECTORS_N = 0

  type, public, extends ( Geometry_F_C_Form ) :: Gravitation_NH_C_Form
    integer ( KDI ) :: &
      N_FIELDS_N = N_FIELDS_N, &
      N_VECTORS_N = N_VECTORS_N
    integer ( KDI ) :: &
      POTENTIAL = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    final :: &
      Finalize
  end type Gravitation_NH_C_Form


contains


  subroutine InitializeAllocate_FS &
               ( FSC, C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Gravitation_NH_C_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      oF, &  !-- oField
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a Gravitation_NH_C' 
    
    Name  =  'Gravitation'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  FSC % N_FIELDS_F

    FSC % POTENTIAL  =  oF + 1

    nFields  =  oF  +  FSC % N_FIELDS_N
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + FSC % N_FIELDS_N ) &
      = [ 'Potential' ]

    !-- Units

    !-- Geometry_F

    call FSC % Geometry_F_C_Form % Initialize &
           ( C, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_FS


  impure elemental subroutine Finalize ( GC )

    type ( Gravitation_NH_C_Form ), intent ( inout ) :: &
      GC

  end subroutine Finalize


end module Gravitation_NH_C__Form
