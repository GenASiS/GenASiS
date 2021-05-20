module Gravitation_G_C__Form

  !-- Gravitation_Galileo_Chart_Form

  use Basics
  use Mathematics

  implicit none
  private

  type, public, extends ( Geometry_F_C_Form ) :: Gravitation_G_C_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
  end type Gravitation_G_C_Form


contains


  subroutine InitializeAllocate_FS &
               ( FSC, C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Gravitation_G_C_Form ), intent ( inout ) :: &
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

    character ( LDL ) :: &
      Name

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a Gravitation_G_C' 
    
    Name  =  'Gravitation'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call FSC % Geometry_F_C_Form % Initialize &
           ( C, &
             FieldOption = FieldOption, &
             VectorOption = VectorOption, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_FS


end module Gravitation_G_C__Form
