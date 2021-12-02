module Gravitation_G__Form

  !-- Gravitation_Galileo__Form

  use Basics
  use Mathematics

  implicit none
  private

  type, public, extends ( Geometry_F_Form ) :: Gravitation_G_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    final :: &
      Finalize
  end type Gravitation_G_Form


contains


  subroutine InitializeAllocate_FS &
               ( FS, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, AssociateFieldsOption, &
                 UnitOption, VectorIndicesOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( Gravitation_G_Form ), intent ( inout ), target :: &
      FS
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption, &
      AssociateFieldsOption
    type ( MeasuredValueForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a Gravitation_G' 
    
    Name  =  'Gravitation'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call FS % Geometry_F_Form % Initialize &
           ( A, &
             FieldOption = FieldOption, &
             VectorOption = VectorOption, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             AssociateFieldsOption = AssociateFieldsOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_FS


  impure elemental subroutine Finalize ( G )

    type ( Gravitation_G_Form ), intent ( inout ) :: &
      G

  end subroutine Finalize


end module Gravitation_G__Form
