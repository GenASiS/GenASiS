module Geometry_F_MH__Form

  !-- Geometry_Flat_ManifoldHeader_Form

  use Basics
  use Manifold_H__Form
  use FieldSet_MH__Form

  implicit none
  private

  type, public, extends ( FieldSet_MH_Form ) :: Geometry_F_MH_Form
  contains
    procedure, public, pass :: &
      Initialize_F
    final :: &
      Finalize
  end type Geometry_F_MH_Form

contains


  subroutine Initialize_F &
               ( GM, M, NameOption, DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption )

    class ( Geometry_F_MH_Form ), intent ( inout ) :: &
      GM
    class ( Manifold_H_Form ), intent ( inout ), target :: &
      M
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption

    character ( LDL ) :: &
      Name

    if ( GM % Type == '' ) &
      GM % Type = 'a Geometry_F_M'

    Name  =  'Geometry'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call GM % FieldSet_MH_Form % Initialize &
           ( M, Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption )

  end subroutine Initialize_F


  impure elemental subroutine Finalize ( GM )

    type ( Geometry_F_MH_Form ), intent ( inout ) :: &
      GM

  end subroutine Finalize


end module Geometry_F_MH__Form
