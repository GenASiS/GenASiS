module Geometry_F_MH__Form

  !-- Geometry_Flat_ManifoldHeader_Form

  use Basics
  use Manifold_H__Form
  use FieldSet_MH__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_FLAT = 19

  type, public, extends ( FieldSet_MH_Form ) :: Geometry_F_MH_Form
    integer ( KDI ) :: &
      N_FIELDS_FLAT = N_FIELDS_FLAT
  contains
    procedure, private, pass :: &
      Initialize_F
    generic, public :: &
      Initialize => Initialize_F
    final :: &
      Finalize
  end type Geometry_F_MH_Form

contains


  subroutine Initialize_F &
               ( GM, M, FieldOption, NameOption, DeviceMemoryOption, &
                 PinnedMemoryOption, DevicesCommunicateOption, UnitOption, &
                 nFieldsOption )

    class ( Geometry_F_MH_Form ), intent ( inout ) :: &
      GM
    class ( Manifold_H_Form ), intent ( inout ), target :: &
      M
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      nFields
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( GM % Type == '' ) &
      GM % Type = 'a Geometry_F'

    Name  =  'Geometry'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    nFields  =  GM % N_FIELDS_FLAT
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    if ( present ( FieldOption ) ) then
      Field  =  FieldOption
    else
      allocate ( Field ( nFields ) )
    end if

    Field ( 1 : GM % N_FIELDS_FLAT ) &
      = [ 'Edge_I_U_1    ', &
          'Edge_I_U_2    ', &
          'Edge_I_U_3    ', &
          'Width_U_1     ', &
          'Width_U_2     ', &
          'Width_U_3     ', &
          'Center_U_1    ', &
          'Center_U_2    ', &
          'Center_U_3    ', &
          'Area_I_D_1    ', &
          'Area_I_D_2    ', &
          'Area_I_D_3    ', &
          'Volume        ', &
          'Metric_F_DD_11', &
          'Metric_F_DD_22', &
          'Metric_F_DD_33', &
          'Metric_F_UU_11', &
          'Metric_F_UU_22', &
          'Metric_F_UU_33' ]

    if ( present ( UnitOption ) ) then
      Unit  =  UnitOption
    else
      allocate ( Unit ( nFields ) )
    end if

    call GM % FieldSet_MH_Form % Initialize &
           ( M, Name, nFields, FieldOption = Field, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption )

  end subroutine Initialize_F


  subroutine Finalize ( GM )

    type ( Geometry_F_MH_Form ), intent ( inout ) :: &
      GM

  end subroutine Finalize


end module Geometry_F_MH__Form
