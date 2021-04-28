module Reconstruction_C__Form

  !-- Reconstruction_Chart_Form

  use Basics
  use Manifolds
  use Fields

  implicit none
  private

  type, public, extends ( FieldSet_C_Form ) :: Reconstruction_C_Form
    integer ( KDI ) :: &
      Order
    logical ( KDL ) :: &
      Streamed
    class ( FieldSet_C_Form ), pointer :: &
      FieldSet_C => null ( )
    class ( Geometry_F_C_Form ), pointer :: &
      Geometry_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_R
    generic, public :: &
      Initialize => InitializeAllocate_R
    procedure, private, pass :: &
      Show_FSC
    final :: &
      Finalize
  end type Reconstruction_C_Form


contains


  subroutine InitializeAllocate_R &
               ( RC, GC, FSC, NameOption, StreamedOption, OrderOption, &
                 IgnorabilityOption )

    class ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GC
    class ( FieldSet_C_Form ), intent ( in ), target :: &
      FSC
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      StreamedOption
    integer ( KDI ), intent ( in ), optional :: &
      OrderOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF, &  !-- iField
      Ignorability
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    Ignorability  =  FSC % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      Ignorability  =  IgnorabilityOption

    RC % Type  =  'a Reconstruction_C'

    Name  =  'Reconstruction_' // trim ( FSC % Name )
    if ( present ( NameOption ) ) &
      Name  =  trim ( NameOption ) // '_' // trim ( FSC % Name )

    associate ( nF  =>  FSC % nFields )
 
    allocate ( Field ( nF ) )
    allocate ( Unit ( nF ) )
    do iS  =  1,  nF
      iF  =  FSC % iaSelected ( iS )
      Field ( iF )  =  FSC % Field ( iF )
      Unit  ( iF )  =  FSC % Unit  ( iF )
    end do !-- iS

    call RC % FieldSet_C_Form % Initialize &
           ( FSC % Chart, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = FSC % Storage_FSC % DeviceMemory, &
             PinnedMemoryOption = FSC % Storage_FSC % PinnedMemory, &
             DevicesCommunicateOption = FSC % GhostExchange_FSC &
                                          % DevicesCommunicate, &
             UnitOption = Unit, &
             nFieldsOption = nF, &
             IgnorabilityOption = Ignorability )

    end associate !-- nF

    RC % FieldSet_C  =>  FSC
    RC % Geometry_C  =>   GC

    RC % Order  =  0
    if ( present ( OrderOption ) ) &
      RC % Order  =  OrderOption

    RC % Streamed  =  .false.
    if ( present ( StreamedOption ) ) &
      RC % Streamed  =  StreamedOption

  end subroutine InitializeAllocate_R


  subroutine Show_FSC ( FSC )

    class ( Reconstruction_C_Form ), intent ( in ) :: &
      FSC

    call FSC % FieldSet_C_Form % Show ( )
    call Show ( FSC % Order, 'Order', FSC % IGNORABILITY )
    call Show ( FSC % Streamed, 'Streamed', FSC % IGNORABILITY )

  end subroutine Show_FSC


  impure elemental subroutine Finalize ( RC )

    type ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC

    nullify ( RC % Geometry_C )
    nullify ( RC % FieldSet_C )

  end subroutine Finalize

  
end module Reconstruction_C__Form
