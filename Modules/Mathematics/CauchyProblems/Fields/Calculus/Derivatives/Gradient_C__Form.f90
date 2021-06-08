module Gradient_C__Form

  !-- Gradient_Chart_Form

  use Basics
  use Manifolds
  use FieldSets
  use Streams
  use Geometries

  implicit none
  private

  type, public, extends ( FieldSet_C_Form ) :: Gradient_C_Form
    class ( FieldSet_C_Form ), pointer :: &
      FieldSet_C  => null ( )
    class ( Geometry_F_C_Form ), pointer :: &
      Geometry_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_G
    generic, public :: &
      Initialize => InitializeAllocate_G
    final :: &
      Finalize
  end type Gradient_C_Form


contains


  subroutine InitializeAllocate_G ( GC, GyC, FSC, NameOption )

    class ( Gradient_C_Form ), intent ( inout ) :: &
      GC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GyC
    class ( FieldSet_C_Form ), intent ( in ), target :: &
      FSC
    character ( * ), intent ( in ), optional :: &
      NameOption
 
    integer ( KDI ) :: &
      iS, &    !-- iSelected
      iF, &    !-- iField
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( GC % Type  ==  '' ) &
      GC % Type  =  'a Gradient_C' 
    
    Name  =  'Grad_' // trim ( FSC % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    GC % FieldSet_C  =>  FSC
    GC % Geometry_C  =>  GyC

    associate &
      ( DeviceMemory  =>  FSC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  FSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  FSC % GhostExchange_FSC % DevicesCommunicate ) 

    nFields  =  FSC % nFields
    
    allocate ( Field ( nFields ) )
    do iS  =  1, nFields
      iF  =  FSC % iaSelected ( iS )
      Field ( iS )  =  FSC % Field ( iF )
    end do !-- iS

    call GC % FieldSet_C_Form % Initialize &
           ( FSC % Chart, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nFields, &
             IgnorabilityOption = FSC % IGNORABILITY )

    end associate !-- DeviceMemory, etc.

  end subroutine InitializeAllocate_G


  impure elemental subroutine Finalize ( GC )

    type ( Gradient_C_Form ), intent ( inout ) :: &
      GC

    nullify ( GC % Geometry_C )
    nullify ( GC % FieldSet_C )

  end subroutine Finalize



end module Gradient_C__Form
