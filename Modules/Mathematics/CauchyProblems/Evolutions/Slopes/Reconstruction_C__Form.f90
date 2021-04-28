module Reconstruction_C__Form

  !-- Reconstruction_Chart_Form

  use Basics
  use Manifolds
  use Fields

  implicit none
  private

  type, public :: Reconstruction_C_Form
    integer ( KDI ) :: &
      IGNORABILITY, &
      Order
    logical ( KDL ) :: &
      Streamed
    character ( LDL ) :: &
      Name
    class ( FieldSet_C_Form ), allocatable :: &
      Output_IL_C, Output_IR_C
    class ( FieldSet_C_Form ), pointer :: &
      FieldSet_C => null ( )
    class ( Geometry_F_C_Form ), pointer :: &
      Geometry_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_R
    generic, public :: &
      Initialize => InitializeAllocate_R
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      Show => Show_RC
    final :: &
      Finalize
  end type Reconstruction_C_Form


contains


  subroutine InitializeAllocate_R &
               ( RC, GC, FSC, NameOption, StreamedOption, OrderOption )

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
      OrderOption

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    RC % IGNORABILITY  =  FSC % IGNORABILITY

    RC % Name  =  'Reconstruction_' // trim ( FSC % Name )
    if ( present ( NameOption ) ) &
      RC % Name  =  trim ( NameOption ) // '_' // trim ( FSC % Name )

    call Show ( 'Initializing a Reconstruction_C', RC % IGNORABILITY )
    call Show ( RC % Name, 'Name', RC % IGNORABILITY )
   
    associate ( nF  =>  FSC % nFields )
 
    allocate ( Field ( nF ) )
    allocate ( Unit ( nF ) )
    do iS  =  1,  nF
      iF  =  FSC % iaSelected ( iS )
      Field ( iF )  =  FSC % Field ( iF )
      Unit  ( iF )  =  FSC % Unit  ( iF )
    end do !-- iS

    allocate ( RC % Output_IL_C )
    allocate ( RC % Output_IR_C )
    call RC % Output_IL_C % Initialize &
           ( FSC % Chart, &
             FieldOption = Field, &
             NameOption = trim ( RC % Name ) // '_IL', &
             DeviceMemoryOption = FSC % Storage_FSC % DeviceMemory, &
             PinnedMemoryOption = FSC % Storage_FSC % PinnedMemory, &
             DevicesCommunicateOption = FSC % GhostExchange_FSC &
                                          % DevicesCommunicate, &
             UnitOption = Unit, &
             nFieldsOption = nF, &
             IgnorabilityOption = FSC % Ignorability )
    call RC % Output_IR_C % Initialize &
           ( FSC % Chart, &
             FieldOption = Field, &
             NameOption = trim ( RC % Name ) // '_IR', &
             DeviceMemoryOption = FSC % Storage_FSC % DeviceMemory, &
             PinnedMemoryOption = FSC % Storage_FSC % PinnedMemory, &
             DevicesCommunicateOption = FSC % GhostExchange_FSC &
                                          % DevicesCommunicate, &
             UnitOption = Unit, &
             nFieldsOption = nF, &
             IgnorabilityOption = FSC % Ignorability )

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


  subroutine Compute ( RC, iD )

    class ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions

    real ( KDR ), dimension ( :, :, : ), pointer :: &
      X, &
      dX_L, dX_R, &
      F_IL, &
      F_IR

    select type ( C  =>  RC % FieldSet_C % Chart )
    class is ( Chart_GS_Form )

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Reconstruction_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine Compute


  subroutine Show_RC ( RC )

    class ( Reconstruction_C_Form ), intent ( in ) :: &
      RC

    call Show ( 'Reconstruction_C Parameters', RC % IGNORABILITY )

    call Show ( RC % Name, 'Name',  RC % IGNORABILITY )
    call Show ( RC % Order, 'Order', RC % IGNORABILITY )
    call Show ( RC % Streamed, 'Streamed', RC % IGNORABILITY )
    call RC % Output_IL_C % Show ( )
    call RC % Output_IR_C % Show ( )

  end subroutine Show_RC


  impure elemental subroutine Finalize ( RC )

    type ( Reconstruction_C_Form ), intent ( inout ) :: &
      RC

    nullify ( RC % Geometry_C )
    nullify ( RC % FieldSet_C )

    if ( allocated ( RC % Output_IR_C ) ) &
      deallocate ( RC % Output_IR_C )
    if ( allocated ( RC % Output_IL_C ) ) &
      deallocate ( RC % Output_IL_C )

    call Show ( 'Finalizing a Reconstruction_C', RC % IGNORABILITY )
    call Show ( RC % Name, 'Name', RC % IGNORABILITY )
   
  end subroutine Finalize

  
end module Reconstruction_C__Form
