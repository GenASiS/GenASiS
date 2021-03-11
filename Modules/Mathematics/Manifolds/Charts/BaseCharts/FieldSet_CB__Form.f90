module FieldSet_CB__Form

  !-- FieldSet_ChartBase__Form

  use Basics
  use ManifoldBasics
  use ChartBasics
  use Chart_BH__Form

  implicit none
  private

  type, public, extends ( FieldSet_CH_Form ) :: FieldSet_CB_Form
    integer ( KDI ) :: &
      nValues  = 0, &
      nFields  = 0, &
      nStreams = 0
    character ( LDL ), dimension ( : ), allocatable :: &
      Field
    class ( StorageForm ), allocatable :: &
      FieldSet, &
      FieldSetStream
    type ( Stream_CH_Pointer ), dimension ( : ), allocatable :: &
      Stream
  contains
    procedure, public, pass :: &
      InitializeAllocate
    generic, public :: &
      Initialize => InitializeAllocate
    procedure, public, pass :: &
      AddStream
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateFieldSet
  end type FieldSet_CB_Form

    integer ( KDI ), private, parameter :: &
      MAX_STREAMS = MANIFOLD % MAX_STREAMS

contains


  subroutine InitializeAllocate &
               ( FSC, FSM, C, NameShort, nFields, FieldOption, PinnedOption, &
                 IgnorabilityOption )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    class ( FieldSet_MH_Form ), intent ( in ), target :: &
      FSM
    class ( Chart_BH_Form ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nFields
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    logical ( KDL ), intent ( in ), optional :: &
      PinnedOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FSC % Type == '' ) &
      FSC % Type = 'a FieldSet_CB' 

    FSC % nValues  =  C % nValues
    FSC % nFields  =  nFields

    allocate ( FSC % Field ( nFields ) )
    FSC % Field = ''
    if ( present ( FieldOption ) ) &
      FSC % Field = FieldOption   

    call FSC % Initialize &
           ( FSM, C, NameShort, PinnedOption, IgnorabilityOption )

    call FSC % AllocateFieldSet ( )

    allocate ( FSC % Stream ( MAX_STREAMS ) )

  end subroutine InitializeAllocate


  subroutine AddStream ( FSC, SC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    class ( Stream_CH_Form ), intent ( in ), target :: &
      SC
    
    integer ( KDI ) :: &
      iS

    associate ( nS  =>  FSC % nStreams )

    do iS  =  1, nS
      if ( associated ( FSC % Stream ( iS ) % Pointer, SC ) ) then
        call Show ( 'Stream already added to ' // FSC % Type, &
                    CONSOLE % WARNING )
        call Show ( FSC % Name, 'FieldSet', CONSOLE % WARNING )
        call Show (  SC % Name, 'Stream',   CONSOLE % WARNING )
        return
      end if
    end do !-- iS

    nS  =  nS + 1
    FSC % Stream ( iS ) % Pointer  =>  SC
    call Show ( 'Adding a Stream to ' // trim ( FSC % Type ), &
                FSC % IGNORABILITY + 1 )
    call Show ( FSC % Name, 'FieldSet', FSC % IGNORABILITY + 1 )
    call Show (  SC % Name, 'Stream',   FSC % IGNORABILITY + 1 )

    end associate !-- nS

  end subroutine AddStream


  impure elemental subroutine Finalize ( FSC )

    type ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

    if ( allocated ( FSC % Stream ) ) &
      deallocate ( FSC % Stream )
    if ( allocated ( FSC % FieldSetStream ) ) &
      deallocate ( FSC % FieldSetStream )
    if ( allocated ( FSC % FieldSet ) ) &
      deallocate ( FSC % FieldSet )
    if ( allocated ( FSC % Field ) ) &
      deallocate ( FSC % Field )

  end subroutine Finalize


  subroutine AllocateFieldSet ( FSC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

    if ( allocated ( FSC % FieldSet ) ) &
      return

    call Show ( 'Allocating a FieldSet',    FSC % IGNORABILITY + 1 )
    call Show ( FSC % NameShort, 'Name',    FSC % IGNORABILITY + 1 )
    call Show ( FSC % Field,     'Field',   FSC % IGNORABILITY + 1 )
    call Show ( FSC % nFields,   'nFields', FSC % IGNORABILITY + 1 )
    call Show ( FSC % nValues,   'nValues', FSC % IGNORABILITY + 1 )
    
    allocate ( FSC % FieldSet )
    call FSC % FieldSet % Initialize &
           ( [ FSC % nValues, FSC % nFields ], &
             VariableOption = FSC % Field, NameOption = FSC % NameShort, &
             PinnedOption = FSC % Pinned )

    allocate ( FSC % FieldSetStream )
    call FSC % FieldSetStream % Initialize ( FSC % FieldSet )

  end subroutine AllocateFieldSet


end module FieldSet_CB__Form
