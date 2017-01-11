!-- MessageTemplate provides an abstraction commonly shared by higher-level,
!   concrete type object with specified datatype.

module Message_Template

  use MPI
  use VariableManagement
  use MessagePassingBasics

  implicit none
  private

  type, public, abstract :: MessageTemplate
    integer ( KDI ) :: &
      Rank, &
      Tag, &
      Handle, &
      Error
    !-- FIXME: the following somehow avoids ICE with gfortan/4.10
    !          GCC bug # 61767
    integer ( KDI ), allocatable :: &
      Dummy
    logical ( KDL ) :: &
      Initialized = .false., &
      AllocatedValue = .false.
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      Wait
  end type MessageTemplate

contains


  subroutine InitializeTemplate ( M, C, Tag, Rank, AllocatedOption )
    
    class ( MessageTemplate ), intent ( inout ), target :: &
      M
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      Tag, &
      Rank
    logical ( KDL ), intent ( in ), optional :: &
      AllocatedOption
    
    M % Rank = Rank
    M % Tag  = Tag
    
    M % Initialized    = .true.
    M % AllocatedValue = .true.
    if ( present ( AllocatedOption ) ) M % AllocatedValue = AllocatedOption
    M % Communicator   => C

  end subroutine InitializeTemplate
  
  
  subroutine Wait ( M )
 
    class ( MessageTemplate ), intent ( inout ) :: &
      M

    call MPI_WAIT ( M % Handle, MPI_STATUS_IGNORE, M % Error )

  end subroutine Wait


end module Message_Template
