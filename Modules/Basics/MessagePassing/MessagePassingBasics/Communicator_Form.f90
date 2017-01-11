!-- CommunicatorForm is a class to abstract the notion of communicator that
!   connect a group of processes in a distributed-memory parallel program. 

module Communicator_Form
  
  use MPI
  use VariableManagement
  use Display

  implicit none
  private

  type, public :: CommunicatorForm
    integer ( KDI ) :: &
      Handle, &
      Size, &
      Rank, &
      Error
    integer ( KDI ), dimension ( : ), allocatable :: &
      RankIndex
    logical ( KDL ) :: &
      Initialized = .false.
    logical ( KDL ), private :: &
      MPI_Initialized = .false.
    character ( LDF ) :: &
      Name
    type ( CommunicatorForm ), pointer :: &
      Parent => null ( )
  contains
    procedure, private, pass :: &
      InitializeWorld
    procedure, private, pass :: &
      InitializeSubcommunicator
    generic :: &
      Initialize => InitializeWorld, InitializeSubcommunicator
    procedure, public, pass :: & 
      Synchronize
    procedure, public, pass :: &
      Abort => Abort_C  !-- avoids conflict with intrinsic "abort"
    final :: &
      Finalize
  end type CommunicatorForm
  
contains

  
  subroutine InitializeWorld ( C )

    class ( CommunicatorForm ), intent ( inout ) :: &
      C
    
    integer ( KDI ) :: &
      iRank, &
      Error
    logical ( KDL ) :: &
      Is_MPI_Initialized
      
    !-- Initialize MPI and the "World" communicator

    call MPI_INITIALIZED ( Is_MPI_Initialized, Error )
    if ( .not. Is_MPI_Initialized ) then
      call MPI_INIT ( C % Error )
      C % MPI_Initialized = .true.
    end if
    
    call MPI_COMM_DUP ( MPI_COMM_WORLD, C % Handle, C % Error )
    call MPI_COMM_SIZE ( C % Handle, C % Size, C % Error )
    call MPI_COMM_RANK ( C % Handle, C % Rank, C % Error )

    allocate ( C % RankIndex ( 0 : C % Size - 1 ) )
    C % RankIndex = [ ( iRank, iRank = 0, C % Size - 1 ) ]

    C % Name = 'World'

    call Show &
           ( 'Initializing a Communicator', CONSOLE % INFO_1, &
             DisplayRankOption = C % Rank )
    call Show &
           ( C % Name, 'Name', CONSOLE % INFO_1, &
             DisplayRankOption = C % Rank )
    call Show &
           ( C % Size, 'Size', CONSOLE % INFO_1, &
             DisplayRankOption = C % Rank )
             
    C % Initialized = .true. 

  end subroutine InitializeWorld
  
  
  subroutine InitializeSubcommunicator ( C, Parent, Ranks, NameOption )
  
    class ( CommunicatorForm ), intent ( inout ) :: &
      C
    type ( CommunicatorForm ), intent ( in ), target :: &
      Parent
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      Ranks
    character ( * ), intent ( in ), optional :: &
      NameOption
    
    integer ( KDI ) :: &
      iRank, &
      OldGroup, NewGroup
    
    !-- Create subcommunicator

    call MPI_COMM_GROUP ( Parent % Handle, OldGroup, C % Error )
    call MPI_GROUP_INCL &
           ( OldGroup, size ( Ranks ), Ranks, NewGroup, C % Error )
    call MPI_COMM_CREATE &
           ( Parent % Handle, NewGroup, C % Handle, C % Error )
    call MPI_GROUP_FREE ( NewGroup, C % Error )
    call MPI_GROUP_FREE ( OldGroup, C % Error )
    call MPI_COMM_SIZE ( C % Handle, C % Size, C % Error )
    call MPI_COMM_RANK ( C % Handle, C % Rank, C % Error )

    allocate ( C % RankIndex ( 0 : C % Size - 1 ) )
    C % RankIndex = [ ( iRank, iRank = 0, C % Size - 1 ) ]
    
    C % Name = ''
    if ( present ( NameOption ) ) C % Name = NameOption

    C % Parent => Parent

    call Show ( 'Initializing a Communicator', CONSOLE % INFO_2 )
    call Show ( C % Parent % Name, 'Parent', CONSOLE % INFO_2 )
    call Show ( C % Name, 'Name', CONSOLE % INFO_2 )
    call Show ( C % Size, 'Size', CONSOLE % INFO_2 )

    C % Initialized = .true. 
  
  end subroutine InitializeSubcommunicator
  
  
  subroutine Synchronize ( C )

    class ( CommunicatorForm ), intent ( inout ) :: &
      C

    call MPI_BARRIER ( C % Handle, C % Error )

  end subroutine Synchronize
  
  
  subroutine Abort_C ( C, AbortCodeOption )

    class ( CommunicatorForm ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ), optional :: &
      AbortCodeOption
    
    integer ( KDI ) :: &
      AbortCode
    
    AbortCode = -1
    if ( present ( AbortCodeOption ) ) AbortCode = AbortCodeOption

    call MPI_ABORT ( C % Handle, AbortCode, C % Error )

  end subroutine Abort_C
  
  
  impure elemental subroutine Finalize ( C )

    type ( CommunicatorForm ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      Error

    if ( C % Name == 'World' ) then
      call Show ( 'Finalizing a Communicator', CONSOLE % INFO_1 )
      call Show ( C % Name, 'Name', CONSOLE % INFO_1 )
      nullify ( C % Parent )
      if ( allocated ( C % RankIndex ) ) deallocate ( C % RankIndex )
      if ( C % MPI_Initialized ) call MPI_FINALIZE ( Error )
    else
      call Show ( 'Finalizing a Communicator', CONSOLE % INFO_2 )
      call Show ( C % Name, 'Name', CONSOLE % INFO_2 )
      nullify ( C % Parent )
      if ( allocated ( C % RankIndex ) ) deallocate ( C % RankIndex )
      call MPI_COMM_FREE ( C % Handle, Error )
    end if
    
    C % Initialized = .false.
        
  end subroutine Finalize


end module Communicator_Form
