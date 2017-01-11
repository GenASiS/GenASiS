!-- Message_I_Form provides the concrete type for Message object for
!   integer datatype

module Message_I__Form

  use MPI
  use VariableManagement
  use MessagePassingBasics
  use Message_Template

  implicit none
  private

  type, public, extends ( MessageTemplate ) :: Message_I_Form 
    integer ( KDI ), dimension ( : ), pointer :: &
      Value => null ( )
  contains
    procedure, public, pass :: &
      InitializeAllocate
    procedure, public, pass :: &
      InitializeAssociate
    generic :: &
      Initialize => InitializeAllocate, InitializeAssociate
    final :: &
      Finalize
  end type Message_I_Form
  
contains


  subroutine InitializeAllocate ( M, C, Tag, Rank, nValues )

    class ( Message_I_Form ), intent ( inout ), target :: &
      M
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      Tag, &
      Rank, &
      nValues

    if ( M % Initialized ) return
    
    call M % InitializeTemplate ( C, Tag, Rank )
    
    allocate ( M % Value ( nValues ) )
  
  end subroutine InitializeAllocate
  
  
  subroutine InitializeAssociate ( M, C, Value, Tag, Rank )

    class ( Message_I_Form ), intent ( inout ), target :: &
      M
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    integer ( KDI ), dimension ( : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( in ) :: &
      Tag, &
      Rank

    if ( M % Initialized ) return
    
    call M % InitializeTemplate ( C, Tag, Rank, AllocatedOption = .false. )
    
    M % Value => Value
  
  end subroutine InitializeAssociate
  
  
  elemental subroutine Finalize ( M )

    type ( Message_I_Form ), intent ( inout ) :: &
      M 
      
    if ( M % AllocatedValue ) then
      if ( associated ( M % Value ) ) deallocate ( M % Value )
    end if
    nullify ( M % Value )

    nullify ( M % Communicator )

  end subroutine Finalize
  

end module Message_I__Form
