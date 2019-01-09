!-- Message_R_Form provides the concrete type for Message object for
!   real datatype

module Message_R__Form

  use MPI
  use Specifiers
  use Devices
  use MessagePassingBasics
  use Message_Template

  implicit none
  private

  type, public, extends ( MessageTemplate ) :: Message_R_Form 
    real ( KDR ), dimension ( : ), pointer :: &
      Value => null ( )
  contains
    procedure, public, pass :: &
      InitializeAllocate
    procedure, public, pass :: &
      InitializeAssociate
    generic :: &
      Initialize => InitializeAllocate, InitializeAssociate
    procedure, private, pass :: &
      AllocateDevice_M
    final :: &
      Finalize
  end type Message_R_Form

contains


  subroutine InitializeAllocate ( M, C, Tag, Rank, nValues )

    class ( Message_R_Form ), intent ( inout ), target :: &
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

    class ( Message_R_Form ), intent ( inout ), target :: &
      M
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( in ) :: &
      Tag, &
      Rank

    if ( M % Initialized ) return
    
    call M % InitializeTemplate ( C, Tag, Rank, AllocatedOption = .false. )

    M % Value => Value
  
  end subroutine InitializeAssociate
  
  
  subroutine AllocateDevice_M ( M )
    
    class ( Message_R_Form ), intent ( inout ), target :: &
      M
    
    call AllocateDevice ( M % Value, M % D_Value )
    
  end subroutine AllocateDevice_M
  
  
  elemental subroutine Finalize ( M )

    type ( Message_R_Form ), intent ( inout ) :: &
      M 
      
    if ( M % AllocatedValue ) then
      if ( associated ( M % Value ) ) deallocate ( M % Value )
    end if
    nullify ( M % Value )

    nullify ( M % Communicator )
    
  end subroutine Finalize
  
  
end module Message_R__Form
