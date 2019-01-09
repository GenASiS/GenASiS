!-- Message_1D_Template provides an abstraction commonly shared by 
!   higher-level, concrete type object with specified datatype to array of 
!   messages.

module Message_1D__Template

  use MPI
  use Specifiers
  use Devices
  use Message_Template
 
  implicit none
  private

  type, public, abstract :: Message_1D_Template
    integer ( KDI ) :: &
      nMessages, &
      Error
    class ( MessageTemplate ), dimension ( : ), pointer :: &
      MessageTemplate => null ( )
  contains
    procedure, public, pass :: &
      AllocateDevice_M_1D
    generic :: &
      AllocateDevice => AllocateDevice_M_1D
    procedure, public, pass :: &
      WaitOne
    procedure, public, pass :: &
      WaitAny
    procedure, public, pass :: &
      WaitAll
    generic :: &
      Wait => WaitOne, WaitAny, WaitAll
  end type Message_1D_Template
  
contains


  subroutine AllocateDevice_M_1D ( M_1D )
    
    class ( Message_1D_Template ), intent ( inout ), target :: &
      M_1D
      
    integer ( KDI ) :: &
      iM    !-- iMessage

    do iM = 1, M_1D % nMessages
      call M_1D % MessageTemplate ( iM ) % AllocateDevice ( )
    end do

  end subroutine AllocateDevice_M_1D
  

  subroutine WaitOne ( M_1D, Index )

    class ( Message_1D_Template ), intent ( inout ) :: &
      M_1D
    integer ( KDI ), intent ( in ) :: &
      Index

    call M_1D % MessageTemplate ( Index ) % Wait ( )

  end subroutine WaitOne


  subroutine WaitAny ( M_1D, AllFinished, Index )
  
    class ( Message_1D_Template ), intent ( inout ) :: &
      M_1D
    logical ( KDL ), intent ( out ) :: &
      AllFinished
    integer ( KDI ), intent ( out ) :: &
      Index
    
    integer ( KDI ) :: &
      iM    !-- iMessage
    integer ( KDI ), dimension ( M_1D % nMessages ) :: &
      Handle
      
    if ( M_1D % nMessages <= 0 ) then
      AllFinished = .true.
      return
    end if
    
    do iM = 1, M_1D % nMessages
      Handle ( iM ) = M_1D % MessageTemplate ( iM ) % Handle
    end do
    
    call MPI_WAITANY &
           ( M_1D % nMessages, Handle, Index, MPI_STATUS_IGNORE, &
             M_1D % Error )
    
    AllFinished = .false.
    if ( Index == MPI_UNDEFINED ) then
      AllFinished = .true.
      return
    else 
      M_1D % MessageTemplate ( Index ) % Handle = Handle ( Index )
    end if
  
  end subroutine WaitAny
  
  
  subroutine WaitAll ( M_1D )
    
    class ( Message_1D_Template ), intent ( inout ) :: &
      M_1D
    
    integer ( KDI ) :: &
      iM    !-- iMessage
    integer ( KDI ), dimension ( M_1D % nMessages ) :: &
      Handle
    
    do iM = 1, M_1D % nMessages
      Handle ( iM ) = M_1D % MessageTemplate ( iM ) % Handle
    end do
          
    if ( M_1D % nMessages <= 0 ) return
    
    call MPI_WAITALL &
           ( M_1D % nMessages, Handle, MPI_STATUSES_IGNORE, M_1D % Error )
    
    if ( M_1D % Error == MPI_SUCCESS ) then
      do iM = 1, M_1D % nMessages
        M_1D % MessageTemplate ( iM ) % Handle = Handle ( iM )
      end do
    end if
    
  end subroutine WaitAll
 

end module Message_1D__Template
