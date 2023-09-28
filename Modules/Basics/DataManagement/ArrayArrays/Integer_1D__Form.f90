!-- Integer_1D_Form allows the construction of an array of 1D integer
!   arrays to form ragged arrays.

module Integer_1D__Form

  use iso_c_binding
  use Specifiers
  use ArrayOperations

  implicit none
  private

  type, public :: Integer_1D_Form
    type ( c_ptr ), private :: &
      D_Value = c_null_ptr
    integer (KDI ) :: &
      ErrorDevice
    integer ( KDI ), dimension ( : ), allocatable :: &
      Value
    logical ( KDL ) :: &
      AllocatedDevice = .false.
  contains
    procedure, private, pass :: &
      Initialize_I_1D
    procedure, private, pass :: &
      Initialize_I_1D_FromValue
    procedure, private, pass :: &
      Initialize_I_1D_Copy
    generic :: &
      Initialize &
        => Initialize_I_1D, Initialize_I_1D_FromValue, Initialize_I_1D_Copy
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_I_1D
    procedure, public, pass :: &
      UpdateDevice => UpdateDevice_I_1D
    procedure, public, pass :: &
      UpdateHost => UpdateHost_I_1D
    final :: &
      Finalize_I_1D
  end type Integer_1D_Form
  
contains


  subroutine Initialize_I_1D ( A, nValues, ClearOption, iLowerBoundOption )
    
    class ( Integer_1D_Form ), intent ( inout ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      nValues
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    integer ( KDI ), intent ( in ), optional :: &
      iLowerBoundOption

    integer ( KDI ) :: &
      iLB
    logical ( KDL ) :: &
      ClearRequested

    if ( nValues < 0 ) return
    
    if ( nValues == 0 ) then
      allocate ( A % Value ( 0 ) )
      return
    end if 
    
    ClearRequested = .false.
    if ( present ( ClearOption ) ) ClearRequested = ClearOption

    iLB = 1
    if ( present ( iLowerBoundOption ) ) iLB = iLowerBoundOption
    
    allocate ( A % Value ( iLB : iLB + nValues - 1 ) )
    
    if ( ClearRequested ) call Clear ( A % Value )

  end subroutine Initialize_I_1D
  
  
  subroutine Initialize_I_1D_FromValue ( A, Value, iLowerBoundOption )
    
    class ( Integer_1D_Form ), intent ( inout ) :: &
      A
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      Value
    integer ( KDI ), intent ( in ), optional :: &
      iLowerBoundOption

    call A % Initialize_I_1D &
           ( size ( Value ), iLowerBoundOption = iLowerBoundOption )
    A % Value = Value 

  end subroutine Initialize_I_1D_FromValue
  
  
  subroutine Initialize_I_1D_Copy ( A, B, iLowerBoundOption )
    
    class ( Integer_1D_Form ), intent ( inout ) :: &
      A
    type (  Integer_1D_Form ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ), optional :: &
      iLowerBoundOption
      
    integer ( KDI ) :: &
      iLB
    
    iLB = lbound ( B % Value, dim = 1 ) 
    if ( present ( iLowerBoundOption ) ) iLB = iLowerBoundOption

    call A % Initialize_I_1D_FromValue ( B % Value, iLowerBoundOption = iLB )

    if ( B % AllocatedDevice ) then
      call A % AllocateDevice ( )
      call Copy ( B % Value, A % Value, UseDeviceOption = B % AllocatedDevice )
    end if
  
  end subroutine Initialize_I_1D_Copy 


  impure elemental subroutine AllocateDevice_I_1D ( A )
  
    class ( Integer_1D_Form ), intent ( inout ) :: &
      A
      
    if ( .not. allocated ( A % Value ) ) &
      return
       
    call AllocateDevice ( size ( A % Value ), A % D_Value )
    A % AllocatedDevice = .true.
    call AssociateHost ( A % D_Value, A % Value )
  
  end subroutine AllocateDevice_I_1D
  
  
  impure elemental subroutine UpdateDevice_I_1D ( A )
  
    class ( Integer_1D_Form ), intent ( inout ) :: &
      A
       
    if ( .not. A % AllocatedDevice ) &
      return
    
    call UpdateDevice &
           ( A % Value, A % D_Value, ErrorOption = A % ErrorDevice )
  
  end subroutine UpdateDevice_I_1D
  
  
  impure elemental subroutine UpdateHost_I_1D ( A )
  
    class ( Integer_1D_Form ), intent ( inout ) :: &
      A
       
    if ( .not. A % AllocatedDevice ) &
      return
    
    call UpdateHost &
           ( A % D_Value, A % Value, ErrorOption = A % ErrorDevice )
  
  end subroutine UpdateHost_I_1D
  
  
  impure elemental subroutine Finalize_I_1D ( A )

    type ( Integer_1D_Form ), intent ( inout ) :: &
      A

    if ( A % AllocatedDevice ) then
      call DisassociateHost ( A % Value ) 
      call DeallocateDevice ( A % D_Value )
    end if

    if ( allocated ( A % Value ) ) &
      deallocate ( A % Value )

  end subroutine Finalize_I_1D
  
  
end module Integer_1D__Form
