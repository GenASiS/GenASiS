!-- Integer_2D_Form allows the construction of an array of 2D integer
!   arrays to form ragged arrays.

module Integer_2D__Form
  
  use iso_c_binding
  use Specifiers
  use Devices
  use ArrayOperations

  implicit none
  private

  type, public :: Integer_2D_Form
    type ( c_ptr ), private :: &
      D_Value = c_null_ptr
    integer ( KDI ) :: &
      ErrorDevice
    integer ( KDI ), dimension ( :, : ), allocatable :: &
      Value
    logical ( KDL ) :: &
      AllocatedDevice = .false.
  contains
    procedure, private, pass :: &
      Initialize_I_2D
    procedure, private, pass :: &
      Initialize_I_2D_FromValue
    procedure, private, pass :: &
      Initialize_I_2D_Copy
    generic :: &
      Initialize &
        => Initialize_I_2D, Initialize_I_2D_FromValue, Initialize_I_2D_Copy
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_I_2D
    procedure, public, pass :: &
      UpdateDevice => UpdateDevice_I_2D
    procedure, public, pass :: &
      UpdateHost => UpdateHost_I_2D
    final :: &
      Finalize_I_2D
  end type Integer_2D_Form
  
contains


  subroutine Initialize_I_2D ( A, nValues, ClearOption, iaLowerBoundOption )
    
    class ( Integer_2D_Form ), intent ( inout ) :: &
      A
    integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
      nValues
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    integer ( KDI ), dimension ( 2 ), intent ( in ), optional :: &
      iaLowerBoundOption

    integer ( KDI ), dimension ( 2 ) :: &
      iaLB
    logical ( KDL ) :: &
      ClearRequested

    if ( any ( nValues < 0 ) ) return
    
    if ( all ( nValues == 0 ) ) then
      allocate ( A % Value ( 0, 0 ) )
      return
    end if 
    
    ClearRequested = .false.
    if ( present ( ClearOption ) ) ClearRequested = ClearOption

    iaLB = 1
    if ( present ( iaLowerBoundOption ) ) iaLB = iaLowerBoundOption
    
    allocate &
      ( A % Value &
          ( iaLB ( 1 ) : iaLB ( 1 ) + nValues ( 1 ) - 1, &
            iaLB ( 2 ) : iaLB ( 2 ) + nValues ( 2 ) - 1 ) )
    
    if ( ClearRequested ) call Clear ( A % Value )

  end subroutine Initialize_I_2D
  
  
  subroutine Initialize_I_2D_FromValue ( A, Value, iaLowerBoundOption )
    
    class ( Integer_2D_Form ), intent ( inout ) :: &
      A
    integer ( KDI ), dimension ( :, : ), intent ( in ) :: &
      Value
    integer ( KDI ), dimension ( 2 ), intent ( in ), optional :: &
      iaLowerBoundOption

    call A % Initialize_I_2D &
           ( shape ( Value ), iaLowerBoundOption = iaLowerBoundOption )
    A % Value = Value 

  end subroutine Initialize_I_2D_FromValue
  
  
  subroutine Initialize_I_2D_Copy ( A, B, iaLowerBoundOption )
    
    class ( Integer_2D_Form ), intent ( inout ) :: &
      A
    type ( Integer_2D_Form ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ), optional :: &
      iaLowerBoundOption
      
    integer ( KDI ), dimension ( 2 ) :: &
      iaLB
    
    iaLB = lbound ( B % Value ) 
    if ( present ( iaLowerBoundOption ) ) iaLB = iaLowerBoundOption

    call A % Initialize_I_2D_FromValue ( B % Value, iaLowerBoundOption = iaLB )
    
    if ( B % AllocatedDevice ) then
      call A % AllocateDevice ( )
      call Copy ( B % Value, A % Value, UseDeviceOption = B % AllocatedDevice )
    end if
  
  end subroutine Initialize_I_2D_Copy 
  
  
  impure elemental subroutine AllocateDevice_I_2D ( A )
  
    class ( Integer_2D_Form ), intent ( inout ) :: &
      A
    
    if ( .not. allocated ( A % Value ) ) &
      return
      
    call AllocateDevice ( size ( A % Value ), A % D_Value )
    A % AllocatedDevice = .true.
    call AssociateHost ( A % D_Value, A % Value )
  
  end subroutine AllocateDevice_I_2D 
  
  
  impure elemental subroutine UpdateDevice_I_2D ( A )
  
    class ( Integer_2D_Form ), intent ( inout ) :: &
      A
       
    if ( .not. A % AllocatedDevice ) &
      return
    
    call UpdateDevice &
           ( A % Value, A % D_Value, ErrorOption = A % ErrorDevice )
  
  end subroutine UpdateDevice_I_2D
  
  
  impure elemental subroutine UpdateHost_I_2D ( A )
  
    class ( Integer_2D_Form ), intent ( inout ) :: &
      A
       
    if ( .not. A % AllocatedDevice ) &
      return
    
    call UpdateHost &
           ( A % D_Value, A % Value, ErrorOption = A % ErrorDevice )
  
  end subroutine UpdateHost_I_2D
  
  
  impure elemental subroutine Finalize_I_2D ( A )

    type ( Integer_2D_Form ), intent ( inout ) :: &
      A

    if ( A % AllocatedDevice ) then
      call DisassociateHost ( A % Value ) 
      call DeallocateDevice ( A % D_Value )
    end if

    if ( allocated ( A % Value ) ) &
      deallocate ( A % Value )

  end subroutine Finalize_I_2D
  
  
end module Integer_2D__Form
