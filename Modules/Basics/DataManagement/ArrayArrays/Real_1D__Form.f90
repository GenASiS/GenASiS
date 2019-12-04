!-- Real_1D_Form allows the construction of an array of 1D real
!   arrays to form ragged arrays.

module Real_1D__Form

  use iso_c_binding
  use Specifiers
  use Devices
  use ArrayOperations

  implicit none
  private

  type, public :: Real_1D_Form
    type ( c_ptr ), private :: &
      D_Value = c_null_ptr
    integer ( KDI ) :: &
      ErrorDevice
    real ( KDR ), dimension ( : ), allocatable :: &
      Value
    logical ( KDL ) :: &
      AllocatedDevice = .false. 
  contains
    procedure, private, pass :: &
      Initialize_R_1D
    procedure, private, pass :: &
      Initialize_R_1D_FromValue
    procedure, private, pass :: &
      Initialize_R_1D_Copy
    generic :: &
      Initialize &
        => Initialize_R_1D, Initialize_R_1D_FromValue, Initialize_R_1D_Copy
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_R_1D
    procedure, public, pass :: &
      UpdateDevice => UpdateDevice_R_1D
    procedure, public, pass :: &
      UpdateHost => UpdateHost_R_1D
    final :: &
      Finalize_R_1D
  end type Real_1D_Form
  
contains


  subroutine Initialize_R_1D ( A, nValues, ClearOption, iLowerBoundOption )
    
    class ( Real_1D_Form ), intent ( inout ) :: &
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

  end subroutine Initialize_R_1D
  
  
  subroutine Initialize_R_1D_FromValue ( A, Value, iLowerBoundOption )
    
    class ( Real_1D_Form ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Value
    integer ( KDI ), intent ( in ), optional :: &
      iLowerBoundOption

    call A % Initialize_R_1D &
           ( size ( Value ), iLowerBoundOption = iLowerBoundOption )
    A % Value = Value 

  end subroutine Initialize_R_1D_FromValue
  
  
  subroutine Initialize_R_1D_Copy ( A, B, iLowerBoundOption )
    
    class ( Real_1D_Form ), intent ( inout ) :: &
      A
    type (  Real_1D_Form ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ), optional :: &
      iLowerBoundOption
      
    integer ( KDI ) :: &
      iLB
    
    iLB = lbound ( B % Value, dim = 1 ) 
    if ( present ( iLowerBoundOption ) ) iLB = iLowerBoundOption

    call A % Initialize_R_1D_FromValue ( B % Value, iLowerBoundOption = iLB )
    
    if ( B % AllocatedDevice ) then
      call A % AllocateDevice ( )
      call Copy ( B % Value, A % Value, UseDeviceOption = B % AllocatedDevice )
    end if
  
  end subroutine Initialize_R_1D_Copy
  
  
  impure elemental subroutine AllocateDevice_R_1D ( A )
  
    class ( Real_1D_Form ), intent ( inout ) :: &
      A
      
    call AllocateDevice ( size ( A % Value ), A % D_Value )
    A % AllocatedDevice = .true.
    call AssociateHost ( A % D_Value, A % Value )
  
  end subroutine AllocateDevice_R_1D
  
  
  impure elemental subroutine UpdateDevice_R_1D ( A )
  
    class ( Real_1D_Form ), intent ( inout ) :: &
      A
    
    if ( .not. A % AllocatedDevice ) &
      return
    
    call UpdateDevice &
           ( A % Value, A % D_Value, ErrorOption = A % ErrorDevice )
  
  end subroutine UpdateDevice_R_1D
  

  impure elemental subroutine UpdateHost_R_1D ( A )
  
    class ( Real_1D_Form ), intent ( inout ) :: &
      A
    
    if ( .not. A % AllocatedDevice ) &
      return
    
    call UpdateHost &
           ( A % D_Value, A % Value, ErrorOption = A % ErrorDevice )
  
  end subroutine UpdateHost_R_1D

  
  impure elemental subroutine Finalize_R_1D ( A )

    type ( Real_1D_Form ), intent ( inout ) :: &
      A

    if ( allocated ( A % Value ) ) &
      deallocate ( A % Value )
    
    if ( A % AllocatedDevice ) &
      call DeallocateDevice ( A % D_Value )

  end subroutine Finalize_R_1D
  
  
end module Real_1D__Form
