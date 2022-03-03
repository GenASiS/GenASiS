!-- Real_3D_Form allows the construction of an array of 3D real
!   arrays to form ragged arrays.

module Real_3D__Form

  use iso_c_binding
  use Specifiers
  use Devices
  use ArrayOperations

  implicit none
  private

  type, public :: Real_3D_Form
    type ( c_ptr ), private :: &
      D_Value = c_null_ptr
    integer ( KDI ) :: &
      ErrorDevice
    real ( KDR ), dimension ( :, :, : ), allocatable :: &
      Value
    logical ( KDL ) :: &
      AllocatedDevice = .false. 
  contains
    procedure, private, pass :: &
      Initialize_R_3D
    procedure, private, pass :: &
      Initialize_R_3D_FromValue
    procedure, private, pass :: &
      Initialize_R_3D_Copy
    generic :: &
      Initialize &
        => Initialize_R_3D, Initialize_R_3D_FromValue, Initialize_R_3D_Copy
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_R_3D 
    procedure, public, pass :: &
      UpdateDevice => UpdateDevice_R_3D
    procedure, public, pass :: &
      UpdateHost => UpdateHost_R_3D
    procedure, public, pass :: &
      Clear => Clear_R_3D
    procedure, private, pass :: &
      MultiplyAdd_R_3D, &
      MultiplyAddInPlace_R_3D
    generic, public :: &
      MultiplyAdd => MultiplyAdd_R_3D, MultiplyAddInPlace_R_3D
    final :: &
      Finalize_R_3D
  end type Real_3D_Form
  
contains


  subroutine Initialize_R_3D ( A, nValues, ClearOption, iaLowerBoundOption )
    
    class ( Real_3D_Form ), intent ( inout ) :: &
      A
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nValues
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      iaLowerBoundOption

    integer ( KDI ), dimension ( 3 ) :: &
      iaLB
    logical ( KDL ) :: &
      ClearRequested

    if ( any ( nValues < 0 ) ) return
    
    if ( all ( nValues == 0 ) ) then
      allocate ( A % Value ( 0, 0, 0 ) )
      return
    end if 
    
    ClearRequested = .false.
    if ( present ( ClearOption ) ) ClearRequested = ClearOption

    iaLB = 1
    if ( present ( iaLowerBoundOption ) ) iaLB = iaLowerBoundOption
    
    allocate &
      ( A % Value &
          ( iaLB ( 1 ) : iaLB ( 1 ) + nValues ( 1 ) - 1, &
            iaLB ( 2 ) : iaLB ( 2 ) + nValues ( 2 ) - 1, &
            iaLB ( 3 ) : iaLB ( 3 ) + nValues ( 3 ) - 1 ) )
    
    if ( ClearRequested ) call Clear ( A % Value )

  end subroutine Initialize_R_3D
  
  
  subroutine Initialize_R_3D_FromValue ( A, Value, iaLowerBoundOption )
    
    class ( Real_3D_Form ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      Value
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      iaLowerBoundOption

    call A % Initialize_R_3D &
           ( shape ( Value ), iaLowerBoundOption = iaLowerBoundOption )
    A % Value = Value 

  end subroutine Initialize_R_3D_FromValue
  
  
  subroutine Initialize_R_3D_Copy ( A, B, iaLowerBoundOption )
    
    class ( Real_3D_Form ), intent ( inout ) :: &
      A
    type (  Real_3D_Form ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ), optional :: &
      iaLowerBoundOption
      
    integer ( KDI ), dimension ( 3 ) :: &
      iaLB
    
    iaLB = lbound ( B % Value ) 
    if ( present ( iaLowerBoundOption ) ) iaLB = iaLowerBoundOption

    call A % Initialize_R_3D_FromValue ( B % Value, iaLowerBoundOption = iaLB )
    
    if ( B % AllocatedDevice ) then
      call A % AllocateDevice ( )
      call Copy ( B % Value, A % Value, UseDeviceOption = B % AllocatedDevice )
    end if
  
  end subroutine Initialize_R_3D_Copy 
  
  
  impure elemental subroutine AllocateDevice_R_3D ( A )
  
    class ( Real_3D_Form ), intent ( inout ) :: &
      A
      
    call AllocateDevice ( size ( A % Value ), A % D_Value )
    A % AllocatedDevice = .true.
    call AssociateHost ( A % D_Value, A % Value )
  
  end subroutine AllocateDevice_R_3D
  
  
  impure elemental subroutine UpdateDevice_R_3D ( A )
  
    class ( Real_3D_Form ), intent ( inout ) :: &
      A
    
    if ( .not. A % AllocatedDevice ) &
      return
    
    call UpdateDevice &
           ( A % Value, A % D_Value, ErrorOption = A % ErrorDevice )
  
  end subroutine UpdateDevice_R_3D


  impure elemental subroutine UpdateHost_R_3D ( A )
  
    class ( Real_3D_Form ), intent ( inout ) :: &
      A
    
    if ( .not. A % AllocatedDevice ) &
      return
    
    call UpdateHost &
           ( A % D_Value, A % Value, ErrorOption = A % ErrorDevice )
  
  end subroutine UpdateHost_R_3D


  impure elemental subroutine Clear_R_3D ( A )
  
    class ( Real_3D_Form ), intent ( inout ) :: &
      A
    
    call Clear ( A % Value, UseDeviceOption = A % AllocatedDevice )
  
  end subroutine Clear_R_3D


  impure elemental subroutine MultiplyAdd_R_3D &
                     ( R_3D_D, R_3D_A, R_3D_B, C, UseDeviceOption )

    class ( Real_3D_Form ), intent ( inout ) :: &
      R_3D_D
    class ( Real_3D_Form ), intent ( in ) :: &
      R_3D_A, &
      R_3D_B
    real ( KDR ), intent ( in ) :: &
      C
    logical ( KDL ), intent ( in ), optional :: &
       UseDeviceOption
    
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  R_3D_D % AllocatedDevice
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    associate &
      ( A  =>  R_3D_A % Value, &
        B  =>  R_3D_B % Value, &
        D  =>  R_3D_D % Value )

    call MultiplyAdd ( A, B, C, D, UseDeviceOption = UseDevice )

    end associate !-- A, etc.

  end subroutine MultiplyAdd_R_3D


  impure elemental subroutine MultiplyAddInPlace_R_3D &
                     ( R_3D_A, R_3D_B, C, UseDeviceOption )

    class ( Real_3D_Form ), intent ( inout ) :: &
      R_3D_A
    class ( Real_3D_Form ), intent ( in ) :: &
      R_3D_B
    real ( KDR ), intent ( in ) :: &
      C
    logical ( KDL ), intent ( in ), optional :: &
       UseDeviceOption
    
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  R_3D_A % AllocatedDevice
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    associate &
      ( A  =>  R_3D_A % Value, &
        B  =>  R_3D_B % Value )

    call MultiplyAdd ( A, B, C, UseDeviceOption = UseDevice )

    end associate !-- A, etc.
        
  end subroutine MultiplyAddInPlace_R_3D


  impure elemental subroutine Finalize_R_3D ( A )

    type ( Real_3D_Form ), intent ( inout ) :: &
      A

    if ( allocated ( A % Value ) ) &
      deallocate ( A % Value )
    
    if ( A % AllocatedDevice ) &
      call DeallocateDevice ( A % D_Value )

  end subroutine Finalize_R_3D
  
  
end module Real_3D__Form
