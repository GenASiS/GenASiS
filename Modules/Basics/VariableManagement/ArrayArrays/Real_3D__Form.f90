!-- Real_3D_Form allows the construction of an array of 3D real
!   arrays to form ragged arrays.

module Real_3D__Form

  use Specifiers
  use ArrayOperations

  implicit none
  private

  type, public :: Real_3D_Form
    real ( KDR ), dimension ( :, :, : ), allocatable :: &
      Value
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
  
  end subroutine Initialize_R_3D_Copy 
  
  
  elemental subroutine Finalize_R_3D ( A )

    type ( Real_3D_Form ), intent ( inout ) :: &
      A

    if ( allocated ( A % Value ) ) deallocate ( A % Value )

  end subroutine Finalize_R_3D
  
  
end module Real_3D__Form
