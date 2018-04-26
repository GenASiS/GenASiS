!-- Character_1D_Form allows the construction of an array of 1D character
!   arrays to form ragged arrays.

module Character_1D__Form

  use Specifiers

  implicit none
  private

  type, public :: Character_1D_Form
    character ( LDL ), dimension ( : ), allocatable :: &
      Value
  contains
    procedure, private, pass :: &
      Initialize_C_1D
    procedure, private, pass :: &
      Initialize_C_1D_FromValue
    generic :: &
      Initialize => Initialize_C_1D, Initialize_C_1D_FromValue
    final :: &
      Finalize_C_1D
  end type Character_1D_Form
	
contains


  subroutine Initialize_C_1D ( A, nValues, iLowerBoundOption )

    class ( Character_1D_Form ), intent ( inout ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      iLowerBoundOption

    integer ( KDI ) :: &
      LowerBound

    if ( nValues < 0 ) return
 
    LowerBound = 1
    if ( present ( iLowerBoundOption ) ) LowerBound = iLowerBoundOption
 
    if ( nValues > 0 )then
      allocate ( A % Value ( LowerBound : LowerBound + nValues - 1 ) )
    else
      allocate ( A % Value ( 0 ) )
    end if

  end subroutine Initialize_C_1D


  subroutine Initialize_C_1D_FromValue ( A, Value, iLowerBoundOption )
    
    class ( Character_1D_Form ), intent ( inout ) :: &
      A
    character ( * ), dimension ( : ), intent ( in ) :: &
      Value
    integer ( KDI ), intent ( in ), optional :: &
      iLowerBoundOption

    call A % Initialize_C_1D &
           ( size ( Value ), iLowerBoundOption = iLowerBoundOption )
    A % Value = Value 

  end subroutine Initialize_C_1D_FromValue
  
  
  elemental subroutine Finalize_C_1D ( A )

    type ( Character_1D_Form ), intent ( inout ) :: &
      A

    if ( allocated ( A % Value ) ) deallocate ( A % Value )

  end subroutine Finalize_C_1D


end module Character_1D__Form
