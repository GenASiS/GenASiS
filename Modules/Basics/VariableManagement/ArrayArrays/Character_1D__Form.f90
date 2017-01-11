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
    procedure, public, pass :: &
      Initialize => Initialize_C_1D
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


  elemental subroutine Finalize_C_1D ( A )

    type ( Character_1D_Form ), intent ( inout ) :: &
      A

    if ( allocated ( A % Value ) ) deallocate ( A % Value )

  end subroutine Finalize_C_1D


end module Character_1D__Form
