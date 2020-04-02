!-- Real_3D_2D_Form allows the construction of an array of 2D arrays of 3D real
!   arrays to form ragged arrays.

module Real_3D_2D__Form

  use Specifiers
  use Real_3D__Form

  implicit none
  private

  type, public :: Real_3D_2D_Form
    type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
      Array
  contains
    procedure, private, pass :: &
      Initialize_R_3D_2D
    generic :: &
      Initialize => Initialize_R_3D_2D
    final :: &
      Finalize_R_3D_2D
  end type Real_3D_2D_Form
  
contains


  subroutine Initialize_R_3D_2D ( A, nArrays )
    
    class ( Real_3D_2D_Form ), intent ( inout ) :: &
      A
    integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
      nArrays

    allocate ( A % Array ( nArrays ( 1 ), nArrays ( 2 ) ) )
    
  end subroutine Initialize_R_3D_2D
  
  
  impure elemental subroutine Finalize_R_3D_2D ( A )

    type ( Real_3D_2D_Form ), intent ( inout ) :: &
      A

    if ( allocated ( A % Array ) ) &
      deallocate ( A % Array )

  end subroutine Finalize_R_3D_2D
  
  
end module Real_3D_2D__Form
