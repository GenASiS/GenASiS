module HomogeneousSpheroid_Form

  use Basics
  use Mathematics
  
  implicit none

  type, public :: HomogeneousSpheroidForm
    type ( GridImageStreamForm ), allocatable :: &
      Stream
    type ( Atlas_SC_CC_Form ), allocatable :: &
      Atlas
    type ( Storage_ASC_Form ), allocatable :: &
      Source, &
      Solution, &
      Reference, &
      Difference
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson
  contains
    procedure, public, pass :: &
      Initialize
   ! procedure, public, pass :: &
   !   ComputeError
   ! final :: &
   !   Finalize
  end type HomogeneousSpheroidForm

contains


  subroutine Initialize ( HS, Name )
    
    class ( HomogeneousSpheroidForm ) :: &
      HS
    character ( * ), intent ( in ) :: &
      Name
  
  end subroutine Initialize 


end module HomogeneousSpheroid_Form
