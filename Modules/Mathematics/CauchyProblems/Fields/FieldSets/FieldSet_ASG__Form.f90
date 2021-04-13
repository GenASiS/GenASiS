module FieldSet_ASG__Form

  !-- FieldSet_AtlasSingleGrid__Form

  use Basics
  use Manifolds
  use FieldSet_GS__Form
  use FieldSet_AH__Form

  implicit none
  private

  type, public, extends ( FieldSet_AH_Form ) :: FieldSet_ASG_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_ASG
    generic, public :: &
      Initialize => InitializeAllocate_ASG
  end type FieldSet_ASG_Form

contains


  subroutine InitializeAllocate_ASG ( FSA )!, ASG )

    class ( FieldSet_ASG_Form ), intent ( inout ) :: &
      FSA
!    class ( Atlas_SG_Form ), intent ( in ) :: &
!      ASG

  end subroutine InitializeAllocate_ASG


end module FieldSet_ASG__Form
