module FieldSet_ASG__Form

  !-- FieldSet_AtlasSingleGrid__Form

  use Basics
  use Manifolds
  use FieldSet_GS__Form
  use FieldSet_AH__Form

  implicit none
  private

  type, public, extends ( FieldSet_AH_Form ) :: FieldSet_ASG_Form
    class ( FieldSet_GS_Form ), pointer :: &
      FieldSet_G => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_ASG
    generic, public :: &
      Initialize => InitializeAllocate_ASG
    final :: &
      Finalize
  end type FieldSet_ASG_Form

contains


  subroutine InitializeAllocate_ASG ( FSA, A, NameOption, IgnorabilityOption )

    class ( FieldSet_ASG_Form ), intent ( inout ) :: &
      FSA
    class ( Atlas_SG_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FSA % Type  ==  '' ) &
      FSA % Type  =  'a FieldSet_ASG'

    call FSA % Initialize_H &
           ( A, &
             NameOption = NameOption, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_ASG


  impure elemental subroutine Finalize ( FSA )

    type ( FieldSet_ASG_Form ), intent ( inout ) :: &
      FSA

    nullify ( FSA % FieldSet_G )

  end subroutine Finalize

end module FieldSet_ASG__Form
