module Relaxation_RM_S__Form

  use Basics
  use Mathematics
  use Relaxation_RM__Template

  implicit none
  private
  
  type, public, extends ( Relaxation_RM_Template ) :: Relaxation_RM_S_Form
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Relaxation_RM_S_Form

contains

  subroutine Initialize ( R, Name )

    class ( Relaxation_RM_S_Form ), intent ( inout ) :: &
      R
    character ( * ), intent ( in ) :: &
      Name

    if ( R % Type == '' ) &
      R % Type  =  'a Relaxation_RM_S'

    call R % InitializeTemplate ( Name )

  end subroutine Initialize


  subroutine Finalize ( R )

    type ( Relaxation_RM_S_Form ), intent ( inout ) :: &
      R

    call R % FinalizeTemplate ( )

  end subroutine Finalize


end module Relaxation_RM_S__Form
