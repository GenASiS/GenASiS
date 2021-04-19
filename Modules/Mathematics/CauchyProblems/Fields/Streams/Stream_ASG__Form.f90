module Stream_ASG__Form

  !-- Stream_AtlasSingleGrid__Form

  use Basics
  use Manifolds
  use Stream_GS__Form
  use Stream_AH__Form

  implicit none
  private

  type, public, extends ( Stream_AH_Form ) :: Stream_ASG_Form
    class ( Stream_GS_Form ), pointer :: &
      Stream_G => null ( )
  contains
   procedure, public, pass :: &
     Initialize
    final :: &
      Finalize
  end type Stream_ASG_Form

contains


  subroutine Initialize ( SA, A, GIS, NameOption, VerboseOption )

    class ( Stream_ASG_Form ), intent ( inout ), target :: &
      SA
    class ( Atlas_SG_Form ), intent ( inout ) :: &
      A
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

!     !-- FIXME: This shouldn't be necessary, but for some reason GCC 10.1.0
!     !          doesn't compile without it
!     class ( Grid_S_Form ), pointer :: &
!       G_Pointer

    if ( SA % Type  ==  '' ) &
      SA % Type  =  'a Stream_ASG'

    call SA % Initialize_H ( A, NameOption )

    allocate ( Stream_GS_Form :: SA % Stream_C ( 1 ) % Element )
    select type ( SG  =>  SA % Stream_C ( 1 ) % Element )
    class is ( Stream_GS_Form )

    select type ( G  =>  A % Chart ( 1 ) % Element )
    class is ( Grid_S_Form )

!     !-- FIXME: See FIXME above
!     G_Pointer  =>  G
    call SG % Initialize ( G, GIS, NameOption, VerboseOption )

    SA % Stream_G  =>  SG

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Stream_ASG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
    end select !-- FSG

  end subroutine Initialize


  impure elemental subroutine Finalize ( SA )

    type ( Stream_ASG_Form ), intent ( inout ) :: &
      SA

    nullify ( SA % Stream_G )

  end subroutine Finalize

end module Stream_ASG__Form
