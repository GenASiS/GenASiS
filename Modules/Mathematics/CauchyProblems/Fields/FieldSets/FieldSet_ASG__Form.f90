module FieldSet_ASG__Form

  !-- FieldSet_AtlasSingleGrid__Form

  use Basics
  use Manifolds
  use FieldSet_CGS__Form
  use FieldSet_AH__Form

  implicit none
  private

  type, public, extends ( FieldSet_AH_Form ) :: FieldSet_ASG_Form
    class ( FieldSet_CGS_Form ), pointer :: &
      FieldSet_G => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type FieldSet_ASG_Form


contains


  subroutine Initialize &
               ( FSA, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption )

    class ( FieldSet_ASG_Form ), intent ( inout ), target :: &
      FSA
    class ( Atlas_SCG_Form ), intent ( in ), target :: &
      A
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    !-- FIXME: This shouldn't be necessary, but for some reason GCC 10.1.0
    !          doesn't compile without it
    class ( Chart_GS_Form ), pointer :: &
      G_Pointer

    if ( FSA % Type  ==  '' ) &
      FSA % Type  =  'a FieldSet_ASG'

    call FSA % Initialize_H ( A, NameOption )

    allocate ( FieldSet_CGS_Form :: FSA % FieldSet_C ( 1 ) % Element )
    select type ( FSG  =>  FSA % FieldSet_C ( 1 ) % Element )
    class is ( FieldSet_CGS_Form )

    select type ( G  =>  A % Chart ( 1 ) % Element )
    class is ( Chart_GS_Form )

    !-- FIXME: See FIXME above
    G_Pointer  =>  G
    call FSG % Initialize &
           ( G_Pointer, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
             nFieldsOption )

    FSA % FieldSet_G  =>  FSG

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'FieldSet_ASG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
    end select !-- FSG

  end subroutine Initialize


  impure elemental subroutine Finalize ( FSA )

    type ( FieldSet_ASG_Form ), intent ( inout ) :: &
      FSA

    nullify ( FSA % FieldSet_G )

  end subroutine Finalize

end module FieldSet_ASG__Form
