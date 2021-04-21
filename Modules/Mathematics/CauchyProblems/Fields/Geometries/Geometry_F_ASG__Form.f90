module Geometry_F_ASG__Form

  !-- Geometry_Flat_AtlasSingleGrid__Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometry_F_GS__Form
  use Geometry_F_AH__Form

  implicit none
  private

  type, public, extends ( Geometry_F_AH_Form ) :: Geometry_F_ASG_Form
    class ( FieldSet_GS_Form ), pointer :: &
      FieldSet_G => null ( )
    class ( Geometry_F_GS_Form ), pointer :: &
      Geometry_G => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Geometry_F_ASG_Form


contains


  subroutine Initialize &
               ( GA, A, NameOption, nFieldsOption, FieldOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption )

    class ( Geometry_F_ASG_Form ), intent ( inout ), target :: &
      GA
    class ( Atlas_SG_Form ), intent ( inout ), target :: &
      A
    character ( * ), intent ( inout ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    if ( GA % Type  ==  '' ) &
      GA % Type  =  'a Geometry_F_ASG'

    call GA % Initialize_H ( A, NameOption )

    allocate ( Geometry_F_GS_Form :: GA % Geometry_C ( 1 ) % Element )
    select type ( GG  =>  GA % Geometry_C ( 1 ) % Element )
    class is ( Geometry_F_GS_Form )

    select type ( G  =>  A % Chart ( 1 ) % Element )
    class is ( Grid_S_Form )

    call GG % Initialize &
           ( G, NameOption, nFieldsOption, FieldOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, UnitOption )

    GA % Geometry_G  =>  GG
    
    select type ( FSG  =>  GG % FieldSet )
    class is ( FieldSet_GS_Form )
      GA % FieldSet_G  =>  FSG
    end select

    call GA % Compute ( )

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Geometry_F_ASG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
    end select !-- GG

  end subroutine Initialize


  subroutine Compute ( GA )

    class ( Geometry_F_ASG_Form ), intent ( inout ) :: &
      GA

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, GA % Atlas % nCharts
      select type ( GC  =>  GA % Geometry_C ( iC ) % Element )
      class is ( Geometry_F_GS_Form )
      call GC % Compute ( )
      end select !-- GC
    end do !-- iC

  end subroutine Compute


  impure elemental subroutine Finalize ( GA )

    type ( Geometry_F_ASG_Form ), intent ( inout ) :: &
      GA

    nullify ( GA % Geometry_G )
    nullify ( GA % FieldSet_G )

  end subroutine Finalize

end module Geometry_F_ASG__Form
