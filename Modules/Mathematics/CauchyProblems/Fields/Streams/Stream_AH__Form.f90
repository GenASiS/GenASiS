module Stream_AH__Form

  !-- Stream_AtlasHeader__Form

  use Basics
  use Manifolds
  use Stream_CH__Form

  implicit none
  private

  type, public :: Stream_AH_Form
    integer ( KDI ) :: &
      IGNORABILITY
    logical ( KDL ) :: &
      Verbose = .false.
    character ( LDL ) :: &
      Type = '', &
      Name
    type ( GridImageStreamForm ), pointer :: &
      GridImageStream => null ( )
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    type ( Stream_C_Element ), dimension ( : ), allocatable :: &
      Stream_C
  contains
    procedure, public, pass :: &
      Initialize_H
  !   procedure, private, pass :: &
  !     Show_FSA
  !   generic, public :: &
  !     Show => Show_FSA
    final :: &
      Finalize
  end type Stream_AH_Form


contains


  subroutine Initialize_H ( SA, A, GIS, NameOption, VerboseOption )

    class ( Stream_AH_Form ), intent ( inout ) :: &
      SA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    SA % IGNORABILITY  =  A % IGNORABILITY

    if ( SA % Type  ==  '' ) &
      SA % Type  =  'a Stream_A'

    SA % Name  =  'Stream'
    if ( present ( NameOption ) ) &
      SA % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( SA % Type ), A % IGNORABILITY )
    call Show ( SA % Name, 'Name', A % IGNORABILITY )

    SA % Verbose  =  .false.
    if ( present ( VerboseOption ) ) &
      SA % Verbose  =  VerboseOption

    SA % GridImageStream  =>  GIS
    SA % Atlas            =>  A

    allocate ( SA % Stream_C ( A % nCharts ) )

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( SA )

    type ( Stream_AH_Form ), intent ( inout ) :: &
      SA

    if ( allocated ( SA % Stream_C ) ) &
      deallocate ( SA % Stream_C )  

    nullify ( SA % Atlas )

    call Show ( 'Finalizing ' // trim ( SA % Type ), SA % IGNORABILITY )
    call Show ( SA % Name, 'Name', SA % IGNORABILITY )

  end subroutine Finalize


end module Stream_AH__Form
