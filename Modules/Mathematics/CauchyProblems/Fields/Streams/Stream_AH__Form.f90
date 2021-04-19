module Stream_AH__Form

  !-- Stream_AtlasHeader__Form

  use Basics
  use Manifolds
  use FieldSets
  use Stream_CH__Form

  implicit none
  private

  type, public :: Stream_AH_Form
    integer ( KDI ) :: &
      IGNORABILITY
    character ( LDL ) :: &
      Type = '', &
      Name
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    type ( Stream_C_Element ), dimension ( : ), allocatable :: &
      Stream_C
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, public, pass :: &
      AddFieldSet
  !   procedure, private, pass :: &
  !     Show_FSA
  !   generic, public :: &
  !     Show => Show_FSA
    final :: &
      Finalize
  end type Stream_AH_Form


contains


  subroutine Initialize_H ( SA, A, NameOption )

    class ( Stream_AH_Form ), intent ( inout ) :: &
      SA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOption

    SA % IGNORABILITY  =  A % IGNORABILITY

    if ( SA % Type  ==  '' ) &
      SA % Type  =  'a Stream_A'

    SA % Name  =  'Stream'
    if ( present ( NameOption ) ) &
      SA % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( SA % Type ), A % IGNORABILITY )
    call Show ( SA % Name, 'Name', A % IGNORABILITY )

    SA % Atlas  =>  A

    allocate ( SA % Stream_C ( A % nCharts ) )

  end subroutine Initialize_H


  subroutine AddFieldSet ( SA, FSA, NameOption, iaSelectedOption )

    class ( Stream_AH_Form ), intent ( inout ) :: &
      SA
    class ( FieldSet_AH_Form ), intent ( in ) :: &
      FSA
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, SA % Atlas % nCharts
      associate &
        (  SC  =>   SA %   Stream_C ( iC ) % Element, &
          FSC  =>  FSA % FieldSet_C ( iC ) % Element )
      call SC % AddFieldSet ( FSC, NameOption, iaSelectedOption )
      end associate !-- SC, etc.
    end do !-- iC

  end subroutine AddFieldSet


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
