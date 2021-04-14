module FieldSet_AH__Form

  !-- FieldSet_AtlasHeader__Form

  use Basics
  use Manifolds
  use FieldSet_CH__Form

  implicit none
  private

  type, public :: FieldSet_AH_Form
    integer ( KDI ) :: &
      IGNORABILITY
    character ( LDL ) :: &
      Type = '', &
      Name
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    type ( FieldSet_C_Element ), dimension ( : ), allocatable :: &
      FieldSet_C
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, private, pass :: &
      Show_FSA
    generic, public :: &
      Show => Show_FSA
    final :: &
      Finalize
  end type FieldSet_AH_Form


contains


  subroutine Initialize_H ( FSA, A, NameOption )

    class ( FieldSet_AH_Form ), intent ( inout ) :: &
      FSA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOption

    FSA % IGNORABILITY  =  A % IGNORABILITY

    if ( FSA % Type  ==  '' ) &
      FSA % Type  =  'a FieldSet_A'

    FSA % Name  =  'Fields'
    if ( present ( NameOption ) ) &
      FSA % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( FSA % Type ), A % IGNORABILITY )
    call Show ( FSA % Name, 'Name', A % IGNORABILITY )

    FSA % Atlas  =>  A

    allocate ( FSA % FieldSet_C ( A % nCharts ) )

  end subroutine Initialize_H


  subroutine Show_FSA ( FSA )

    class ( FieldSet_AH_Form ), intent ( in ) :: &
      FSA

   integer ( KDI ) :: &
     iC  !-- iC
   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( FSA % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FSA % IGNORABILITY )

    associate ( A  =>  FSA % Atlas )

    call Show ( FSA % Name, 'Name',  FSA % IGNORABILITY )
    call Show (   A % Name, 'Atlas', FSA % IGNORABILITY )

    do iC  =  1, A % nCharts
      if ( allocated ( FSA % FieldSet_C ( iC ) % Element ) ) then
        associate ( FSC  =>  FSA % FieldSet_C ( iC ) % Element )
        call FSC % Show ( )
        end associate !-- C
      end if  
    end do !-- iC

    end associate  !-- A

  end subroutine Show_FSA


  impure elemental subroutine Finalize ( FSA )

    type ( FieldSet_AH_Form ), intent ( inout ) :: &
      FSA

    if ( allocated ( FSA % FieldSet_C ) ) &
      deallocate ( FSA % FieldSet_C )  

    nullify ( FSA % Atlas )

    call Show ( 'Finalizing ' // trim ( FSA % Type ), FSA % IGNORABILITY )
    call Show ( FSA % Name, 'Name', FSA % IGNORABILITY )

  end subroutine Finalize


end module FieldSet_AH__Form
