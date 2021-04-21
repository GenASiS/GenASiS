module Geometry_F_AH__Form

  !-- Geometry_Flat_AtlasHeader_Form

  use Basics
  use Manifolds
  use Streams
  use Geometry_F_CH__Form

  implicit none
  private

  type, public :: Geometry_F_AH_Form
    integer ( KDI ) :: &
      IGNORABILITY
    character ( LDL ) :: &
      Type = '', &
      Name
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    type ( Geometry_C_Element ), dimension ( : ), allocatable :: &
      Geometry_C
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, public, pass ( GA ) :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_GA
    final :: &
      Finalize
  end type Geometry_F_AH_Form


contains


  subroutine Initialize_H ( GA, A, NameOption )

    class ( Geometry_F_AH_Form ), intent ( inout ) :: &
      GA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOption

    GA % IGNORABILITY  =  A % IGNORABILITY

    if ( GA % Type  ==  '' ) &
      GA % Type  =  'a Geometry_F_A'

    GA % Name  =  'Geometry'
    if ( present ( NameOption ) ) &
      GA % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( GA % Type ), A % IGNORABILITY )
    call Show ( GA % Name, 'Name', A % IGNORABILITY )

    GA % Atlas  =>  A

    allocate ( GA % Geometry_C ( A % nCharts ) )

  end subroutine Initialize_H


  subroutine SetStream ( SA, GA )

    class ( Stream_AH_Form ), intent ( inout ) :: &
      SA
    class ( Geometry_F_AH_Form ), intent ( in ) :: &
      GA

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, GA % Atlas % nCharts
      associate &
        ( SC  =>  SA %   Stream_C ( iC ) % Element, &
          GC  =>  GA % Geometry_C ( iC ) % Element )
      call GC % SetStream ( SC )
      end associate !-- SC, etc.
    end do !-- iC

  end subroutine SetStream


  subroutine Show_GA ( GA )

    class ( Geometry_F_AH_Form ), intent ( inout ) :: &
      GA

   integer ( KDI ) :: &
     iC  !-- iC
   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( GA % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', GA % IGNORABILITY )

    associate ( A  =>  GA % Atlas )

    call Show ( GA % Name, 'Name',  GA % IGNORABILITY )
    call Show (  A % Name, 'Atlas', GA % IGNORABILITY )

    do iC  =  1, A % nCharts
      if ( allocated ( GA % Geometry_C ( iC ) % Element ) ) then
        associate ( GC  =>  GA % Geometry_C ( iC ) % Element )
        call GC % Show ( )
        end associate !-- GC
      end if  
    end do !-- iC

    end associate  !-- A

  end subroutine Show_GA


  impure elemental subroutine Finalize ( GA )

    type ( Geometry_F_AH_Form ), intent ( inout ) :: &
      GA

    if ( allocated ( GA % Geometry_C ) ) &
      deallocate ( GA % Geometry_C )  

    nullify ( GA % Atlas )

    call Show ( 'Finalizing ' // trim ( GA % Type ), GA % IGNORABILITY )
    call Show ( GA % Name, 'Name', GA % IGNORABILITY )

  end subroutine Finalize


end module Geometry_F_AH__Form
