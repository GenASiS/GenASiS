module FluxSet_Form

  use Basics
  use FieldSets

  implicit none
  private

  type, public :: FluxSetForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimer = 0
    character ( LDL ) :: &
      Type = '', &
      Name
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      Timer
    procedure, public, pass ( FS ) :: &
      Compute
    final :: &
      Finalize
  end type FluxSetForm

  type, public :: FluxSetElement
    class ( FluxSetForm ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type FluxSetElement


contains


  subroutine Initialize_H ( FS, NameOption, IgnorabilityOption )

    class ( FluxSetForm ), intent ( inout ) :: &
      FS
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FS % IGNORABILITY  =  CONSOLE % INFO_1
    if ( present ( IgnorabilityOption ) ) &
      FS % IGNORABILITY  =  IgnorabilityOption

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a FluxSet' 
    
    FS % Name  =  'Fluxes'
    if ( present ( NameOption ) ) &
      FS % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( FS % Type ), FS % IGNORABILITY )
    call Show ( FS % Name, 'Name', FS % IGNORABILITY )

  end subroutine Initialize_H


  subroutine Show_FS ( FS )

    class ( FluxSetForm ), intent ( in ) :: &
      FS

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( FS % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FS % IGNORABILITY )

    call Show ( FS % Name, 'Name',  FS % IGNORABILITY )

  end subroutine Show_FS


  function Timer ( FS, LevelOption ) result ( T )

    class ( FluxSetForm ), intent ( inout ) :: &
      FS
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  FS % iTimer )

    if ( iT == 0 ) then
      TimerName  =  FS % Name
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer


  subroutine Compute ( FS_FS, FS_SS, FS, FS_CS, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_FS, &  !-- FluxSet
      FS_SS     !-- StressSet
    class ( FluxSetForm ), intent ( in ) :: &
      FS
    class ( FieldSetForm ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

  end subroutine Compute


  impure elemental subroutine Finalize ( FS )

    type ( FluxSetForm ), intent ( inout ) :: &
      FS

    call Show ( 'Finalizing ' // trim ( FS % Type ), FS % IGNORABILITY )
    call Show ( FS % Name, 'Name', FS % IGNORABILITY )
   
  end subroutine Finalize


  impure elemental subroutine Finalize_E ( FSE )
    
    type ( FluxSetElement ), intent ( inout ) :: &
      FSE

    if ( allocated ( FSE % Element ) ) &
      deallocate ( FSE % Element )

  end subroutine Finalize_E


end module FluxSet_Form
