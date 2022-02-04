module DivergenceContribution_Form

  use Basics
  use FieldSets

  implicit none
  private

  type, public :: DivergenceContributionForm
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
      Show => Show_DC
    procedure, public, pass :: &
      Timer
    procedure, public, pass ( DC ) :: &
      ComputeFluxes
    final :: &
      Finalize
  end type DivergenceContributionForm

  type, public :: DivergenceContributionElement
    class ( DivergenceContributionForm ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type DivergenceContributionElement


contains


  subroutine Initialize_H ( DC, NameOption, IgnorabilityOption )

    class ( DivergenceContributionForm ), intent ( inout ) :: &
      DC
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    DC % IGNORABILITY  =  CONSOLE % INFO_1
    if ( present ( IgnorabilityOption ) ) &
      DC % IGNORABILITY  =  IgnorabilityOption

    if ( DC % Type  ==  '' ) &
      DC % Type  =  'a DivergenceContribution' 
    
    DC % Name  =  'DivergenceContribution'
    if ( present ( NameOption ) ) &
      DC % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( DC % Type ), DC % IGNORABILITY )
    call Show ( DC % Name, 'Name', DC % IGNORABILITY )

  end subroutine Initialize_H


  subroutine Show_DC ( DC )

    class ( DivergenceContributionForm ), intent ( in ) :: &
      DC

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( DC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', DC % IGNORABILITY )

    call Show ( DC % Name, 'Name',  DC % IGNORABILITY )

  end subroutine Show_DC


  function Timer ( DC, LevelOption ) result ( T )

    class ( DivergenceContributionForm ), intent ( inout ) :: &
      DC
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  DC % iTimer )

    if ( iT == 0 ) then
      TimerName  =  DC % Name
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer


  subroutine ComputeFluxes ( FS_F, DC, FS_CS, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_F  !-- DivergenceContribution
    class ( DivergenceContributionForm ), intent ( in ) :: &
      DC
    class ( FieldSetForm ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

  end subroutine ComputeFluxes


  impure elemental subroutine Finalize ( DC )

    type ( DivergenceContributionForm ), intent ( inout ) :: &
      DC

    call Show ( 'Finalizing ' // trim ( DC % Type ), DC % IGNORABILITY )
    call Show ( DC % Name, 'Name', DC % IGNORABILITY )
   
  end subroutine Finalize


  impure elemental subroutine Finalize_E ( DCE )
    
    type ( DivergenceContributionElement ), intent ( inout ) :: &
      DCE

    if ( allocated ( DCE % Element ) ) &
      deallocate ( DCE % Element )

  end subroutine Finalize_E


end module DivergenceContribution_Form
