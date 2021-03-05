!-- FieldHeader_C__Form handles metadata for a set of related fields on 
!   a Chart.

module FieldHeader_C__Form

  !-- FieldHeader_Chart__Form

  use Basics
  use ManifoldBasics
  use ChartHeader_Form

  implicit none
  private

  type, public :: FieldHeader_C_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
    logical ( KDL ) :: &
      Pinned
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameShort = ''
    class ( ChartHeaderForm ), pointer :: &
      Chart => null ( )
    class ( FieldHeader_M_Form ), pointer :: &
      Field_M => null ( )
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    final :: &
      Finalize
  end type FieldHeader_C_Form


contains


  subroutine Initialize_H &
               ( FC, FM, C, NameShort, PinnedOption, &
                 IgnorabilityOption )

    class ( FieldHeader_C_Form ), intent ( inout ) :: &
      FC
    class ( FieldHeader_M_Form ), intent ( in ), target :: &
      FM
    class ( ChartHeaderForm ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      PinnedOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FC % IGNORABILITY = C % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FC % IGNORABILITY = IgnorabilityOption

    if ( FC % Type == '' ) &
      FC % Type = 'a Field_C' 
    
    FC % Pinned = .false.
    if ( present ( PinnedOption ) ) &
      FC % Pinned = PinnedOption
    
    FC % Name = trim ( NameShort ) // '_' // trim ( C % Name ) 

    call Show ( 'Initializing ' // trim ( FC % Type ), FC % IGNORABILITY )
    call Show ( FC % Name, 'Name', FC % IGNORABILITY )
   
    FC % NameShort = NameShort
    call Show ( FC % NameShort, 'NameShort', FC % IGNORABILITY )

    FC % Chart    =>   C
    FC % Field_M  =>  FM 

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( FC )

    type ( FieldHeader_C_Form ), intent ( inout ) :: &
      FC

    nullify ( FC % Field_M )
    nullify ( FC % Chart )

    call Show ( 'Finalizing ' // trim ( FC % Type ), FC % IGNORABILITY )
    call Show ( FC % Name, 'Name', FC % IGNORABILITY )
   
  end subroutine Finalize


end module FieldHeader_C__Form
