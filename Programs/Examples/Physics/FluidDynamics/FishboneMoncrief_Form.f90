module FishboneMoncrief_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_CE_Form ) :: FishboneMoncriefForm
  contains
    procedure, private, pass :: &
      Initialize_H
    final :: &
      Finalize
  end type FishboneMoncriefForm


contains


  subroutine Initialize_H ( U, NameOption )

    class ( FishboneMoncriefForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a FishboneMoncrief'

    Name  =  'FishboneMoncrief'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    ! call InitializeUniverse ( U, Name )
    ! call InitializeDiagnostics ( U )

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( FM )
    
    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

  end subroutine Finalize


end module FishboneMoncrief_Form
