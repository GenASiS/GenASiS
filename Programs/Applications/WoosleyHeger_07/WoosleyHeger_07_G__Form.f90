module WoosleyHeger_07_G__Form

  !-- WoosleyHeger_07_Grey__Form

  use GenASiS
  use WoosleyHeger_07__Template

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Template ) :: WoosleyHeger_07_G_Form
  contains
    procedure, private, pass :: &
      Initialize_WH
    generic, public :: &
      Initialize => Initialize_WH
  end type WoosleyHeger_07_G_Form

!    private :: &
!      InitializeRadiationCentralCore

contains


  subroutine Initialize_WH ( WH, Name )

    class ( WoosleyHeger_07_G_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    if ( WH % Type == '' ) &
      WH % Type = 'a WoosleyHeger_07_G'

!    call InitializeRadiationCentralCore ( WH, Name )
    call WH % SetFluid ( )

  end subroutine Initialize_WH


  ! subroutine InitializeRadiationCentralCore ( WH, Name )

  !   class ( WoosleyHeger_07_G_Form ), intent ( inout ), target :: &
  !     WH
  !   character ( * ), intent ( in )  :: &
  !     Name

  !   character ( LDL ) :: &
  !     GeometryType

  !   GeometryType = 'NEWTONIAN'
  !   call PROGRAM_HEADER % GetParameter ( GeometryType, 'GeometryType' )

  !   call WH % Initialize &
  !          ( RadiationName = [ 'None' ], RadiationType = [ 'NONE' ], &
  !            MomentsType = 'NONE', FluidType = 'HEAVY_NUCLEUS', &
  !            GeometryType = GeometryType, Name = Name, &
  !            ShockThresholdOption = 1.0_KDR, nWriteOption = 30 )

  ! end subroutine InitializeRadiationCentralCore


end module WoosleyHeger_07_G__Form
