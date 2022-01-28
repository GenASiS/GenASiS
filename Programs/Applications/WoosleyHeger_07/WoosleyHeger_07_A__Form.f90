module WoosleyHeger_07_A__Form

  !-- WoosleyHeger_07_Adiabatic__Form

  use GenASiS
  use WoosleyHeger_07__Form

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Form ) :: WoosleyHeger_07_A_Form
  contains
    procedure, private, pass :: &
      Initialize_H
!     procedure, private, pass :: &
!       Initialize_WH
!     generic, public :: &
!       Initialize => Initialize_WH
    final :: &
      Finalize
  end type WoosleyHeger_07_A_Form

    private :: &
      InitializeUniverse, &
      SetInitial

!     private :: &
!       InitializeRadiationCentralCore, &
!       SetProblem

contains


  subroutine Initialize_H ( U, NameOption )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a WoosleyHeger_07_A'

    Name  =  'WoosleyHeger_07_A'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeUniverse ( U, Name )

  end subroutine Initialize_H


!   subroutine Initialize_WH ( WH, Name )

!     class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
!       WH
!     character ( * ), intent ( in )  :: &
!       Name

!     if ( WH % Type == '' ) &
!       WH % Type = 'a WoosleyHeger_07_A'

!     call InitializeRadiationCentralCore ( WH, Name )
!     call SetProblem ( WH )

!   end subroutine Initialize_WH


  subroutine Finalize ( WH )

    type ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      WH

  end subroutine Finalize


  subroutine InitializeUniverse ( WH, Name )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ) :: &
      FinishTime

    FinishTime  =  1.0_KDR  *  UNIT % SECOND

    call WH % Initialize &
           ( FluidType = 'HEAVY_NUCLEUS', &
             GravitationType = 'NEWTON_SG', &
             NameOption = Name, &
             FinishTimeOption = FinishTime, &
             nCellsPolarOption = 128, &
             nWriteOption = 30 )

    WH % Integrator % SetInitial  =>  SetInitial
    WH % Integrator % System      =>  WH

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( WH  =>  I % System )
      class is ( WoosleyHeger_07_Form )

    call WH % SetFluid ( )

    end select !-- WH

  end subroutine SetInitial


!   subroutine InitializeRadiationCentralCore ( WH, Name )

!     class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
!       WH
!     character ( * ), intent ( in )  :: &
!       Name

!     logical ( KDL ) :: &
!       UseDevice
!     character ( LDL ) :: &
!       GeometryType

!     GeometryType = 'NEWTONIAN'
!     call PROGRAM_HEADER % GetParameter ( GeometryType, 'GeometryType' )
    
!     UseDevice = ( OffloadEnabled ( ) .and. GetNumberOfDevices ( ) >= 1 )
!     call PROGRAM_HEADER % GetParameter ( UseDevice, 'UseDevice' )

!     call WH % Initialize &
!            ( RadiationName = [ 'None' ], RadiationType = [ 'NONE' ], &
!              MomentsType = 'NONE', FluidType = 'HEAVY_NUCLEUS', &
!              GeometryType = GeometryType, Name = Name, &
!              ShockThresholdOption = 1.0_KDR, nWriteOption = 30,  &
!              RadiationUseDeviceOption = UseDevice, &
!              FluidUseDeviceOption = UseDevice, &
!              GeometryUseDeviceOption = UseDevice )

!   end subroutine InitializeRadiationCentralCore


!   subroutine SetProblem ( WH )

!     class ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
!       WH

!     call WH % SetFluid ( )

!   end subroutine SetProblem


end module WoosleyHeger_07_A__Form
