module WoosleyHeger_07_A__Form

  !-- WoosleyHeger_07_Adiabatic__Form

  use GenASiS
  use WoosleyHeger_07__Template

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Template ) :: WoosleyHeger_07_A_Form
  contains
    procedure, private, pass :: &
      Initialize_WH
    generic, public :: &
      Initialize => Initialize_WH
    final :: &
      Finalize
  end type WoosleyHeger_07_A_Form

    private :: &
      InitializeRadiationCentralCore, &
      SetProblem

contains


  subroutine Initialize_WH ( WH, Name )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    if ( WH % Type == '' ) &
      WH % Type = 'a WoosleyHeger_07_A'

    call InitializeRadiationCentralCore ( WH, Name )
    call SetProblem ( WH )

  end subroutine Initialize_WH


  subroutine Finalize ( WH )

    type ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      WH

  end subroutine Finalize


  subroutine InitializeRadiationCentralCore ( WH, Name )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    logical ( KDL ) :: &
      UseDevice
    character ( LDL ) :: &
      GeometryType

    GeometryType = 'NEWTONIAN'
    call PROGRAM_HEADER % GetParameter ( GeometryType, 'GeometryType' )
    
    UseDevice = ( OffloadEnabled ( ) .and. GetNumberOfDevices ( ) >= 1 )
    call PROGRAM_HEADER % GetParameter ( UseDevice, 'UseDevice' )

    call WH % Initialize &
           ( RadiationName = [ 'None' ], RadiationType = [ 'NONE' ], &
             MomentsType = 'NONE', FluidType = 'HEAVY_NUCLEUS', &
             GeometryType = GeometryType, Name = Name, &
             ShockThresholdOption = 1.0_KDR, nWriteOption = 30,  &
             RadiationUseDeviceOption = UseDevice, &
             FluidUseDeviceOption = UseDevice, &
             GeometryUseDeviceOption = UseDevice )

  end subroutine InitializeRadiationCentralCore


  subroutine SetProblem ( WH )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      WH

    call WH % SetFluid ( )

  end subroutine SetProblem


end module WoosleyHeger_07_A__Form
