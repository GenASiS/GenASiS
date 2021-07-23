module Universe_F_CC__Form

  !-- Universe_Fluid_CentralCore__Form

  use Basics
  use Mathematics
  use Universe_F_C__Form

  implicit none
  private

  type, public, extends ( Universe_F_C_Form ) :: Universe_F_CC_Form
    real ( KDR ) :: &
      GravityFactor
  contains
    procedure, private, pass :: &
      Initialize_F_CC
    generic, public :: &
      Initialize => Initialize_F_CC
    procedure, private, pass :: &
      InitializeAtlas
  end type Universe_F_CC_Form

contains


  subroutine Initialize_F_CC &
               ( U, FluidType, GravitationType, NameOption, FinishTimeOption, &
                 RadiusMaxOption, RadiusCoreOption, RadialRatioOption, &
                 nCellsPolarOption, nWriteOption )

    class ( Universe_F_CC_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType
    character ( * ), intent ( in ), optional :: &
      NameOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    if ( U % Type == '' ) &
      U % Type = 'a Universe_F_CC'

    call U % Initialize_F_C &
           ( FluidType, GravitationType, &
             NameOption = NameOption, &
             FinishTimeOption = FinishTimeOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption, &
             nWriteOption = nWriteOption ) 

  end subroutine Initialize_F_CC


  subroutine InitializeAtlas &
               ( U, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

      class ( Universe_F_CC_Form ), intent ( inout ) :: &
        U
      real ( KDR ), intent ( in ), optional :: &
        RadiusMaxOption, &
        RadiusCoreOption, &
        RadiusExcisionOption, &
        RadialRatioOption
      integer ( KDI ), intent ( in ), optional :: &
        nCellsPolarOption

    real ( KDR ) :: &
      RadiusMax, &
      RadiusCore, &
      RadialRatio

    associate ( I  =>  U % Integrator )

    allocate ( Atlas_SCG_CC_Form :: I % X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_CC_Form )

    allocate ( U % Units_F ( 1 ) )

    RadiusMax  =  10.0_KDR
    if ( present ( RadiusMaxOption ) ) &
      RadiusMax  =  RadiusMaxOption

    RadiusCore  =  10.0_KDR / 8.0_KDR
    if ( present ( RadiusCoreOption ) ) &
      RadiusCore  =  RadiusCoreOption

    call PS % Initialize &
           ( RadiusMax = RadiusMax, &
             RadiusCore = RadiusCore, &
             CommunicatorOption = PROGRAM_HEADER % Communicator, &
             NameOption = 'PositionSpace' )

    ! if ( FC % Dimensionless ) then

    !   call PS % CreateChart_CC ( )

    ! else

    !   RadiusCore   =   16.0_KDR  *  UNIT % KILOMETER
    !   RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
    !   RadialRatio  =  2.4_KDR

    !   call PS % CreateChart_CC &
    !          ( CoordinateUnitOption = FC % Units % Coordinate_PS, &
    !            RadiusCoreOption = RadiusCore, &
    !            RadiusMaxOption = RadiusMax, &
    !            RadialRatioOption = RadialRatio, &
    !            nCellsPolarOption = nCellsPolarOption )

    !   FC % RadiusPolarMomentum  =  8.0_KDR  *  UNIT % KILOMETER
    !   call PROGRAM_HEADER % GetParameter &
    !          ( FC % RadiusPolarMomentum, 'RadiusPolarMomentum' )

    ! end if !-- Dimensionless

    end select !-- PS
    end associate !-- I

  end subroutine InitializeAtlas


end module Universe_F_CC__Form
