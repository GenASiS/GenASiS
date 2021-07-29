module Universe_F_CC__Form

  !-- Universe_Fluid_CentralCore__Form

  use Basics
  use Mathematics
  use Fluids
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

      private :: &
        Set_T_CheckpointInterval

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

    associate ( I  =>  U % Integrator )
    I % Set_T_CheckpointInterval => Set_T_CheckpointInterval
    end associate !-- I

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

    allocate ( Atlas_SCG_CC_Form :: U % PositionSpace_SA )
    select type ( PS_SA  =>  U % PositionSpace_SA )
      class is ( Atlas_SCG_CC_Form )
    call PS_SA % Initialize ( PS )
    end select !-- PS_SA

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


  subroutine Set_T_CheckpointInterval ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    ! integer ( KDI ) :: &
    !   iRadius
    ! real ( KDR ) :: &
    !   GravitationalConstant, &
    !   BaryonMass, &
    !   VelocityMax, &
    !   VelocityMaxRadius, &
    !   NumberDensityAve, &
    !   NumberDensityMax, &
    !   TimeScaleVelocity, &
    !   TimeScaleDensity
    ! type ( CollectiveOperation_R_Form ), allocatable :: &
    !   CO
    ! class ( GeometryFlatForm ), pointer :: &
    !   G_SA
    ! class ( Fluid_D_Form ), pointer :: &
    !   F_SA

    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( F => I % CurrentSet_X )
      class is ( Fluid_D_Form )

    ! select type ( FC => I % Universe )
    ! class is ( FluidCentralCoreForm )

    ! GravitationalConstant  =  CONSTANT % GRAVITATIONAL
    !            BaryonMass  =  CONSTANT % ATOMIC_MASS_UNIT 
    ! if ( FC % Dimensionless ) then
    !   GravitationalConstant  =  1.0_KDR
    !              BaryonMass  =  1.0_KDR
    ! end if

    ! G_SA  =>  FC % PositionSpace_SA % Geometry ( )
    ! F_SA  =>  FC % Fluid_ASC_SA % Fluid_D ( )

    ! !-- Velocity

    ! iRadius &
    !   =  maxloc ( abs ( F_SA % Value ( :,  F_SA % VELOCITY_U ( 1 ) ) ), &
    !               dim = 1 )
    ! VelocityMaxRadius &
    !   =  G_SA % Value ( iRadius, G_SA % CENTER_U ( 1 ) )
    ! VelocityMax  &
    !   =  max ( maxval ( abs ( F_SA % Value ( :, F_SA % VELOCITY_U ( 1 ) ) ) ), &
    !            sqrt ( tiny ( 0.0_KDR ) ) )  

    ! !-- Density

    ! select type ( TI => FA % TallyInterior )
    ! class is ( Tally_F_D_Form )
    !   NumberDensityAve  =  TI % Value ( TI % BARYON_NUMBER ) &
    !                  / ( 4.0_KDR / 3.0_KDR  *  CONSTANT % PI  &
    !                      *  VelocityMaxRadius ** 3 )
    ! end select !-- TI

    ! NumberDensityMax  &
    !   =  max ( maxval ( F_SA % Value ( :, F_SA % COMOVING_BARYON_DENSITY ) ), &
    !            sqrt ( tiny ( 0.0_KDR ) ) )

    ! !-- Time scales and CheckpointTimeInterval

    ! TimeScaleVelocity &
    !   =  VelocityMaxRadius  /  VelocityMax
    ! TimeScaleDensity &
    !   =  ( GravitationalConstant * BaryonMass &
    !        * min ( NumberDensityAve, NumberDensityMax ) ) &
    !      ** ( -0.5_KDR )

    ! allocate ( CO )
    ! associate ( C => I % PositionSpace % Communicator ) 
    ! call CO % Initialize &
    !        ( C, nOutgoing = [ 1 ], nIncoming = [ 1 ], &
    !          RootOption = CONSOLE % DisplayRank )
    ! end associate !-- C

    ! CO % Outgoing % Value ( 1 )  &
    !   =  min ( TimeScaleVelocity, TimeScaleDensity )  /  I % nWrite

    ! call CO % Broadcast ( )

    ! I % CheckpointTimeInterval  =  CO % Incoming % Value ( 1 )

    ! !-- Display

    ! call Show ( 'Time Scales', I % IGNORABILITY )
    ! call Show ( VelocityMaxRadius, FC % Units % Coordinate_PS ( 1 ), &
    !             'VelocityMaxRadius', I % IGNORABILITY )
    ! call Show ( VelocityMax, FC % Units % Velocity_U ( 1 ), &
    !             'VelocityMax', I % IGNORABILITY )
    ! call Show ( NumberDensityAve, FC % Units % NumberDensity, &
    !             'NumberDensityAve', I % IGNORABILITY )
    ! call Show ( NumberDensityMax, FC % Units % NumberDensity, &
    !             'NumberDensityMax', I % IGNORABILITY )
    ! call Show ( BaryonMass * NumberDensityAve, FC % Units % MassDensity, &
    !             'MassDensityAve', I % IGNORABILITY )
    ! call Show ( BaryonMass * NumberDensityMax, FC % Units % MassDensity, &
    !             'MassDensityMax', I % IGNORABILITY )
    ! call Show ( TimeScaleVelocity, I % TimeUnit, &
    !             'TimeScaleVelocity', I % IGNORABILITY )
    ! call Show ( TimeScaleDensity, I % TimeUnit, &
    !             'TimeScaleDensity', I % IGNORABILITY )

    !-- Cleanup

    ! end select !-- FC
    end select !-- F
    end select !-- I

    ! nullify ( G_SA, F_SA )

  end subroutine Set_T_CheckpointInterval


end module Universe_F_CC__Form
