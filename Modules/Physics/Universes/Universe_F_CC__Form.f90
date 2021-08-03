module Universe_F_CC__Form

  !-- Universe_Fluid_CentralCore__Form

  use Basics
  use Mathematics
  use Gravitations
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
    final :: &
      Finalize
    procedure, public, pass :: &
      SetBoundaryConditions
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
    I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
    end associate !-- I

  end subroutine Initialize_F_CC


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_CC_Form ), intent ( inout ) :: &
      U

  end subroutine Finalize


  subroutine SetBoundaryConditions ( U )

    class ( Universe_F_CC_Form ), intent ( inout ) :: &
      U
    
    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    associate &
      ( F  =>  I % CurrentSet_X )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
    call F % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    end associate !-- F

    end select !-- I

  end subroutine SetBoundaryConditions


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

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      oC, &  !-- oCell
      oI     !-- oIncoming
    !   iRadius
    real ( KDR ) :: &
      Constant_G
    !   VelocityMax, &
    !   VelocityMaxRadius, &
    !   NumberDensityAve, &
    !   NumberDensityMax, &
    !   TimeScaleVelocity, &
    !   TimeScaleDensity
    real ( KDR ), dimension ( : ), allocatable :: &
      Radius, &
      Density, &
      Velocity
    real ( KDR ), dimension ( :, : ), pointer :: &
      Outgoing_2D, &
      Incoming_2D
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO

    select type ( U  =>  I % System )
      class is ( Universe_F_C_Form )
    select type ( F_SA  =>  U % SA_Fluid % FieldSet_SA )
      class is ( Fluid_D_Form )
    select type ( G_SA  =>  U % SA_Gravitation % FieldSet_SA )
      class is ( Gravitation_G_Form )
    select type ( A_SA  =>  F_SA % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C_SA    =>  A_SA % Chart_GS, &
        G_SA_V  =>  G_SA % Storage_GS % Value, &
        F_SA_V  =>  F_SA % Storage_GS % Value )
    associate &
      ( nGL  =>  C_SA % nGhostLayers ( 1 ), &
        nC   =>  C_SA % nCells ( 1 ), &
        nCB  =>  C_SA % nCellsBrick ( 1 ), &
        nP   =>  C_SA % Communicator % Size, &
        nF   =>  3, &
         R   =>  G_SA_V ( :, G_SA % CENTER_U_1 ), &
         M   =>  F_SA_V ( :, F_SA % BARYON_MASS ), &
         N   =>  F_SA_V ( :, F_SA % BARYON_DENSITY_C ), &
         V   =>  F_SA_V ( :, F_SA % VELOCITY_U_1 ) )

    !-- Gather spherically averaged density and velocity
    !   (assume decomposition in spherical shells)

    allocate ( CO )
    call CO % Initialize &
           ( C_SA % Communicator, &
             nOutgoing = [ nF * nCB ], nIncoming = [ nF * nC ] )

    oC  =  C_SA % nGhostLayers ( 1 )
    Outgoing_2D ( 1 : nCB, 1 : 3 )  &
      =>  CO % Outgoing % Value
    Outgoing_2D ( 1 : nCB, 1 )  &
      =   R ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 2 )  &
      =   M ( nGL + 1 : nGL + nCB )  *  N ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 3 )  &
      =   V ( nGL + 1 : nGL + nCB )

call Show ( Outgoing_2D, '>>> Outgoing_2D' )

    call CO % Gather ( )

    allocate ( Radius ( nC ), Density ( nC ), Velocity ( nC ) )
    do iP  =  0,  nP - 1
      oC  =  iP * nCB
      oI  =  oC * nF
      Incoming_2D ( 1 : nCB, 1 : nF )  &
        =>  CO % Incoming % Value ( oI + 1 : oI + nCB * nF ) 
      Radius   ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 1 )
      Density  ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 2 )
      Velocity ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 3 )
    end do !-- iP

call Show ( Radius,  '>>> Radius' )
call Show ( Density, '>>> Density' )
call Show ( Velocity, '>>> Velocity' )

    !-- FIXME: nontrivial units
    Constant_G  =  1.0_KDR

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

    deallocate ( Velocity, Density, Radius )
    deallocate ( CO )

    end associate !-- nC, etc.
    end associate !-- C_SA
    end select !-- A_SA
    end select !-- G_SA
    end select !-- F_SA
    end select !-- U

  end subroutine Set_T_CheckpointInterval


end module Universe_F_CC__Form
