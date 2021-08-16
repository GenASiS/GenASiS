module Universe_F_C__Form

  !-- Universe_Fluid_Central__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_H__Form

  implicit none
  private

  type, public, extends ( Universe_H_Form ) :: Universe_F_C_Form
    real ( KDR ) :: &
      GravityFactor = 0.0_KDR
    class ( Atlas_SCG_Form ), allocatable :: &
      PositionSpace_SA  !-- SphericalAverage
    type ( StreamForm ), allocatable :: &
      Stream_SA
    type ( SphericalAverageForm ), allocatable :: &
      SA_Gravitation, &
      SA_Fluid
    type ( Units_F_Form ), dimension ( : ), allocatable :: &
      Units_F
  contains
    procedure, public, pass :: &
      Initialize_F_C
    procedure, public, pass :: &
      Show => Show_U
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_F_C
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_F_C
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeGravitation
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SetBoundaryConditions
    procedure, public, pass :: &
      InitializeStep
    procedure, private, pass :: &
      InitializeAtlas
    procedure, public, nopass :: &
      Analyze_C
    procedure, public, nopass :: &
      Write_C
  end type Universe_F_C_Form

    private :: &
      SetSlope_N_SG


contains


  subroutine Initialize_F_C &
               ( U, FluidType, GravitationType, NameOption, FinishTimeOption, &
                 RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption, nWriteOption )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
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
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    character ( LDL ) :: &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a Universe_F_C'

    Name  =  'Universe'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call U % Universe_H_Form % Initialize ( NameOption = Name )

    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadiusExcisionOption = RadiusExcisionOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )
    call U % InitializeGravitation &
           ( GravitationType )
    call U % InitializeFluid &
           ( FluidType )
    call U % SetBoundaryConditions &
           ( )
    call U % InitializeStep &
           ( )

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    if ( .not. allocated ( I % dT_Label ) ) then
      if ( any ( trim ( GravitationType ) == [ 'NEWTON_SG' ] ) ) then
        allocate ( I % dT_Label ( 2 ) )
        I % dT_Label ( 1 )  =  'Fluid advection'
        I % dT_Label ( 2 )  =  'Gravitation acceleration'
        U % GravityFactor = 0.7_KDR
    !    if ( present ( GravityFactorOption ) ) &
    !      U % GravityFactor = GravityFactorOption
        call PROGRAM_HEADER % GetParameter &
               ( U % GravityFactor, 'GravityFactor' )
      else
        allocate ( I % dT_Label ( 1 ) )
        I % dT_Label ( 1 )  =  'Fluid advection'
      end if
    end if

    call I % Initialize &
           ( Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
!             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

    !-- SphericalAverage Stream

    allocate ( U % Stream_SA )
    associate &
      (   A_SA  =>  U % PositionSpace_SA, &
          S_SA  =>  U % Stream_SA, &
        GIS     =>  I % GridImageStream, &
          S     =>  I % Checkpoint_X )

    call S_SA % Initialize &
           ( A_SA, GIS, NameOption = trim ( S % Name ) // '_SA' )

    select type ( G_SA  =>  U % SA_Gravitation % FieldSet_SA )
      class is ( Geometry_F_Form )
    select type ( F_SA  =>  U % SA_Fluid % FieldSet_SA )
      class is ( Fluid_D_Form )
    call G_SA % SetStream ( S_SA )
    call F_SA % SetStream ( S_SA )
    end select !-- F_SA
    end select !-- G_SA

    end associate !-- A_SA, etc.

    I % Analyze  =>  Analyze_C
    I % Write    =>  Write_C

    end select !-- I

  end subroutine Initialize_F_C


  subroutine Show_U ( U )

    class ( Universe_F_C_Form ), intent ( in ) :: &
      U

    call U % Universe_H_Form % Show ( )

    call Show ( 'Universe_F_C Proper Parameters' )

    if ( U % GravityFactor  >  0.0_KDR ) &
      call Show ( U % GravityFactor, 'GravityFactor' )

    call U % PositionSpace_SA % Show ( )
    call U % SA_Gravitation % FieldSet_SA % Show ( )
    call U % SA_Fluid % FieldSet_SA % Show ( )
    call U % Stream_SA % Show ( )

  end subroutine Show_U


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_F ) ) &
      deallocate ( U % Units_F )
    if ( allocated ( U % SA_Fluid ) ) &
      deallocate ( U % SA_Fluid )
    if ( allocated ( U % SA_Gravitation ) ) &
      deallocate ( U % SA_Gravitation )
    if ( allocated ( U % Stream_SA ) ) &
      deallocate ( U % Stream_SA )
    if ( allocated ( U % PositionSpace_SA ) ) &
      deallocate ( U % PositionSpace_SA )
    
  end subroutine Finalize


  subroutine AllocateIntegrator_F_C ( U )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    allocate ( Integrator_CS_Form :: U % Integrator )

    if ( allocated ( U % dT_Label ) ) then
      associate ( I => U % Integrator )
      allocate ( I % dT_Label, source = U % dT_Label )
      end associate !-- I
    end if

  end subroutine AllocateIntegrator_F_C


  subroutine InitializePositionSpace &
               ( U, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    ! if ( .not. U % Dimensionless ) then
    !   U % Units % Time &
    !     =  UNIT % SECOND
    !   U % Units % Length &
    !     =  UNIT % KILOMETER
    !   U % Units % Coordinate_PS  &
    !     =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]
    ! end if

    call U % InitializeAtlas &
           ( RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadiusExcisionOption = RadiusExcisionOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )

    ! select type ( PS => U % Integrator % PositionSpace )
    ! class is ( Atlas_SC_Form )

    ! allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    ! select type ( GA => PS % Geometry_ASC )
    ! class is ( Geometry_ASC_Form )

    ! call U % InitializeGeometry &
    !        ( GA, PS, GeometryType, &
    !          UsePinnedMemoryOption = GeometryUseDeviceOption, &
    !          CentralMassOption = CentralMassOption )

    ! call PS % SetGeometry ( GA )
    
    ! if ( present ( GeometryUseDeviceOption ) ) then
    !   if ( GeometryUseDeviceOption ) &
    !     call GA % AllocateDevice ( )
    ! end if

    ! U % UseCoarsening = .true.
    ! call PROGRAM_HEADER % GetParameter ( U % UseCoarsening, 'UseCoarsening' )
    ! if ( U % UseCoarsening ) &
    !   call PS % SetCoarsening ( )

    ! end select !-- GA
    ! end select !-- PS

  end subroutine InitializePositionSpace


  subroutine InitializeGravitation ( U, GravitationType )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ) :: &
      GravitationType

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaAverage

    associate ( I  =>  U % Integrator )

    select case ( trim ( GravitationType ) )
    case ( 'GALILEO' )

      allocate ( Gravitation_G_Form  ::  I % Geometry_X )
      select type ( G  =>  I % Geometry_X )
        class is ( Gravitation_G_Form )
      call G % Initialize &
             ( I % X, &
               DeviceMemoryOption = U % DeviceMemory, &
               PinnedMemoryOption = U % PinnedMemory, &
               DevicesCommunicateOption = U % DevicesCommunicate )

      allocate ( iaAverage ( 0 ) )

      allocate ( U % SA_Gravitation )
      associate &
        ( SA     =>  U % SA_Gravitation, &
           A_SA  =>  U % PositionSpace_SA )
      allocate ( Gravitation_G_Form :: SA % FieldSet_SA )
      select type ( G_SA  =>  SA % FieldSet_SA )
        type is ( Gravitation_G_Form )
      call G_SA % Initialize ( A_SA, NameOption = trim ( G % Name ) // '_SA' )
      call SA % Initialize ( G, G, A_SA, iaAverageOption = iaAverage )
      end select !-- G_SA
      end associate !-- SA, etc.

      end select !-- G

    case ( 'NEWTON_SG' )

      allocate ( Gravitation_N_SG_Form  ::  I % Geometry_X )
      select type ( G  =>  I % Geometry_X )
        class is ( Gravitation_N_SG_Form )
      call G % Initialize &
             ( I % X, &
               DeviceMemoryOption = U % DeviceMemory, &
               PinnedMemoryOption = U % PinnedMemory, &
               DevicesCommunicateOption = U % DevicesCommunicate )

      allocate &
        ( iaAverage, source = [ G % POTENTIAL, G % POTENTIAL_GRADIENT_D ] )

      allocate ( U % SA_Gravitation )
      associate &
        ( SA     =>  U % SA_Gravitation, &
           A_SA  =>  U % PositionSpace_SA )
      allocate ( Gravitation_N_H_Form :: SA % FieldSet_SA )
      select type ( G_SA  =>  SA % FieldSet_SA )
        type is ( Gravitation_N_H_Form )
      call G_SA % Initialize ( A_SA, NameOption = trim ( G % Name ) // '_SA' )
      call SA % Initialize ( G, G, A_SA, iaAverageOption = iaAverage )
      end select !-- G_SA
      end associate !-- SA, etc.

      end select !-- G

    case default
      call Show ( 'GravitationType not recognized', CONSOLE % ERROR )
      call Show ( GravitationType, 'GravitationType', CONSOLE % ERROR )
      call Show ( 'Universe_F_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeGravitation', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GravitationType

    end associate !-- I

  end subroutine InitializeGravitation


  subroutine InitializeFluid ( U, FluidType )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( G  =>  I % Geometry_X )

    ! if ( .not. FC % Dimensionless ) then

    !   FC % Units % BaryonMass     =  UNIT % ATOMIC_MASS_UNIT
    !   FC % Units % NumberDensity  =  UNIT % NUMBER_DENSITY_NUCLEAR
    !   FC % Units % MassDensity    =  UNIT % MASS_DENSITY_CGS
    !   FC % Units % EnergyDensity  =  UNIT % ENERGY_DENSITY_NUCLEAR
    !   FC % Units % Temperature    =  UNIT % MEGA_ELECTRON_VOLT

    !   FC % Units % Velocity_U  &
    !     =  FC % Units % Coordinate_PS  /  FC % Units % Time

    !   FC % Units % MomentumDensity_D  &
    !     =  FC % Units % BaryonMass  *  FC % Units % NumberDensity  &
    !        *  FC % Units % Velocity_U
    !   FC % Units % MomentumDensity_D ( 2 )  &
    !     =  FC % Units % MomentumDensity_D ( 2 )  &
    !        *  FC % Units % Coordinate_PS ( 1 ) ** 2
    !   FC % Units % MomentumDensity_D ( 3 )  &
    !     =  FC % Units % MomentumDensity_D ( 3 )  &
    !        *  FC % Units % Coordinate_PS ( 1 ) ** 2

    !   FC % Units % MomentumDensity_U  &
    !     =  FC % Units % BaryonMass  *  FC % Units % NumberDensity  &
    !        *  FC % Units % Velocity_U
    !   FC % Units % MomentumDensity_U ( 2 )  &
    !     =  FC % Units % MomentumDensity_U ( 2 )  &
    !        /  FC % Units % Coordinate_PS ( 1 ) ** 2
    !   FC % Units % MomentumDensity_U ( 3 )  &
    !     =  FC % Units % MomentumDensity_U ( 3 )  &
    !        /  FC % Units % Coordinate_PS ( 1 ) ** 2

    !   FC % Units % Number           =  UNIT % SOLAR_BARYON_NUMBER
    !   FC % Units % Energy           =  UNIT % ENERGY_SOLAR_MASS
    !   FC % Units % Momentum         =  UNIT % MOMENTUM_SOLAR_MASS
    !   FC % Units % AngularMomentum  =  UNIT % SOLAR_KERR_PARAMETER

    ! end if

    select case ( trim ( FluidType ) )
    case ( 'DUST' )

      allocate ( Fluid_D_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_D_Form )
      call F % Initialize ( G, U % Units_F )

      allocate ( U % SA_Fluid )
      associate &
        ( SA     =>  U % SA_Fluid, &
           A_SA  =>  U % PositionSpace_SA )
      allocate ( Fluid_D_Form :: SA % FieldSet_SA )
      select type ( F_SA  =>  SA % FieldSet_SA )
        type is ( Fluid_D_Form )
      select type ( G_SA  =>  U % SA_Gravitation % FieldSet_SA )
        class is ( Geometry_F_Form )
      call F_SA % Initialize &
             ( G_SA, U % Units_F, NameOption = trim ( F % Name ) // '_SA' )
      call SA % Initialize ( G, F, A_SA, iaAverageOption = F % iaBalanced )
      end select !-- G_SA
      end select !-- F_SA
      end associate !-- SA, etc.

      end select !-- F
      
    case default
      call Show ( 'FluidType not recognized', CONSOLE % ERROR )
      call Show ( FluidType, 'FluidType', CONSOLE % ERROR )
      call Show ( 'Universe_F_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeFluid', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FluidType

    end associate !-- G
    end select !-- I

  end subroutine InitializeFluid


  subroutine SetBoundaryConditions ( U )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    call Show ( 'SetBoundaryConditions should be overridden', &
                CONSOLE % WARNING )
    call Show ( 'Universe_F_C__Form', 'module', CONSOLE % WARNING )
    call Show ( 'SetBoundaryConditions', 'subroutine', CONSOLE % WARNING )

  end subroutine SetBoundaryConditions


  subroutine InitializeStep ( U )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )

    allocate ( Step_RK_CS_Form :: I % Step_X )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_Form )

    select type ( G  =>  I % Geometry_X )
    class is ( Gravitation_N_SG_Form )
      S % SetSlope  =>  SetSlope_N_SG
    end select !-- G

    call S % Initialize ( F )

    end select !-- S

    end associate !-- F
    end select !-- I

  end subroutine InitializeStep


  subroutine InitializeAtlas &
               ( U, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

      class ( Universe_F_C_Form ), intent ( inout ) :: &
        U
      real ( KDR ), intent ( in ), optional :: &
        RadiusMaxOption, &
        RadiusCoreOption, &
        RadiusExcisionOption, &
        RadialRatioOption
      integer ( KDI ), intent ( in ), optional :: &
        nCellsPolarOption

      call Show ( 'InitializeAtlas should be overridden', CONSOLE % WARNING )
      call Show ( 'Universe_F_C__Form', 'module', CONSOLE % WARNING )
      call Show ( 'InitializeAtlas', 'subroutine', CONSOLE % WARNING )

  end subroutine InitializeAtlas


  subroutine Analyze_C ( I, T_A )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ) :: &
      T_A

    call I % Analyze_H ( T_A )

    !-- Spherical average


    select type ( U  =>  I % System )
      class is ( Universe_F_C_Form )
    select type ( F_SA  =>  U % SA_Fluid % FieldSet_SA )
      class is ( Fluid_D_Form )

    call F_SA % ComputeFromInitial ( )  !-- Ensure BARYON_MASS set

    call U % SA_Gravitation % Compute ( )
    call U % SA_Fluid % Compute ( )

    call F_SA % ComputeFromBalanced ( )

    end select !-- F_SA
    end select !-- U

  end subroutine Analyze_C


  subroutine Write_C ( I, T_W )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ) :: &
      T_W

    type ( TimerForm ), pointer :: &
      T_PS

    !-- Reproduce and expand Write_H functionality rather than call Write_H 
    !   functionality to avoid multiple GIS % Open calls.

    select type ( U  =>  I % System )
      class is ( Universe_F_C_Form )
    associate &
      ( GIS     =>  I % GridImageStream, &
          S_PS  =>  I % Checkpoint_X, &
          S_SA  =>  U % Stream_SA )
    T_PS  =>  S_PS % TimerWrite ( LevelOption = T_W % Level + 1 )
    call T_PS % Start ( )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call S_PS % Write &
           ( TimeOption  =  I % T  /  I % Unit_T, &
             CycleNumberOption  =  I % iCycle )
    call S_SA % Write &
           ( TimeOption  =  I % T  /  I % Unit_T, &
             CycleNumberOption  =  I % iCycle )
    call GIS % Close ( )

    call T_PS % Stop ( )
    end associate !-- GIS, etc.
    end select !-- U

  end subroutine Write_C


  subroutine SetSlope_N_SG ( S, K, iS_Option )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    integer ( KDI ) :: &
      iEnergy_B
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_B
    real ( KDR ) :: &
      Constant_G
    character ( 1 ) :: &
      StageNumber

    allocate ( Slope_DFV_N_Form :: K )
    select type ( K )
      class is ( Slope_DFV_N_Form )
    select type ( S )
      class is ( Step_RK_CS_Form )
    select type ( F  =>  S % CurrentSet )
      class is ( Fluid_D_Form ) 

    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_B ( 1 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_B ( 2 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_B ( 3 ) )

    !-- FIXME: more general fluid
    iEnergy_B  =  0
    
    !-- FIXME: nontrivial units
    Constant_G  =  1.0_KDR

    if ( present ( iS_Option ) ) then
      write ( StageNumber, fmt = '(i1.1)' ) iS_Option
      call K % Initialize &
             ( S % RiemannSolver, &
               Constant_G  =  Constant_G, &
               iVelocity_F = F % VELOCITY_U, &
               iMomentum_B = iMomentum_B, &
               iBaryonMass_F = F % BARYON_MASS, &
               iBaryonDensity_F = F % BARYON_DENSITY_B, &
               iEnergy_B = iEnergy_B, &
               SuffixOption = StageNumber )
    else
      call K % Initialize &
             ( S % RiemannSolver, &
               Constant_G  =  Constant_G, &
               iVelocity_F = F % VELOCITY_U, &
               iMomentum_B = iMomentum_B, &
               iBaryonMass_F = F % BARYON_MASS, &
               iBaryonDensity_F = F % BARYON_DENSITY_B, &
               iEnergy_B = iEnergy_B )
    end if

    end select !-- F
    end select !-- S
    end select !-- K

  end subroutine SetSlope_N_SG


end module Universe_F_C__Form
