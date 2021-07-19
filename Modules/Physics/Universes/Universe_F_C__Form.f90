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
    type ( Units_F_Form ), dimension ( : ), allocatable :: &
      Units_F
  contains
    procedure, private, pass :: &
      Initialize_F_C
    generic, public :: &
      Initialize => Initialize_F_C
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_F_C
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_F_C
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, private, pass :: &
      InitializeAtlas
  end type Universe_F_C_Form


contains


  subroutine Initialize_F_C &
               ( U, FluidType, GravitationType, NameOption, RadiusMaxOption, &
                 RadiusCoreOption, RadiusExcisionOption, RadialRatioOption, &
                 nCellsPolarOption )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType
    character ( * ), intent ( in ), optional :: &
      NameOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    character ( LDL ) :: &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a FluidCentral'

    Name  =  'FluidCentral'
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

  end subroutine Initialize_F_C


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_F ) ) &
      deallocate ( U % Units_F )

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


end module Universe_F_C__Form
