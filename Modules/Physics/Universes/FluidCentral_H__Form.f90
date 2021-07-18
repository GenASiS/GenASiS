module FluidCentral_H__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_H__Form

  implicit none
  private

  type, public, extends ( Universe_H_Form ) :: FluidCentral_H_Form
    type ( Units_F_Form ), dimension ( : ), allocatable :: &
      Units_F
  contains
    procedure, private, pass :: &
      Initialize_FC
    generic, public :: &
      Initialize => Initialize_FC
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_FC
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_FC
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, private, pass :: &
      InitializeAtlas
  end type FluidCentral_H_Form


contains


  subroutine Initialize_FC &
               ( FC, FluidType, GravitationType, NameOption, RadiusMaxOption, &
                 RadiusCoreOption, RadiusExcisionOption, RadialRatioOption, &
                 nCellsPolarOption )

    class ( FluidCentral_H_Form ), intent ( inout ) :: &
      FC
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

    if ( FC % Type == '' ) &
      FC % Type = 'a FluidCentral'

    Name  =  'FluidCentral'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call FC % Universe_H_Form % Initialize ( NameOption = Name )

    call FC % AllocateIntegrator &
           ( )
    call FC % InitializePositionSpace &
           ( RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadiusExcisionOption = RadiusExcisionOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )

  end subroutine Initialize_FC


  impure elemental subroutine Finalize ( FC )

    type ( FluidCentral_H_Form ), intent ( inout ) :: &
      FC

    if ( allocated ( FC % Units_F ) ) &
      deallocate ( FC % Units_F )

  end subroutine Finalize


  subroutine AllocateIntegrator_FC ( FC )

    class ( FluidCentral_H_Form ), intent ( inout ) :: &
      FC

    allocate ( Integrator_CS_Form :: FC % Integrator )

    if ( allocated ( FC % dT_Label ) ) then
      associate ( I => FC % Integrator )
      allocate ( I % dT_Label, source = FC % dT_Label )
      end associate !-- I
    end if

  end subroutine AllocateIntegrator_FC


  subroutine InitializePositionSpace &
               ( FC, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

    class ( FluidCentral_H_Form ), intent ( inout ) :: &
      FC
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    ! if ( .not. FC % Dimensionless ) then
    !   FC % Units % Time &
    !     =  UNIT % SECOND
    !   FC % Units % Length &
    !     =  UNIT % KILOMETER
    !   FC % Units % Coordinate_PS  &
    !     =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]
    ! end if

    call FC % InitializeAtlas &
           ( RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadiusExcisionOption = RadiusExcisionOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )

    ! select type ( PS => FC % Integrator % PositionSpace )
    ! class is ( Atlas_SC_Form )

    ! allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    ! select type ( GA => PS % Geometry_ASC )
    ! class is ( Geometry_ASC_Form )

    ! call FC % InitializeGeometry &
    !        ( GA, PS, GeometryType, &
    !          UsePinnedMemoryOption = GeometryUseDeviceOption, &
    !          CentralMassOption = CentralMassOption )

    ! call PS % SetGeometry ( GA )
    
    ! if ( present ( GeometryUseDeviceOption ) ) then
    !   if ( GeometryUseDeviceOption ) &
    !     call GA % AllocateDevice ( )
    ! end if

    ! FC % UseCoarsening = .true.
    ! call PROGRAM_HEADER % GetParameter ( FC % UseCoarsening, 'UseCoarsening' )
    ! if ( FC % UseCoarsening ) &
    !   call PS % SetCoarsening ( )

    ! end select !-- GA
    ! end select !-- PS

  end subroutine InitializePositionSpace


  subroutine InitializeAtlas &
               ( FC, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

      class ( FluidCentral_H_Form ), intent ( inout ) :: &
        FC
      real ( KDR ), intent ( in ), optional :: &
        RadiusMaxOption, &
        RadiusCoreOption, &
        RadiusExcisionOption, &
        RadialRatioOption
      integer ( KDI ), intent ( in ), optional :: &
        nCellsPolarOption

      call Show ( 'InitializeAtlas should be overridden', CONSOLE % WARNING )
      call Show ( 'FluidCentral_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'InitializeAtlas', 'subroutine', CONSOLE % WARNING )

  end subroutine InitializeAtlas


end module FluidCentral_H__Form
