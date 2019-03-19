module FluidSymmetricCurvilinear_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use Universe_Template

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: &
    FluidSymmetricCurvilinearForm
  contains
    procedure, private, pass :: &
      Initialize_FSC
    generic, public :: &
      Initialize => Initialize_FSC
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_FSC
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_FSC
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      InitializeStep
  end type FluidSymmetricCurvilinearForm

contains


  subroutine Initialize_FSC &
               ( FSC, FluidType, Name, FinishTimeOption, CourantFactorOption, &
                 RadiusMaxOption, nCellsRadiusOption, nWriteOption )

    class ( FluidSymmetricCurvilinearForm ), intent ( inout ) :: &
      FSC
    character ( * ), intent ( in )  :: &
      FluidType, &
      Name
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      RadiusMaxOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption, &
      nWriteOption

    if ( FSC % Type == '' ) &
      FSC % Type = 'a FluidSymmetricCurvilinear'

    call FSC % InitializeTemplate ( Name )

    call FSC % AllocateIntegrator &
           ( )
    call FSC % InitializePositionSpace &
           ( RadiusMaxOption = RadiusMaxOption, &
             nCellsRadiusOption = nCellsRadiusOption )
    call FSC % InitializeFluid &
           ( FluidType )
    call FSC % InitializeStep &
           ( Name ) 

    select type ( I => FSC % Integrator )
    class is ( Integrator_C_PS_Form )
      call I % Initialize &
             ( FSC, Name, TimeUnitOption = FSC % Units % Time, &
               FinishTimeOption = FinishTimeOption, &
               CourantFactorOption = CourantFactorOption, &
               nWriteOption = nWriteOption )
    end select !-- I

  end subroutine Initialize_FSC


  impure elemental subroutine Finalize ( FSC )

    type ( FluidSymmetricCurvilinearForm ), intent ( inout ) :: &
      FSC

    call FSC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine AllocateIntegrator_FSC ( FSC )

    class ( FluidSymmetricCurvilinearForm ), intent ( inout ) :: &
      FSC

    allocate ( Integrator_C_PS_Form :: FSC % Integrator )

  end subroutine AllocateIntegrator_FSC


  subroutine InitializePositionSpace &
               ( FSC, RadiusMaxOption, nCellsRadiusOption )

    class ( FluidSymmetricCurvilinearForm ), intent ( inout ) :: &
      FSC
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption

    associate ( I => FSC % Integrator )

    allocate ( Atlas_SC_SC_Form :: I % PositionSpace )
    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )
    
    call PS % CreateChart_SC &
           ( RadiusMaxOption = RadiusMaxOption, &
             nCellsRadiusOption = nCellsRadiusOption )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    call GA % Initialize ( PS, 'GALILEAN' )
    call PS % SetGeometry ( GA )

    end select !-- GA
    end select !-- PS
    end associate !-- I

  end subroutine InitializePositionSpace


  subroutine InitializeFluid ( FSC, FluidType )

    class ( FluidSymmetricCurvilinearForm ), intent ( inout ) :: &
      FSC
    character ( * ), intent ( in )  :: &
      FluidType

    select type ( I => FSC % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Fluid_ASC_Form :: I % Current_ASC )
    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
      call FA % Initialize ( PS, FluidType, FSC % Units )
    end select !-- FA

    end select !-- PS
    end select !-- I

  end subroutine InitializeFluid


  subroutine InitializeStep ( FSC, Name )

    class ( FluidSymmetricCurvilinearForm ), intent ( inout ) :: &
      FSC
    character ( * ), intent ( in )  :: &
      Name

    select type ( I => FSC % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    allocate ( Step_RK2_C_ASC_Form :: I % Step )
    select type ( S => I % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( I, FA, Name )
    S % ApplySources % Pointer => ApplyCurvilinear_F

    end select !-- S
    end select !-- FA
    end select !-- I

  end subroutine InitializeStep


end module FluidSymmetricCurvilinear_Form
