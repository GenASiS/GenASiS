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
  end type FluidSymmetricCurvilinearForm

contains


  subroutine Initialize_FSC &
               ( FSC, Name, FluidType, FinishTimeOption, RadiusMaxOption, &
                 nCellsRadiusOption )

    class ( FluidSymmetricCurvilinearForm ), intent ( inout ) :: &
      FSC
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption

    if ( FSC % Type == '' ) &
      FSC % Type = 'a FluidSymmetricCurvilinear'

    call FSC % InitializeTemplate ( Name )

    call FSC % AllocateIntegrator &
           ( )
    call FSC % InitializePositionSpace &
           ( RadiusMaxOption = RadiusMaxOption, &
             nCellsRadiusOption = nCellsRadiusOption )

!     !-- Fluid

!     allocate ( Fluid_ASC_Form :: FSC % Current_ASC )
!     select type ( FA => FSC % Current_ASC )  !-- FluidAtlas
!     class is ( Fluid_ASC_Form )
!     call FA % Initialize ( PS, FluidType )


!     !-- Step

!     allocate ( Step_RK2_C_ASC_Form :: FSC % Step )
!     select type ( S => FSC % Step )
!     class is ( Step_RK2_C_ASC_Form )
!     call S % Initialize ( FSC, FA, Name )
!     S % ApplySources % Pointer => ApplyCurvilinear_F
!     end select !-- S


!     !-- Template

!     call FSC % InitializeTemplate_C_PS &
!            ( Name, FinishTimeOption = FinishTimeOption )


!     !-- Cleanup

!     end select !-- FA
!     end select !-- PS

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


end module FluidSymmetricCurvilinear_Form
