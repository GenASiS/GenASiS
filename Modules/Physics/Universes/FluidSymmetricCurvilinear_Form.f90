module FluidSymmetricCurvilinear_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: &
    FluidSymmetricCurvilinearForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type FluidSymmetricCurvilinearForm

contains


  subroutine Initialize &
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


    !-- PositionSpace

    allocate ( Atlas_SC_SC_Form :: FSC % PositionSpace )
    select type ( PS => FSC % PositionSpace )
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


    !-- Fluid

    allocate ( Fluid_ASC_Form :: FSC % Current_ASC )
    select type ( FA => FSC % Current_ASC )  !-- FluidAtlas
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, FluidType )


    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: FSC % Step )
    select type ( S => FSC % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FSC, FA, Name )
    S % ApplySources % Pointer => ApplyCurvilinear_F
    end select !-- S


    !-- Template

    call FSC % InitializeTemplate_C_PS &
           ( Name, FinishTimeOption = FinishTimeOption )


    !-- Cleanup

    end select !-- FA
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( FSC )

    type ( FluidSymmetricCurvilinearForm ), intent ( inout ) :: &
      FSC

    call FSC % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


end module FluidSymmetricCurvilinear_Form
