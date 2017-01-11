module Fluid_CSL__Form

  !-- Fluid_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use Fluid_D__Form
  use Fluid_P_P__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: Fluid_CSL_Form
    type ( MeasuredValueForm ) :: &
      MassDensityUnit, &
      EnergyDensityUnit, &
      NumberDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      VelocityUnit
    character ( LDF ) :: &
      FluidType = ''
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Fluid_D
    procedure, public, pass :: &
      Fluid_P_P
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Fluid_CSL_Form

contains


  subroutine Initialize &
               ( FC, C, FluidType, VelocityUnit, MassDensityUnit, &
                 EnergyDensityUnit, NumberDensityUnit, TemperatureUnit, &
                 nValues, NameOutputOption )

    class ( Fluid_CSL_Form ), intent ( inout ) :: &
      FC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    character ( * ), intent ( in ) :: &
      FluidType
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit, &
      EnergyDensityUnit, &
      NumberDensityUnit, &
      TemperatureUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    if ( FC % Type == '' ) &
      FC % Type = 'a Fluid_CSL'
    FC % FluidType = FluidType

    FC % MassDensityUnit   = MassDensityUnit
    FC % EnergyDensityUnit = EnergyDensityUnit
    FC % NumberDensityUnit = NumberDensityUnit
    FC % TemperatureUnit   = TemperatureUnit
    FC % VelocityUnit      = VelocityUnit

    call FC % InitializeTemplate_CSL &
           ( C, nValues, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  function Fluid_D ( FC ) result ( F )

    class ( Fluid_CSL_Form ), intent ( in ), target :: &
      FC
    class ( Fluid_D_Form ), pointer :: &
      F
      
    class ( VariableGroupForm ), pointer :: &
      Field

    F => null ( )

    Field => FC % Field
    select type ( Field )
    class is ( Fluid_D_Form )
    F => Field
    end select !-- Field

  end function Fluid_D


  function Fluid_P_P ( FC ) result ( F )

    class ( Fluid_CSL_Form ), intent ( in ), target :: &
      FC
    class ( Fluid_P_P_Form ), pointer :: &
      F
      
    class ( VariableGroupForm ), pointer :: &
      Field

    F => null ( )

    Field => FC % Field
    select type ( Field )
    class is ( Fluid_P_P_Form )
    F => Field
    end select !-- Field

  end function Fluid_P_P


  impure elemental subroutine Finalize ( FC )

    type ( Fluid_CSL_Form ), intent ( inout ) :: &
      FC

    call FC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC, NameOption )

    class ( Fluid_CSL_Form ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in ), optional :: &
      NameOption

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % FluidType ) )
    case ( 'DUST' )
      allocate ( Fluid_D_Form :: FC % Field )
      select type ( F => FC % Field )
      type is ( Fluid_D_Form )
        call F % Initialize &
               ( FC % VelocityUnit, FC % MassDensityUnit, FC % nValues, &
                 NameOption = NameOption )
        call F % SetOutput ( FC % FieldOutput )
      end select !-- F
    case ( 'POLYTROPIC' )
      allocate ( Fluid_P_P_Form :: FC % Field )
      select type ( F => FC % Field )
      type is ( Fluid_P_P_Form )
        call F % Initialize &
               ( FC % VelocityUnit, FC % MassDensityUnit, &
                 FC % EnergyDensityUnit, FC % TemperatureUnit, FC % nValues, &
                 NameOption = NameOption )
        call F % SetOutput ( FC % FieldOutput )
      end select !-- F
    case default
      call Show ( 'FluidType not recognized', CONSOLE % ERROR )
      call Show ( FC % FluidType, 'FluidType', CONSOLE % ERROR )
      call Show ( 'Fluid_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FluidType

  end subroutine SetField


end module Fluid_CSL__Form
