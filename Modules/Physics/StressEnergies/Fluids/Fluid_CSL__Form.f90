module Fluid_CSL__Form

  !-- Fluid_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use FluidFeatures_Template
  use Fluid_D__Form
  use Fluid_P_I__Form
  use Fluid_P_HN__Form
  use Sources_F__Form
  use Sources_F_CSL__Form
  use FluidFeatures_CSL__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: Fluid_CSL_Form
    real ( KDR ) :: &
      LimiterParameter, &
      BaryonMassReference
    type ( MeasuredValueForm ) :: &
      BaryonMassUnit, &
      NumberDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U_Unit, &
      MomentumDensity_D_Unit
    logical ( KDL ) :: &
      UseLimiter
    character ( LDF ) :: &
      FluidType = '', &
      RiemannSolverType = ''
    class ( Sources_F_CSL_Form ), pointer :: &
      Sources_CSL => null ( )
    class ( FluidFeatures_CSL_Form ), pointer :: &
      Features_CSL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Fluid_D
    procedure, public, pass :: &
      Fluid_P_I
    procedure, public, pass :: &
      Fluid_P_HN
    procedure, public, pass :: &
      SetSources
    procedure, public, pass :: &
      SetFeatures
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Fluid_CSL_Form

contains


  subroutine Initialize &
               ( FC, C, NameShort, FluidType, RiemannSolverType, UseLimiter, &
                 Velocity_U_Unit, MomentumDensity_D_Unit, BaryonMassUnit, &
                 NumberDensityUnit, EnergyDensityUnit, TemperatureUnit, &
                 BaryonMassReference, LimiterParameter, nValues, &
                 IgnorabilityOption )

    class ( Fluid_CSL_Form ), intent ( inout ) :: &
      FC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      FluidType, &
      RiemannSolverType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit, &
      MomentumDensity_D_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
      BaryonMassUnit, &
      NumberDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    real ( KDR ), intent ( in ) :: &
      BaryonMassReference, &
      LimiterParameter
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FC % Type == '' ) &
      FC % Type = 'a Fluid_CSL'
    FC % FluidType           = FluidType
    FC % RiemannSolverType   = RiemannSolverType
    FC % UseLimiter          = UseLimiter
    FC % LimiterParameter    = LimiterParameter
    FC % BaryonMassReference = BaryonMassReference

    FC % NumberDensityUnit      = NumberDensityUnit
    FC % EnergyDensityUnit      = EnergyDensityUnit
    FC % TemperatureUnit        = TemperatureUnit
    FC % Velocity_U_Unit        = Velocity_U_Unit
    FC % MomentumDensity_D_Unit = MomentumDensity_D_Unit

    call FC % InitializeTemplate_CSL &
           ( C, NameShort, nValues, IgnorabilityOption )

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


  function Fluid_P_I ( FC ) result ( F )

    class ( Fluid_CSL_Form ), intent ( in ), target :: &
      FC
    class ( Fluid_P_I_Form ), pointer :: &
      F
      
    class ( VariableGroupForm ), pointer :: &
      Field

    F => null ( )

    Field => FC % Field
    select type ( Field )
    class is ( Fluid_P_I_Form )
    F => Field
    end select !-- Field

  end function Fluid_P_I


  function Fluid_P_HN ( FC ) result ( F )

    class ( Fluid_CSL_Form ), intent ( in ), target :: &
      FC
    class ( Fluid_P_HN_Form ), pointer :: &
      F
      
    class ( VariableGroupForm ), pointer :: &
      Field

    F => null ( )

    Field => FC % Field
    select type ( Field )
    class is ( Fluid_P_HN_Form )
    F => Field
    end select !-- Field

  end function Fluid_P_HN


  subroutine SetSources ( FC, SFC )

    class ( Fluid_CSL_Form ), intent ( inout ) :: &
      FC
    class ( Sources_F_CSL_Form ), intent ( in ), target :: &
      SFC

    class ( Fluid_D_Form ), pointer :: &
      F

    FC % Sources_CSL => SFC

    F => FC % Fluid_D ( )

    call Show ( 'Setting Sources', F % IGNORABILITY - 1 )
    call Show ( FC % Name, 'Name', F % IGNORABILITY - 1 )
    call Show ( SFC % Name, 'Sources', F % IGNORABILITY - 1 )

    select type ( SF => SFC % Field )
    class is ( Sources_F_Form )
      call F % SetSources ( SF )
    end select !-- SF

    nullify ( F )

  end subroutine SetSources


  subroutine SetFeatures ( FC, FFC )

    class ( Fluid_CSL_Form ), intent ( inout ) :: &
      FC
    class ( FluidFeatures_CSL_Form ), intent ( in ), target :: &
      FFC

    class ( Fluid_D_Form ), pointer :: &
      F

    FC % Features_CSL => FFC

    F => FC % Fluid_D ( )

    call Show ( 'Setting Features', F % IGNORABILITY - 1 )
    call Show ( FC % Name, 'Name', F % IGNORABILITY - 1 )
    call Show ( FFC % Name, 'Sources', F % IGNORABILITY - 1 )

    select type ( FF => FFC % Field )
    class is ( FluidFeaturesTemplate )
    call F % SetFeatures ( FF )
    end select !-- FF

    nullify ( F )

  end subroutine SetFeatures


  impure elemental subroutine Finalize ( FC )

    type ( Fluid_CSL_Form ), intent ( inout ) :: &
      FC

    call FC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( Fluid_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % FluidType ) )
    case ( 'DUST' )
      allocate ( Fluid_D_Form :: FC % Field )
      select type ( F => FC % Field )
      type is ( Fluid_D_Form )
        call F % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_D_Unit, &
                 FC % BaryonMassUnit, FC % NumberDensityUnit, &
                 FC % BaryonMassReference, FC % LimiterParameter, &
                 FC % nValues, NameOption = FC % NameShort )
        call F % SetPrimitiveConserved ( )
        call F % SetOutput ( FC % FieldOutput )
      end select !-- F
    case ( 'IDEAL' )
      allocate ( Fluid_P_I_Form :: FC % Field )
      select type ( F => FC % Field )
      type is ( Fluid_P_I_Form )
        call F % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_D_Unit, &
                 FC % BaryonMassUnit, FC % NumberDensityUnit, &
                 FC % EnergyDensityUnit, FC % TemperatureUnit, &
                 FC % BaryonMassReference, FC % LimiterParameter, &
                 FC % nValues, NameOption = FC % NameShort )
        call F % SetPrimitiveConserved ( )
        call F % SetOutput ( FC % FieldOutput )
      end select !-- F
    case ( 'HEAVY_NUCLEUS' )
      allocate ( Fluid_P_HN_Form :: FC % Field )
      select type ( F => FC % Field )
      type is ( Fluid_P_HN_Form )
        call F % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_D_Unit, &
                 FC % BaryonMassUnit, FC % NumberDensityUnit, &
                 FC % EnergyDensityUnit, FC % TemperatureUnit, &
                 FC % BaryonMassReference, FC % LimiterParameter, &
                 FC % nValues, NameOption = FC % NameShort )
        call F % SetPrimitiveConserved ( )
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
