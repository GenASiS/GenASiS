module FluidFeatures_CSL__Form

  !-- FluidFeatures_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use FluidFeatures_P__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: FluidFeatures_CSL_Form
    real ( KDR ) :: &
      ShockThreshold
    character ( LDF ) :: &
      FluidType = ''
    class ( Field_CSL_Template ), pointer :: &
      Fluid_CSL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type FluidFeatures_CSL_Form

contains


  subroutine Initialize &
               ( FFC, Fluid_CSL, NameShort, FluidType, ShockThreshold, &
                 nValues, IgnorabilityOption )

    class ( FluidFeatures_CSL_Form ), intent ( inout ) :: &
      FFC
    class ( Field_CSL_Template ), intent ( in ), target :: &
      Fluid_CSL
    character ( * ), intent ( in ) :: &
      NameShort, &
      FluidType
    real ( KDR ), intent ( in ) :: &
      ShockThreshold
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FFC % Type == '' ) &
      FFC % Type = 'a FluidFeatures_CSL'
    FFC % FluidType = FluidType

    FFC % ShockThreshold = ShockThreshold
    FFC % Fluid_CSL => Fluid_CSL

    call FFC % InitializeTemplate_CSL &
           ( Fluid_CSL % Chart, NameShort, nValues, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( FFC )

    type ( FluidFeatures_CSL_Form ), intent ( inout ) :: &
      FFC

    nullify ( FFC % Fluid_CSL )

    call FFC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( FluidFeatures_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % FluidType ) )
    case ( 'POLYTROPIC', 'NON_RELATIVISTIC', 'MEAN_HEAVY_NUCLEUS' )
      allocate ( FluidFeatures_P_Form :: FC % Field )
      select type ( FF => FC % Field )
      type is ( FluidFeatures_P_Form )
        call FF % Initialize &
               ( FC % Fluid_CSL % Field, FC % Chart, FC % ShockThreshold, &
                 FC % nValues, NameOption = FC % NameShort )
        call FF % SetOutput ( FC % FieldOutput )
      end select !-- F
    case default
      call Show ( 'FluidType not recognized', CONSOLE % ERROR )
      call Show ( FC % FluidType, 'FluidType', CONSOLE % ERROR )
      call Show ( 'FluidFeatures_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FluidType

  end subroutine SetField


end module FluidFeatures_CSL__Form
