module Interactions_G_G_CSL__Form

  use GenASiS
  use Interactions_G_G__Form

  implicit none
  private

  type, public, extends ( Interactions_CSL_Template ) :: &
    Interactions_G_G_CSL_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AllocateField
  end type Interactions_G_G_CSL_Form

contains


  subroutine Initialize &
               ( IC, C, NameShort, InteractionsType, LengthUnit, &
                 EnergyDensityUnit, TemperatureUnit, nValues, &
                 IgnorabilityOption  )

    class ( Interactions_G_G_CSL_Form ), intent ( inout ) :: &
      IC
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      InteractionsType
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( IC % Type == '' ) &
      IC % Type = 'an Interactions_G_G_CSL'

    call IC % IntializeTemplate_I_CSL &
           ( C, NameShort, InteractionsType, LengthUnit, &
             EnergyDensityUnit, TemperatureUnit, nValues, &
             IgnorabilityOption )

  end subroutine Initialize


  subroutine AllocateField ( IC )

    class ( Interactions_G_G_CSL_Form ), intent ( inout ) :: &
      IC

    select case ( trim ( IC % InteractionsType ) )
    case ( 'GENERIC_GREY' )
      allocate ( Interactions_G_G_Form :: IC % Field )
    case default
      call Show ( 'InteractionsType not recognized', CONSOLE % ERROR )
      call Show ( IC % InteractionsType, 'InteractionsType', &
                  CONSOLE % ERROR )
      call Show ( 'Interactions_G_G_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine AllocateField


end module Interactions_G_G_CSL__Form
