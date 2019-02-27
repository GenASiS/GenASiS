module InteractionsExamples_CSL__Form

  !-- InteractionsExamples_ChartSingleLevel_Form

  use GenASiS
  use Interactions_C__Form
  use Interactions_MWV_1__Form

  implicit none
  private

  type, public, extends ( Interactions_CSL_Template ) :: &
    InteractionsExamples_CSL_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AllocateField
  end type InteractionsExamples_CSL_Form

contains


  subroutine Initialize &
               ( IC, C, NameShort, InteractionsType, MomentsType, LengthUnit, &
                 EnergyDensityUnit, TemperatureUnit, nValues, &
                 IgnorabilityOption  )

    class ( InteractionsExamples_CSL_Form ), intent ( inout ) :: &
      IC
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      InteractionsType, &
      MomentsType
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( IC % Type == '' ) &
      IC % Type = 'an InteractionsExamples_CSL'

    call IC % IntializeTemplate_I_CSL &
           ( C, NameShort, InteractionsType, MomentsType, LengthUnit, &
             EnergyDensityUnit, TemperatureUnit, nValues, &
             IgnorabilityOption )

  end subroutine Initialize


  subroutine AllocateField ( IC )

    class ( InteractionsExamples_CSL_Form ), intent ( inout ) :: &
      IC

    select case ( trim ( IC % InteractionsType ) )
    case ( 'CONSTANT' )
      allocate ( Interactions_C_Form :: IC % Field )
    case ( 'MARSHAK_WAVE_VAYTET_1' )
      allocate ( Interactions_MWV_1_Form :: IC % Field )
    case default
      call Show ( 'InteractionsType not recognized', CONSOLE % ERROR )
      call Show ( IC % InteractionsType, 'InteractionsType', &
                  CONSOLE % ERROR )
      call Show ( 'InteractionsExamples_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine AllocateField


end module InteractionsExamples_CSL__Form
