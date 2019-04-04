module InteractionsExamples_CSL__Form

  !-- InteractionsExamples_ChartSingleLevel_Form

  use GenASiS
  use Interactions_C__Form
  use Interactions_MWV_1__Form
  use Interactions_MWV_2__Form
  use Interactions_MWV_3__Form

  implicit none
  private

  type, public, extends ( Interactions_CSL_Template ) :: &
    InteractionsExamples_CSL_Form
  contains
    final :: &
      Finalize
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      AllocateField
  end type InteractionsExamples_CSL_Form

contains

  
  impure elemental subroutine Finalize ( IC )

    type ( InteractionsExamples_CSL_Form ), intent ( inout ) :: &
      IC

    call IC % FinalizeTemplate_I_CSL ( )

  end subroutine Finalize


  subroutine SetType ( IC )

    class ( InteractionsExamples_CSL_Form ), intent ( inout ) :: &
      IC

    IC % Type = 'an InteractionsExamples_CSL'

  end subroutine SetType


  subroutine AllocateField ( IC )

    class ( InteractionsExamples_CSL_Form ), intent ( inout ) :: &
      IC

    select case ( trim ( IC % InteractionsType ) )
    case ( 'CONSTANT' )
      allocate ( Interactions_C_Form :: IC % Field )
    case ( 'MARSHAK_WAVE_VAYTET_1' )
      allocate ( Interactions_MWV_1_Form :: IC % Field )
    case ( 'MARSHAK_WAVE_VAYTET_2' )
      allocate ( Interactions_MWV_2_Form :: IC % Field )
    case ( 'MARSHAK_WAVE_VAYTET_3' )
      allocate ( Interactions_MWV_3_Form :: IC % Field )
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
