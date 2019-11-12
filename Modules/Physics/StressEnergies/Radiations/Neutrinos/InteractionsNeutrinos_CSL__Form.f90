module InteractionsNeutrinos_CSL__Form

  !-- InteractionsNeutrinos_ChartSingleLevel_Form

  use Basics
  use RadiationBasics
  use Interactions_B85__Form
  use Interactions_OCO__Form

  implicit none
  private

  type, public, extends ( Interactions_CSL_Template ) :: &
    InteractionsNeutrinos_CSL_Form
  contains
    final :: &
      Finalize
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      AllocateField
  end type InteractionsNeutrinos_CSL_Form

contains


  impure elemental subroutine Finalize ( IC )

    type ( InteractionsNeutrinos_CSL_Form ), intent ( inout ) :: &
      IC

    call IC % FinalizeTemplate_I_CSL ( )

  end subroutine Finalize


  subroutine SetType ( IC )

    class ( InteractionsNeutrinos_CSL_Form ), intent ( inout ) :: &
      IC

    IC % Type = 'an InteractionsNeutrinos_CSL'

  end subroutine SetType


  subroutine AllocateField ( IC )

    class ( InteractionsNeutrinos_CSL_Form ), intent ( inout ) :: &
      IC

    select case ( trim ( IC % InteractionsType ) )
    case ( 'BRUENN_1985' )
      allocate ( Interactions_B85_Form :: IC % Field )
    case ( 'O_CONNOR_OTT' )
      allocate ( Interactions_OCO_Form :: IC % Field )
    case default
      call Show ( 'InteractionsType not recognized', CONSOLE % ERROR )
      call Show ( IC % InteractionsType, 'InteractionsType', &
                  CONSOLE % ERROR )
      call Show ( 'InteractionsNeutrinos_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine AllocateField


end module InteractionsNeutrinos_CSL__Form
