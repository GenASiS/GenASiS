module InteractionsNeutrinos_ASC__Form

  !-- InteractionsNeutrinos_AtlasSingleChart_Form

  use RadiationBasics
  use InteractionsNeutrinos_CSL__Form

  implicit none
  private

  type, public, extends ( Interactions_ASC_Template ) :: &
    InteractionsNeutrinos_ASC_Form
  contains
    final :: &
      Finalize
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      AllocateField
  end type InteractionsNeutrinos_ASC_Form

contains


  impure elemental subroutine Finalize ( IA )

    type ( InteractionsNeutrinos_ASC_Form ), intent ( inout ) :: &
      IA

    call IA % FinalizeTemplate_I_ASC ( )

  end subroutine Finalize


  subroutine SetType ( IA )

    class ( InteractionsNeutrinos_ASC_Form ), intent ( inout ) :: &
      IA

    if ( IA % Type == '' ) &
      IA % Type = 'an InteractionsNeutrinos_ASC'

  end subroutine SetType


  subroutine AllocateField ( IA )

    class ( InteractionsNeutrinos_ASC_Form ), intent ( inout ) :: &
      IA

    allocate ( InteractionsNeutrinos_CSL_Form :: IA % Chart )

  end subroutine AllocateField


end module InteractionsNeutrinos_ASC__Form
