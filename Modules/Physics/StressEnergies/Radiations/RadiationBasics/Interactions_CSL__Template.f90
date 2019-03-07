module Interactions_CSL__Template

  !-- Interactions_ChartSingleLevel_Template

  use Basics
  use Mathematics
  use StressEnergyBasics
  use Interactions_Template

  implicit none
  private

  type, public, extends ( Field_CSL_Template ), abstract :: &
    Interactions_CSL_Template
      class ( StressEnergyUnitsForm ), pointer :: &
        Units => null ( )
      character ( LDL ) :: &
        InteractionsType = '', &
        MomentsType = ''
  contains
    procedure ( I ), public, pass, deferred :: &
      Initialize
    procedure, public, pass :: &
      IntializeTemplate_I_CSL
    procedure, public, pass :: &
      Interactions
    procedure, public, pass :: &
      FinalizeTemplate_I_CSL
    procedure, public, pass :: &
      SetField
    procedure ( AF ), public, pass, deferred :: &
      AllocateField
  end type Interactions_CSL_Template

  abstract interface

    subroutine I ( IC, C, NameShort, InteractionsType, MomentsType, Units, &
                   nValues, IgnorabilityOption )
      use Basics
      use Mathematics
      use StressEnergyBasics
      import Interactions_CSL_Template
      class ( Interactions_CSL_Template ), intent ( inout ) :: &
        IC
      class ( ChartHeader_SL_Form ), intent ( in ) :: &
        C
      character ( * ), intent ( in ) :: &
        NameShort, &
        InteractionsType, &
        MomentsType
      class ( StressEnergyUnitsForm ), intent ( in ) :: &
        Units
      integer ( KDI ), intent ( in ) :: &
        nValues
      integer ( KDI ), intent ( in ), optional :: &
        IgnorabilityOption
    end subroutine I

    subroutine AF ( IC )
      import Interactions_CSL_Template
      class ( Interactions_CSL_Template ), intent ( inout ) :: &
        IC
    end subroutine AF

  end interface

contains


  subroutine IntializeTemplate_I_CSL &
               ( IC, C, NameShort, InteractionsType, MomentsType, Units, &
                 nValues, IgnorabilityOption )

    class ( Interactions_CSL_Template ), intent ( inout ) :: &
      IC
    class ( ChartHeader_SL_Form ), intent ( in ) :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      InteractionsType, &
      MomentsType
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( IC % Type == '' ) &
      IC % Type = 'an Interactions_CSL'
    IC % InteractionsType = InteractionsType
    IC % MomentsType      = MomentsType

    IC % Units => Units

    call IC % InitializeTemplate_CSL &
           ( C, NameShort, nValues, IgnorabilityOption )

  end subroutine IntializeTemplate_I_CSL


  function Interactions ( IC ) result ( I )

    class ( Interactions_CSL_Template ), intent ( in ), target :: &
      IC
    class ( InteractionsTemplate ), pointer :: &
      I
      
    class ( StorageForm ), pointer :: &
      Field

    I => null ( )

    Field => IC % Field
    select type ( Field )
    class is ( InteractionsTemplate )
    I => Field
    end select !-- Field

  end function Interactions


  impure elemental subroutine FinalizeTemplate_I_CSL ( IC )

    class ( Interactions_CSL_Template ), intent ( inout ) :: &
      IC

    call IC % FinalizeTemplate_CSL ( )

  end subroutine FinalizeTemplate_I_CSL


  subroutine SetField ( FC )

    class ( Interactions_CSL_Template ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    call FC % AllocateField ( )

    select type ( I => FC % Field )
    class is ( InteractionsTemplate )
      call I % Initialize &
             ( FC % MomentsType, FC % Units, FC % nValues, &
               NameOption = FC % NameShort )
      call I % SetOutput ( FC % FieldOutput )
    end select !-- I

  end subroutine SetField


end module Interactions_CSL__Template
