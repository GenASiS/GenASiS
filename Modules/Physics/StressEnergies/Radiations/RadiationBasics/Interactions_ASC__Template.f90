module Interactions_ASC__Template

  !-- Interactions_AtlasSingleChart_Template

  use Basics
  use Mathematics
  use StressEnergyBasics
  use Fluids
  use Interactions_Template
  use Interactions_CSL__Template

  implicit none
  private

  type, public, extends ( Field_ASC_Template ), abstract :: &
    Interactions_ASC_Template
      class ( StressEnergyUnitsForm ), pointer :: &
        Units => null ( )
      character ( LDL ) :: &
        InteractionsType = '', &
        MomentsType = ''
  contains
    procedure ( I ), public, pass, deferred :: &
      Initialize
    procedure, public, pass :: &
      InitializeTemplate_I_ASC
    procedure, private, pass :: &
      Interactions_CSL
    generic, public :: &
      Interactions => Interactions_CSL
    procedure, public, pass :: &
      ComputeTimeScale
    procedure, public, pass :: &
      FinalizeTemplate_I_ASC
    procedure, public, pass :: &
      SetField
    procedure ( AF ), public, pass, deferred :: &
      AllocateField
  end type Interactions_ASC_Template

  abstract interface

    subroutine I ( IA, A, InteractionsType, MomentsType, Units, &
                   NameShortOption, IgnorabilityOption )
      use Basics
      use Mathematics
      use StressEnergyBasics
      import Interactions_ASC_Template
      class ( Interactions_ASC_Template ), intent ( inout ) :: &
        IA
      class ( Atlas_SC_Template ), intent ( in ) :: &
        A
      character ( * ), intent ( in ) :: &
        InteractionsType, &
        MomentsType
      class ( StressEnergyUnitsForm ), intent ( in ) :: &
        Units
      character ( * ), intent ( in ), optional :: &
        NameShortOption
      integer ( KDI ), intent ( in ), optional :: &
        IgnorabilityOption
    end subroutine I

    subroutine AF ( IA )
      import Interactions_ASC_Template
      class ( Interactions_ASC_Template ), intent ( inout ) :: &
        IA
    end subroutine AF

  end interface

    private :: &
      ComputeTimeScaleKernel_CSL_G

contains


  subroutine InitializeTemplate_I_ASC &
               ( IA, A, InteractionsType, MomentsType, Units, NameShortOption, &
                 IgnorabilityOption )

    class ( Interactions_ASC_Template ), intent ( inout ) :: &
      IA
    class ( Atlas_SC_Template ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      InteractionsType, &
      MomentsType
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      NameShort

    if ( IA % Type == '' ) &
      IA % Type = 'an Interactions_ASC'
    IA % InteractionsType = InteractionsType    
    IA % MomentsType      = MomentsType

    IA % Units => Units

    NameShort = 'Interactions'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call IA % InitializeTemplate_ASC ( A, NameShort, IgnorabilityOption )

    call Show ( IA % InteractionsType, 'InteractionsType', IA % IGNORABILITY )
    call Show ( IA % MomentsType, 'MomentsType', IA % IGNORABILITY )

  end subroutine InitializeTemplate_I_ASC


  function Interactions_CSL ( IA ) result ( I )

    class ( Interactions_ASC_Template ), intent ( in ) :: &
      IA
    class ( InteractionsTemplate ), pointer :: &
      I

    select type ( IC => IA % Chart )
    class is ( Interactions_CSL_Template )
      I => IC % Interactions ( )
    class default
      call Show ( 'Interactions Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Interactions_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'Interactions_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- IC

  end function Interactions_CSL


  subroutine ComputeTimeScale ( IA, RA, TimeScale )

    class ( Interactions_ASC_Template ), intent ( in ) :: &
      IA
    class ( Current_ASC_Template ), intent ( in ) :: &
      RA
    real ( KDR ), intent ( out ) :: &
      TimeScale

    class ( CurrentTemplate ), pointer :: &
      R
    class ( InteractionsTemplate ), pointer :: &
      I

    select type ( A => IA % Atlas )
    class is ( Atlas_SC_Form )
    select type ( C => A % Chart )
    class is ( Chart_SL_Template )

    R => RA % Current ( )
    I => IA % Interactions ( )
    call I % ComputeTimeScale ( R )

    associate ( F => I % Fluid )
    select type ( SF => F % Sources )
    class is ( Sources_F_Form )

    call ComputeTimeScaleKernel_CSL_G &
           ( C % IsProperCell, SF % Value ( :, SF % RADIATION_TIME ), &
             C % nDimensions, TimeScale )

    end select !-- SF
    end associate !-- F
    end select !-- C
    end select !-- A
    nullify ( R, I )

  end subroutine ComputeTimeScale


  impure elemental subroutine FinalizeTemplate_I_ASC ( IA )

    class ( Interactions_ASC_Template ), intent ( inout ) :: &
      IA

    nullify ( IA % Units )

    call IA % FinalizeTemplate_ASC ( )

  end subroutine FinalizeTemplate_I_ASC


  subroutine SetField ( FA )

    class ( Interactions_ASC_Template ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    call FA % AllocateField ( )

    select type ( FC => FA % Chart )
    class is ( Interactions_CSL_Template )
      call FC % Initialize &
             ( C, FA % NameShort, FA % InteractionsType, FA % MomentsType, &
               FA % Units, nValues, IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


  subroutine ComputeTimeScaleKernel_CSL_G &
               ( IsProperCell, TS, nDimensions, TimeScale )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      TS
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    real ( KDR ), intent ( out ) :: &
      TimeScale

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( TS )

    select case ( nDimensions )
    case ( 1, 2 )
      TimeScale = minval ( TS, mask = IsProperCell )
    case ( 3 )
      TimeScale = huge ( 1.0_KDR )
      !$OMP parallel do private ( iV ) &
      !$OMP reduction ( min : TimeScale )
      do iV = 1, nV
        if ( IsProperCell ( iV ) ) &
          TimeScale = min ( TimeScale, TS ( iV ) )
      end do
      !$OMP end parallel do
    end select !-- nDimensions

  end subroutine ComputeTimeScaleKernel_CSL_G


end module Interactions_ASC__Template
