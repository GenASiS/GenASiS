module Interactions_BSLL_ASC_CSLD__Form

  use Basics
  use Mathematics
  use Interactions_Template
  use Interactions_F__Form
  use Interactions_ASC__Form

  implicit none
  private

  type, public, extends ( Field_BSLL_ASC_CSLD_Template ) :: &
    Interactions_BSLL_ASC_CSLD_Form
      type ( MeasuredValueForm ) :: &
        LengthUnit, &
        EnergyDensityUnit
      character ( LDF ) :: &
        InteractionsType = ''
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Interactions
    procedure, public, pass :: &
      Interactions_F
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Interactions_BSLL_ASC_CSLD_Form

contains


  subroutine Initialize &
               ( IB, B, InteractionsType, NameShortOption, LengthUnitOption, &
                 EnergyDensityUnitOption )

    class ( Interactions_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in )  :: &
      InteractionsType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      LengthUnitOption, &
      EnergyDensityUnitOption

    character ( LDL ) :: &
      NameShort
    ! class ( GeometryFlatForm ), pointer :: &
    !   GF

    if ( IB % Type == '' ) &
      IB % Type = 'an Interactions_BSLL_ASC_CSLD'
    IB % InteractionsType = InteractionsType

    ! select type ( B )
    ! class is ( Bundle_SLL_ASC_CSLD_Form )

    ! IB % nFibers = B % nFibers

    ! IB % nBaseValues &
    !   = B % Base_CSLD % nProperCells  +  B % Base_CSLD % nGhostCells
    ! IB % ChartBase => B % Base_CSLD

    ! allocate ( IB % iaBaseCell ( size ( B % iaBaseCell ) ) )
    ! IB % iaBaseCell = B % iaBaseCell

    ! GF => B % GeometryFiber ( )
    ! associate ( Energy => GF % Value ( :, GF % CENTER ( 1 ) ) )
    ! IB % nEnergyValues = size ( Energy )
    ! allocate ( IB % Energy ( IB % nEnergyValues ) )
    ! IB % Energy = Energy
    ! end associate !-- Energy

    ! IB % EnergyUnit = GF % Unit ( GF % CENTER ( 1 ) )

    ! associate ( AF => B % FiberMaster )
    ! select type ( CF => AF % Chart )
    ! class is ( Chart_SLL_Form )
    !   IB % ChartFiber => CF
    ! end select !--C
    ! end associate !-- AF

    ! class default
    !   call Show ( 'Bundle type not recognized', CONSOLE % ERROR )
    !   call Show ( 'Interactions_BSLL_ASC_CSLD__Form', 'module', &
    !               CONSOLE % ERROR )
    !   call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
    !   call PROGRAM_HEADER % Abort ( )
    ! end select !-- B

    if ( present ( LengthUnitOption ) ) &
      IB % LengthUnit = LengthUnitOption
    if ( present ( EnergyDensityUnitOption ) ) &
      IB % EnergyDensityUnit = EnergyDensityUnitOption

    NameShort = 'Interactions'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call IB % InitializeTemplate_BSLL_ASC_CSLD ( B, NameShort )

    ! nullify ( GF )

  end subroutine Initialize


  function Interactions ( IB, iFiber ) result ( IF )

    class ( Interactions_BSLL_ASC_CSLD_Form ), intent ( in ), target :: &
      IB
    integer ( KDI ), intent ( in ) :: &
      iFiber
    class ( InteractionsTemplate ), pointer :: &
      IF

    select type ( IA => IB % Fiber % Atlas ( iFiber ) % Element )
    class is ( Field_ASC_Template )
      select type ( IC => IA % Chart )
      class is ( Field_CSL_Template )   
        select type ( I => IC % Field )
        class is ( InteractionsTemplate )
          IF => I
        end select !-- I
      end select !-- IC
    end select !-- IA

  end function Interactions


  function Interactions_F ( IB, iFiber ) result ( IF )

    class ( Interactions_BSLL_ASC_CSLD_Form ), intent ( in ), target :: &
      IB
    integer ( KDI ), intent ( in ) :: &
      iFiber
    class ( Interactions_F_Form ), pointer :: &
      IF

    select type ( IA => IB % Fiber % Atlas ( iFiber ) % Element )
    class is ( Field_ASC_Template )
      select type ( IC => IA % Chart )
      class is ( Field_CSL_Template )   
        select type ( I => IC % Field )
        class is ( Interactions_F_Form )
          IF => I
        end select !-- I
      end select !-- IC
    end select !-- IA

  end function Interactions_F


  impure elemental subroutine Finalize ( IB )

    type ( Interactions_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB

    call IB % FinalizeTemplate_BSLL_ASC_CSLD ( )

  end subroutine Finalize


  subroutine SetField ( FB )

    class ( Interactions_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      FB

    integer ( KDI ) :: &
      iF  !-- iFiber

    select type ( B => FB % Bundle )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    allocate ( FB % Fiber )
    associate ( FBF => FB % Fiber )
    call FBF % Initialize ( B % nFibers )

    do iF = 1, B % nFibers
      allocate ( Interactions_ASC_Form :: FBF % Atlas ( iF ) % Element )
      select type ( IA => FBF % Atlas ( iF ) % Element )
      class is ( Interactions_ASC_Form )

      select type ( AF => B % Fiber % Atlas ( iF ) % Element )
      class is ( Atlas_SC_Form )

      call IA % Initialize &
             ( AF, FB % InteractionsType, &
               NameShortOption = FB % NameShort, & 
               LengthUnitOption = FB % LengthUnit, &
               EnergyDensityUnitOption = FB % EnergyDensityUnit )

      end select !-- AF
      end select !-- IA
    end do !-- iF
    
    end associate !-- FBF
    end select !-- B

  end subroutine SetField


end module Interactions_BSLL_ASC_CSLD__Form
