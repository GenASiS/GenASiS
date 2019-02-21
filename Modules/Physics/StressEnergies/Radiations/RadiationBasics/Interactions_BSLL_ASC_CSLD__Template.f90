module Interactions_BSLL_ASC_CSLD__Template

  !-- Interactions_Generic_Spectral_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Form

  use Basics
  use Mathematics
  use Interactions_Template
  use Interactions_ASC__Template

  implicit none
  private

  type, public, extends ( Field_BSLL_ASC_CSLD_Template ), abstract :: &
    Interactions_BSLL_ASC_CSLD_Template
      type ( MeasuredValueForm ) :: &
        LengthUnit, &
        EnergyDensityUnit, &
        TemperatureUnit
      character ( LDF ) :: &
        InteractionsType = ''
  contains
    procedure, public, pass :: &
      InitializeTemplate_I_BSLL_ASC_CSLD
    procedure, public, pass :: &
      Interactions
    procedure, public, pass :: &
      FinalizeTemplate_I_BSLL_ASC_CSLD
    procedure, private, pass :: &
      SetField
    procedure ( AF ), public, pass, deferred :: &
      AllocateField
  end type Interactions_BSLL_ASC_CSLD_Template

  abstract interface
    subroutine AF ( IB, iF )
      use Basics
      import Interactions_BSLL_ASC_CSLD_Template
      class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
        IB
      integer ( KDI ), intent ( in ) :: &
        iF  !-- iFiber
    end subroutine AF
  end interface

contains


  subroutine InitializeTemplate_I_BSLL_ASC_CSLD &
               ( IB, B, InteractionsType, NameShortOption, LengthUnitOption, &
                 EnergyDensityUnitOption, TemperatureUnitOption )

    class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      IB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in )  :: &
      InteractionsType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      LengthUnitOption, &
      EnergyDensityUnitOption, &
      TemperatureUnitOption

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
    if ( present ( TemperatureUnitOption ) ) &
      IB % TemperatureUnit = TemperatureUnitOption

    NameShort = 'Interactions'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call IB % InitializeTemplate_BSLL_ASC_CSLD ( B, NameShort )

    call Show ( IB % InteractionsType, 'InteractionsType', IB % IGNORABILITY )

    ! nullify ( GF )

  end subroutine InitializeTemplate_I_BSLL_ASC_CSLD


  function Interactions ( IB, iFiber ) result ( IF )

    class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
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


  impure elemental subroutine FinalizeTemplate_I_BSLL_ASC_CSLD ( IB )

    class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      IB

    call IB % FinalizeTemplate_BSLL_ASC_CSLD ( )

  end subroutine FinalizeTemplate_I_BSLL_ASC_CSLD


  subroutine SetField ( FB )

    class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB

    integer ( KDI ) :: &
      iF  !-- iFiber

    select type ( B => FB % Bundle )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    allocate ( FB % Fiber )
    associate ( FBF => FB % Fiber )
    call FBF % Initialize ( B % nFibers )

    do iF = 1, B % nFibers

      call FB % AllocateField ( iF )

      select type ( IA => FBF % Atlas ( iF ) % Element )
      class is ( Interactions_ASC_Template )

      select type ( AF => B % Fiber % Atlas ( iF ) % Element )
      class is ( Atlas_SC_Form )

      call IA % Initialize &
             ( AF, FB % InteractionsType, MomentsType = 'SPECTRAL', &
               NameShortOption = FB % NameShort, & 
               LengthUnitOption = FB % LengthUnit, &
               EnergyDensityUnitOption = FB % EnergyDensityUnit, &
               TemperatureUnitOption = FB % TemperatureUnit )

      end select !-- AF
      end select !-- IA
    end do !-- iF
    
    end associate !-- FBF
    end select !-- B

  end subroutine SetField


end module Interactions_BSLL_ASC_CSLD__Template
