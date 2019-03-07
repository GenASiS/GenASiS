module Interactions_BSLL_ASC_CSLD__Template

  !-- Interactions_Generic_Spectral_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use Fluids
  use Interactions_Template
  use Interactions_ASC__Template

  implicit none
  private

  type, public, extends ( Field_BSLL_ASC_CSLD_Template ), abstract :: &
    Interactions_BSLL_ASC_CSLD_Template
      type ( StressEnergyUnitsForm ) :: &
        UnitsSpectral
      class ( StressEnergyUnitsForm ), pointer :: &
        Units => null ( )
      character ( LDF ) :: &
        InteractionsType = ''
  contains
    procedure ( I ), public, pass, deferred :: &
      Initialize
    procedure, public, pass :: &
      InitializeTemplate_I_BSLL_ASC_CSLD
    procedure, public, pass :: &
      Interactions
    procedure, public, pass :: &
      ComputeTimeScale
    procedure, public, pass :: &
      FinalizeTemplate_I_BSLL_ASC_CSLD
    procedure, private, pass :: &
      SetField
    procedure ( AF ), public, pass, deferred :: &
      AllocateField
  end type Interactions_BSLL_ASC_CSLD_Template

  abstract interface

    subroutine I ( IB, B, InteractionsType, Units, NameShortOption )
      use Basics
      use Mathematics
      use StressEnergyBasics
      import Interactions_BSLL_ASC_CSLD_Template
      class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
        IB
      class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ) :: &
        B
      character ( * ), intent ( in ) :: &
        InteractionsType
      class ( StressEnergyUnitsForm ), intent ( in ) :: &
        Units
      character ( * ), intent ( in ), optional :: &
        NameShortOption
    end subroutine I

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
               ( IB, B, InteractionsType, Units, NameShortOption )

    class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      IB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ) :: &
      B
    character ( * ), intent ( in )  :: &
      InteractionsType
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    character ( * ), intent ( in ), optional :: &
      NameShortOption

    type ( MeasuredValueForm ) :: &
      ParticleEnergyUnit
    character ( LDL ) :: &
      NameShort

    if ( IB % Type == '' ) &
      IB % Type = 'an Interactions_BSLL_ASC_CSLD'
    IB % InteractionsType = InteractionsType

    IB % Units => Units

    associate ( AF => B % FiberMaster )
    select type ( CF => AF % Chart )
    class is ( Chart_SLL_Form )
      ParticleEnergyUnit  =  CF % CoordinateUnit ( 1 )
    end select !--CF
    end associate !-- AF
    IB % UnitsSpectral  =  Units
    IB % UnitsSpectral % EnergyDensity &
      =  Units % EnergyDensity  *  ParticleEnergyUnit ** (-3)
    IB % UnitsSpectral % MomentumDensity_U &
      =  Units % MomentumDensity_U  *  ParticleEnergyUnit ** (-3)
    IB % UnitsSpectral % MomentumDensity_D &
      =  Units % MomentumDensity_D  *  ParticleEnergyUnit ** (-3)

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


  subroutine ComputeTimeScale ( IB, RB, TimeScale )

    class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( in ) :: &
      IB
    class ( Current_BSLL_ASC_CSLD_Template ), intent ( in ) :: &
      RB
    real ( KDR ), intent ( out ) :: &
      TimeScale

    integer ( KDI ) :: &
      iF  !-- iFiber
    class ( CurrentTemplate ), pointer :: &
      R
    class ( StorageForm ), pointer :: &
      I

    associate ( MS => RB % Bundle_SLL_ASC_CSLD )

    TimeScale = huge ( 1.0_KDR )

    do iF = 1, IB % nFibers
      associate ( iBC => MS % iaBaseCell ( iF ) )

      I => IB % FieldFiber ( iF )
      R => RB % CurrentFiber ( iF )

      select type ( I )
      class is ( InteractionsTemplate )

      associate ( F => I % Fluid )
      select type ( SF => F % Sources )
      class is ( Sources_F_Form )

        call I % ComputeTimeScale ( R )
        TimeScale = min ( TimeScale, SF % Value ( iBC, SF % RADIATION_TIME ) )

      end select !-- SF
      end associate !-- F
      end select !-- I

      end associate !-- iBC
    end do !-- iF

    end associate !-- MS
    nullify ( R, I )

  end subroutine ComputeTimeScale


  impure elemental subroutine FinalizeTemplate_I_BSLL_ASC_CSLD ( IB )

    class ( Interactions_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      IB

    nullify ( IB % Units )

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
               Units = FB % UnitsSpectral, NameShortOption = FB % NameShort )

      end select !-- AF
      end select !-- IA
    end do !-- iF
    
    end associate !-- FBF
    end select !-- B

  end subroutine SetField


end module Interactions_BSLL_ASC_CSLD__Template
