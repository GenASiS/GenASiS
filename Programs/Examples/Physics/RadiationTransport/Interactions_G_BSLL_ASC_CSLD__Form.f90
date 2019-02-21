module Interactions_G_BSLL_ASC_CSLD__Form

  !-- Interactions_Generic_Spectral_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Form

  use GenASiS
  use Interactions_G__Form
  use Interactions_G_ASC__Form

  implicit none
  private

  type, public, extends ( Interactions_BSLL_ASC_CSLD_Template ) :: &
    Interactions_G_BSLL_ASC_CSLD_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Set => Set_G_BSLL_ASC_CSLD
    procedure, public, pass :: &
      AllocateField
  end type Interactions_G_BSLL_ASC_CSLD_Form

contains


  subroutine Initialize &
               ( IB, B, InteractionsType, NameShortOption, LengthUnitOption, &
                 EnergyDensityUnitOption, TemperatureUnitOption )

    class ( Interactions_G_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in ) :: &
      InteractionsType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      LengthUnitOption, &
      EnergyDensityUnitOption, &
      TemperatureUnitOption

    if ( IB % Type == '' ) &
      IB % Type = 'an Interactions_G_BSLL_ASC_CSLD'

    call IB % InitializeTemplate_I_BSLL_ASC_CSLD &
           ( B, InteractionsType, NameShortOption, LengthUnitOption, &
             EnergyDensityUnitOption, TemperatureUnitOption )

  end subroutine Initialize


  subroutine Set_G_BSLL_ASC_CSLD ( IB, FA, OpacityAbsorption )

    class ( Interactions_G_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA
    real ( KDR ), intent ( in ) :: &
      OpacityAbsorption

    integer ( KDI ) :: &
      iF  !-- iFiber
    class ( GeometryFlatForm ), pointer :: &
      GF
    class ( Fluid_P_I_Form ), pointer :: &
      Fluid
    class ( InteractionsTemplate ), pointer :: &
      I

    select type ( B => IB % Bundle )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    Fluid => FA % Fluid_P_I ( )

    GF => B % GeometryFiber ( )
    associate ( Energy => GF % Value ( :, GF % CENTER_U ( 1 ) ) )

    do iF = 1, B % nFibers
      I => IB % Interactions ( iF )
      select type ( I )
      type is ( Interactions_G_Form )
        call I % Set &
               ( Fluid = Fluid, &
                 Energy = Energy, &
                 OpacityAbsorption = OpacityAbsorption, &
                 iBaseCell = B % iaBaseCell ( iF ) )
      end select !-- I
      nullify ( I )
    end do !-- iF

    end associate !-- Energy
    end select !-- B
    nullify ( GF )

  end subroutine Set_G_BSLL_ASC_CSLD


  subroutine AllocateField ( IB, iF )

    class ( Interactions_G_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    integer ( KDI ), intent ( in ) :: &
      iF  !-- iFiber

    allocate ( Interactions_G_ASC_Form :: &
                 IB % Fiber % Atlas ( iF ) % Element )

  end subroutine AllocateField


end module Interactions_G_BSLL_ASC_CSLD__Form
