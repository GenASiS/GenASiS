module InteractionsExamples_BSLL_ASC_CSLD__Form

  !-- InteractionsExamples_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Form

  use GenASiS
  use Interactions_C__Form
  use Interactions_MWV_1__Form
  use Interactions_MWV_2__Form
  use Interactions_MWV_3__Form
  use InteractionsExamples_ASC__Form

  implicit none
  private

  type, public, extends ( Interactions_BSLL_ASC_CSLD_Template ) :: &
    InteractionsExamples_BSLL_ASC_CSLD_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Set_C_Spectral
    procedure, private, pass :: &
      Set_MWV_1_Spectral
    procedure, private, pass :: &
      Set_MWV_2_Spectral
    procedure, private, pass :: &
      Set_MWV_3_Spectral
    generic, public :: &
      Set_MWV_Spectral &
        => Set_MWV_1_Spectral, Set_MWV_2_Spectral, Set_MWV_3_Spectral
    procedure, public, pass :: &
      AllocateField
  end type InteractionsExamples_BSLL_ASC_CSLD_Form

contains


  subroutine Initialize &
               ( IB, B, InteractionsType, Units, NameShortOption )

    class ( InteractionsExamples_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ) :: &
      B
    character ( * ), intent ( in ) :: &
      InteractionsType
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
    character ( * ), intent ( in ), optional :: &
      NameShortOption

    if ( IB % Type == '' ) &
      IB % Type = 'an InteractionsExamples_BSLL_ASC_CSLD'

    call IB % InitializeTemplate_I_BSLL_ASC_CSLD &
           ( B, InteractionsType, Units, NameShortOption )

  end subroutine Initialize


  subroutine Set_C_Spectral ( IB, FA, OpacityAbsorption )

    class ( InteractionsExamples_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
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
      type is ( Interactions_C_Form )
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

  end subroutine Set_C_Spectral


  subroutine Set_MWV_1_Spectral ( IB, FA, SpecificOpacity )

    class ( InteractionsExamples_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity

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
    associate &
      (    Energy => GF % Value ( :, GF % CENTER_U ( 1 ) ), &
        d3_Energy => GF % Value ( :, GF % VOLUME ) )

    do iF = 1, B % nFibers
      I => IB % Interactions ( iF )
      select type ( I )
      type is ( Interactions_MWV_1_Form )
        call I % Set &
               ( Fluid = Fluid, &
                 Energy = Energy, &
                 d3_Energy = d3_Energy, &
                 SpecificOpacity = SpecificOpacity, &
                 iBaseCell = B % iaBaseCell ( iF ) )
      end select !-- I
      nullify ( I )
    end do !-- iF

    end associate !-- Energy, etc.
    end select !-- B
    nullify ( GF )

  end subroutine Set_MWV_1_Spectral


  subroutine Set_MWV_2_Spectral &
               ( IB, FA, SpecificOpacity, SpecificOpacityFloor, EnergyMax )

    class ( InteractionsExamples_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      SpecificOpacityFloor, &
      EnergyMax

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
    associate &
      (    Energy => GF % Value ( :, GF % CENTER_U ( 1 ) ), &
        d3_Energy => GF % Value ( :, GF % VOLUME ) )

    do iF = 1, B % nFibers
      I => IB % Interactions ( iF )
      select type ( I )
      type is ( Interactions_MWV_2_Form )
        call I % Set &
               ( Fluid = Fluid, &
                 Energy = Energy, &
                 d3_Energy = d3_Energy, &
                 SpecificOpacity = SpecificOpacity, &
                 SpecificOpacityFloor = SpecificOpacityFloor, &
                 EnergyMax = EnergyMax, &
                 iBaseCell = B % iaBaseCell ( iF ) )
      end select !-- I
      nullify ( I )
    end do !-- iF

    end associate !-- Energy, etc.
    end select !-- B
    nullify ( GF )

  end subroutine Set_MWV_2_Spectral


  subroutine Set_MWV_3_Spectral &
               ( IB, FA, SpecificOpacity, SpecificOpacityFloor, EnergyMax, &
                 TemperatureScale )

    class ( InteractionsExamples_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      SpecificOpacityFloor, &
      EnergyMax, &
      TemperatureScale

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
    associate &
      (    Energy => GF % Value ( :, GF % CENTER_U ( 1 ) ), &
        d3_Energy => GF % Value ( :, GF % VOLUME ) )

    do iF = 1, B % nFibers
      I => IB % Interactions ( iF )
      select type ( I )
      type is ( Interactions_MWV_3_Form )
        call I % Set &
               ( Fluid = Fluid, &
                 Energy = Energy, &
                 d3_Energy = d3_Energy, &
                 SpecificOpacity = SpecificOpacity, &
                 SpecificOpacityFloor = SpecificOpacityFloor, &
                 EnergyMax = EnergyMax, &
                 TemperatureScale = TemperatureScale, &
                 iBaseCell = B % iaBaseCell ( iF ) )
      end select !-- I
      nullify ( I )
    end do !-- iF

    end associate !-- Energy, etc.
    end select !-- B
    nullify ( GF )

  end subroutine Set_MWV_3_Spectral


  subroutine AllocateField ( IB, iF )

    class ( InteractionsExamples_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    integer ( KDI ), intent ( in ) :: &
      iF  !-- iFiber

    allocate ( InteractionsExamples_ASC_Form :: &
                 IB % Fiber % Atlas ( iF ) % Element )

  end subroutine AllocateField


end module InteractionsExamples_BSLL_ASC_CSLD__Form
