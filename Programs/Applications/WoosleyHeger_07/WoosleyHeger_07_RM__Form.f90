module WoosleyHeger_07_RM__Form

  !-- WoosleyHeger_07_RadiationMoments_Form

  use GenASiS
  use WoosleyHeger_07__Template

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Template ) :: WoosleyHeger_07_RM_Form
  contains
    procedure, private, pass :: &
      Initialize_WH
    generic, public :: &
      Initialize => Initialize_WH
    final :: &
      Finalize
  end type WoosleyHeger_07_RM_Form

    private :: &
      InitializeRadiationCentralCore, &
      SetProblem

      private :: &
        InitializeInteractions

    character ( LDF ), private :: &
      InteractionsType

contains


  subroutine Initialize_WH ( WH, MomentsType, Name )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    if ( WH % Type == '' ) &
      WH % Type = 'a WoosleyHeger_07_RM'

    call InitializeRadiationCentralCore ( WH, MomentsType, Name )
    call SetProblem ( WH, MomentsType )

  end subroutine Initialize_WH


  subroutine Finalize ( WH )

    type ( WoosleyHeger_07_RM_Form ), intent ( inout ) :: &
      WH

  end subroutine Finalize


  subroutine InitializeRadiationCentralCore ( WH, MomentsType, Name )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in ) :: &
      MomentsType, &
      Name

    real ( KDR ) :: &
      MaxEnergy, &
      MinWidthEnergy, &
      EnergyScale
    character ( LDL ) :: &
      GeometryType
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType

    GeometryType = 'NEWTONIAN'
    call PROGRAM_HEADER % GetParameter ( GeometryType, 'GeometryType' )

    allocate ( RadiationName ( 2 ) )
    allocate ( RadiationType ( 2 ) )
    RadiationName  =  [ 'Nu_E   ', 'NuBar_E' ]
    RadiationType  =  [ 'NEUTRINOS_E    ', 'NEUTRINOS_E_BAR' ]

    MinWidthEnergy  =  2.0_KDR    *  UNIT % MEGA_ELECTRON_VOLT  !-- Geometric
    MaxEnergy       =  310.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT

    EnergyScale  =  10.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT  !-- Compactified

    call WH % Initialize &
           ( RadiationName = RadiationName, RadiationType = RadiationType, &
             MomentsType = MomentsType, FluidType = 'HEAVY_NUCLEUS', &
             GeometryType = GeometryType, Name = Name, &
             ShockThresholdOption = 1.0_KDR, &
             MinWidthEnergyOption = MinWidthEnergy, &
             MaxEnergyOption = MaxEnergy, &
             EnergyScaleOption = EnergyScale, &
             nWriteOption = 30 )

  end subroutine InitializeRadiationCentralCore


  subroutine SetProblem ( WH, MomentsType )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ) :: &
      WH
    character ( * ), intent ( in )  :: &
      MomentsType

    call InitializeInteractions ( WH, MomentsType )
    call WH % SetFluid ( )

  end subroutine SetProblem


  subroutine InitializeInteractions ( WH, MomentsType )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ) :: &
      WH
    character ( * ), intent ( in )  :: &
      MomentsType

    integer ( KDI ) :: &
      iR  !-- iRadiation
    logical ( KDL ) :: &
      Include_NES, &
      IncludePairs

    InteractionsType = 'O_CONNOR_OTT'
    call PROGRAM_HEADER % GetParameter ( InteractionsType, 'InteractionsType' )
    
    Include_NES   =  .true.
    IncludePairs  =  .true.
    call PROGRAM_HEADER % GetParameter ( Include_NES,   'Include_NES' )
    call PROGRAM_HEADER % GetParameter ( IncludePairs, 'IncludePairs' )

    select type ( I => WH % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

      ! select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
      ! class is ( RadiationMoments_ASC_Form )

      ! select type ( FA => I % Current_ASC )
      ! class is ( Fluid_ASC_Form )

      ! select type ( PS => I % PositionSpace )
      ! class is ( Atlas_SC_Form )

      ! allocate ( InteractionsExamples_ASC_Form :: MW % Interactions_ASC )
      ! select type ( IA => MW % Interactions_ASC )
      ! class is ( InteractionsExamples_ASC_Form )
      ! call IA % Initialize &
      !        ( PS, InteractionsType = InteractionsType, &
      !          MomentsType = MomentsType, Units = MW % Units )
      ! select case ( trim ( InteractionsType ) )
      ! case ( 'MARSHAK_WAVE_VAYTET_1' )
      !   call IA % Set_MWV_Grey &
      !          ( FA, SpecificOpacity = MW % SpecificOpacity )
      ! case ( 'MARSHAK_WAVE_VAYTET_2' )
      !   call IA % Set_MWV_Grey &
      !          ( FA, SpecificOpacity = MW % SpecificOpacity, &
      !            EnergyMax = MW % EnergyMax )
      ! case ( 'MARSHAK_WAVE_VAYTET_3' )
      !   call IA % Set_MWV_Grey &
      !          ( FA, SpecificOpacity = MW % SpecificOpacity, &
      !            EnergyMax = MW % EnergyMax, &
      !            TemperatureScale = MW % Temperature )
      ! end select !-- InteractionsType
      ! call RMA % SetInteractions ( IA )
      ! end select !-- IA

      ! end select !-- PS
      ! end select !-- FA
      ! end select !-- RMA

    class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

      select type ( FA => I % Current_ASC )
      class is ( Fluid_ASC_Form )

      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )

      allocate ( InteractionsNeutrinos_BSLL_ASC_CSLD_Form :: &
                   WH % Interactions_BSLL_ASC_CSLD )
      select type ( IB => WH % Interactions_BSLL_ASC_CSLD )
      class is ( InteractionsNeutrinos_BSLL_ASC_CSLD_Form )
      call IB % Initialize &
             ( MS, InteractionsType = InteractionsType, Units = WH % Units )
      select case ( trim ( InteractionsType ) )
      case ( 'O_CONNOR_OTT' )
        call IB % Set ( FA )
      end select !-- InteractionsType

      do iR = 1, I % N_CURRENTS_1D
        select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( iR ) % Element )
        class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
          call RMB % SetInteractions ( IB )
        end select !-- RMB
      end do

      end select !-- IB

      end select !-- MS
      end select !-- FA

    end select !-- Integrator

  end subroutine InitializeInteractions


end module WoosleyHeger_07_RM__Form
