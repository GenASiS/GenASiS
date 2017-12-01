module RadiationMoments_CSL__Form

  !-- RadiationMoments_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use Interactions_Template
  use RadiationMoments_Form
  use PhotonMoments_G__Form
  use PhotonMoments_S__Form
  use NeutrinoMoments_G__Form
  use NeutrinoMoments_S__Form
  use Sources_RM__Form
  use Sources_RM_CSL__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: RadiationMoments_CSL_Form
    real ( KDR ) :: &
      LimiterParameter
    type ( MeasuredValueForm ) :: &
      EnergyDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    logical ( KDL ) :: &
      UseLimiter
    character ( LDF ) :: &
      RadiationMomentsType = '', &
      RiemannSolverType = ''
    class ( Sources_RM_CSL_Form ), pointer :: &
      Sources_CSL => null ( )
    class ( Field_CSL_Template ), pointer :: &
      Interactions_CSL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      RadiationMoments
    procedure, public, pass :: &
      PhotonMoments_G
    procedure, public, pass :: &
      PhotonMoments_S
    procedure, public, pass :: &
      NeutrinoMoments_G
    procedure, public, pass :: &
      NeutrinoMoments_S
    procedure, public, pass :: &
      SetSources
    procedure, public, pass :: &
      SetInteractions
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type RadiationMoments_CSL_Form

contains


  subroutine Initialize &
               ( RMC, C, NameShort, RadiationMomentsType, RiemannSolverType, &
                 UseLimiter, Velocity_U_Unit, MomentumDensity_U_Unit, &
                 MomentumDensity_D_Unit, EnergyDensityUnit, TemperatureUnit, &
                 LimiterParameter, nValues, IgnorabilityOption )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      RMC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      RadiationMomentsType, &
      RiemannSolverType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit, &
      TemperatureUnit
    real ( KDR ), intent ( in ) :: &
      LimiterParameter
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( RMC % Type == '' ) &
      RMC % Type = 'a RadiationMoments_CSL'
    RMC % RadiationMomentsType = RadiationMomentsType
    RMC % RiemannSolverType    = RiemannSolverType
    RMC % UseLimiter           = UseLimiter
    RMC % LimiterParameter     = LimiterParameter

    RMC % EnergyDensityUnit      = EnergyDensityUnit
    RMC % TemperatureUnit        = TemperatureUnit
    RMC % Velocity_U_Unit        = Velocity_U_Unit
    RMC % MomentumDensity_U_Unit = MomentumDensity_U_Unit
    RMC % MomentumDensity_D_Unit = MomentumDensity_D_Unit

    call RMC % InitializeTemplate_CSL &
           ( C, NameShort, nValues, IgnorabilityOption )

  end subroutine Initialize


  function RadiationMoments ( RMC ) result ( RM )

    class ( RadiationMoments_CSL_Form ), intent ( in ), target :: &
      RMC
    class ( RadiationMomentsForm ), pointer :: &
      RM
      
    class ( StorageForm ), pointer :: &
      Field

    RM => null ( )

    Field => RMC % Field
    select type ( Field )
    class is ( RadiationMomentsForm )
    RM => Field
    end select !-- Field

  end function RadiationMoments


  function PhotonMoments_G ( RMC ) result ( PM )

    class ( RadiationMoments_CSL_Form ), intent ( in ), target :: &
      RMC
    class ( PhotonMoments_G_Form ), pointer :: &
      PM
      
    class ( StorageForm ), pointer :: &
      Field

    PM => null ( )

    Field => RMC % Field
    select type ( Field )
    class is ( PhotonMoments_G_Form )
    PM => Field
    end select !-- Field

  end function PhotonMoments_G


  function PhotonMoments_S ( RMC ) result ( PM )

    class ( RadiationMoments_CSL_Form ), intent ( in ), target :: &
      RMC
    class ( PhotonMoments_S_Form ), pointer :: &
      PM
      
    class ( StorageForm ), pointer :: &
      Field

    PM => null ( )

    Field => RMC % Field
    select type ( Field )
    class is ( PhotonMoments_S_Form )
    PM => Field
    end select !-- Field

  end function PhotonMoments_S


  function NeutrinoMoments_G ( RMC ) result ( NM )

    class ( RadiationMoments_CSL_Form ), intent ( in ), target :: &
      RMC
    class ( NeutrinoMoments_G_Form ), pointer :: &
      NM

    class ( VariableGroupForm ), pointer :: &
      Field

    NM => null ( )

    Field => RMC % Field
    select type ( Field )
    class is ( NeutrinoMoments_G_Form )
    NM => Field
    end select !-- Field

  end function NeutrinoMoments_G


  function NeutrinoMoments_S ( RMC ) result ( NM )

    class ( RadiationMoments_CSL_Form ), intent ( in ), target :: &
      RMC
    class ( NeutrinoMoments_S_Form ), pointer :: &
      NM

    class ( VariableGroupForm ), pointer :: &
      Field

    NM => null ( )

    Field => RMC % Field
    select type ( Field )
    class is ( NeutrinoMoments_S_Form )
    NM => Field
    end select !-- Field

  end function NeutrinoMoments_S


  subroutine SetSources ( RMC, SRMC )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      RMC
    class ( Sources_RM_CSL_Form ), intent ( in ), target :: &
      SRMC

    class ( RadiationMomentsForm ), pointer :: &
      RM

    RMC % Sources_CSL => SRMC

    RM => RMC % RadiationMoments ( )
    select type ( SRM => SRMC % Field )
    class is ( Sources_RM_Form )
      call RM % SetSources ( SRM )
    end select !-- SRM

    nullify ( RM )

  end subroutine SetSources


  subroutine SetInteractions ( RMC, IC )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      RMC
    class ( Field_CSL_Template ), intent ( in ), target :: &
      IC

    class ( RadiationMomentsForm ), pointer :: &
      RM

    RMC % Interactions_CSL => IC

    RM => RMC % RadiationMoments ( )
    select type ( I => IC % Field )
    class is ( InteractionsTemplate )
    call RM % SetInteractions ( I )
    end select !-- I

    nullify ( RM )

  end subroutine SetInteractions


  impure elemental subroutine Finalize ( RMC )

    type ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      RMC

    nullify ( RMC % Interactions_CSL )
    nullify ( RMC % Sources_CSL )

    call RMC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % RadiationMomentsType ) )
    case ( 'GENERIC' )
      allocate ( RadiationMomentsForm :: FC % Field )
      select type ( RM => FC % Field )
      type is ( RadiationMomentsForm )
        call RM % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_U_Unit, &
                 FC % MomentumDensity_D_Unit, FC % EnergyDensityUnit, &
                 FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call RM % SetPrimitiveConserved ( )
        call RM % SetOutput ( FC % FieldOutput )
      end select !-- RM
    case ( 'PHOTONS_GREY' )
      allocate ( PhotonMoments_G_Form :: FC % Field )
      select type ( PM => FC % Field )
      type is ( PhotonMoments_G_Form )
        call PM % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_U_Unit, &
                 FC % MomentumDensity_D_Unit, FC % EnergyDensityUnit, &
                 FC % TemperatureUnit, FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call PM % SetPrimitiveConserved ( )
        call PM % SetOutput ( FC % FieldOutput )
      end select !-- PM
    case ( 'PHOTONS_SPECTRAL' )
      allocate ( PhotonMoments_S_Form :: FC % Field )
      select type ( PM => FC % Field )
      type is ( PhotonMoments_S_Form )
        call PM % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_U_Unit, &
                 FC % MomentumDensity_D_Unit, FC % EnergyDensityUnit, &
                 FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call PM % SetPrimitiveConserved ( )
        call PM % SetOutput ( FC % FieldOutput )
      end select !-- PM
    case ( 'NEUTRINOS_GREY' )
      allocate ( NeutrinoMoments_G_Form :: FC % Field )
      select type ( NM => FC % Field )
      type is ( NeutrinoMoments_G_Form )
        call Show( 'NEUTRINOS_GREY' )
        call NM % InitializeAllocate_NM &
               ( FC % RadiationMomentsType, FC % RiemannSolverType, &
                 FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_U_Unit, &
                 FC % MomentumDensity_D_Unit, FC % EnergyDensityUnit, &
                 FC % TemperatureUnit, FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call NM % SetPrimitiveConserved ( )
        call NM % SetOutput ( FC % FieldOutput )
      end select !-- NM
    case ( 'NEUTRINOS_SPECTRAL' )
      allocate ( NeutrinoMoments_S_Form :: FC % Field )
      select type ( NM => FC % Field )
      type is ( NeutrinoMoments_S_Form )
        call NM % Initialize &
               ( FC % RiemannSolverType, FC % UseLimiter, &
                 FC % Velocity_U_Unit, FC % MomentumDensity_U_Unit, &
                 FC % MomentumDensity_D_Unit, FC % EnergyDensityUnit, &
                 FC % LimiterParameter, FC % nValues, &
                 NameOption = FC % NameShort )
        call NM % SetPrimitiveConserved ( )
        call NM % SetOutput ( FC % FieldOutput )
      end select !-- NM
    case default
      call Show ( 'RadiationMomentsType not recognized', CONSOLE % ERROR )
      call Show ( FC % RadiationMomentsType, 'RadiationMomentsType', &
                  CONSOLE % ERROR )
      call Show ( 'RadiationMoments_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- RadiationMomentsType

  end subroutine SetField


end module RadiationMoments_CSL__Form
