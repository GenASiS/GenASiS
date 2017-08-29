module RadiationMoments_ASC__Form

  !-- RadiationMoments_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use PhotonMoments_G__Form
  use NeutrinoMoments_G__Form
  use Sources_RM_CSL__Form
  use Sources_RM_ASC__Form
  use RadiationMoments_CSL__Form

  implicit none
  private
  
  type, public, extends ( Current_ASC_Template ) :: RadiationMoments_ASC_Form
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
    type ( Sources_RM_ASC_Form ), allocatable :: &
      Sources_ASC
    class ( Field_ASC_Template ), pointer :: &
      Interactions_ASC => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      RadiationMoments_CSL
    generic, public :: &
      RadiationMoments => RadiationMoments_CSL
    procedure, private, pass :: &
      PhotonMoments_G_CSL
    generic, public :: &
      PhotonMoments_G => PhotonMoments_G_CSL
    procedure, private, pass :: &
      NeutrinoMoments_G_CSL
    generic, public :: &
      NeutrinoMoments_G => NeutrinoMoments_G_CSL
    procedure, public, pass :: &
      SetInteractions
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type RadiationMoments_ASC_Form

contains


  subroutine Initialize &
               ( RMA, A, RadiationMomentsType, NameShortOption, &
                 RiemannSolverTypeOption, UseLimiterOption, &
                 AllocateSourcesOption, Velocity_U_UnitOption, &
                 MomentumDensity_U_UnitOption, MomentumDensity_D_UnitOption, &
                 EnergyDensityUnitOption, TemperatureUnitOption, &
                 EnergyUnitOption, MomentumUnitOption, &
                 AngularMomentumUnitOption, TimeUnitOption, &
                 LimiterParameterOption, IgnorabilityOption )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      RMA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      RadiationMomentsType
    character ( * ), intent ( in ), optional :: &
      NameShortOption, &
      RiemannSolverTypeOption
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterOption, &
      AllocateSourcesOption
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ), optional :: &
      Velocity_U_UnitOption, &
      MomentumDensity_U_UnitOption, &
      MomentumDensity_D_UnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      EnergyDensityUnitOption, &
      TemperatureUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption, &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      LimiterParameterOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

!    integer ( KDI ) :: &
!      iB  !-- iBoundary
    character ( LDL ) :: &
      NameShort
    logical ( KDL ) :: &
      AllocateSources

    if ( RMA % Type == '' ) &
      RMA % Type = 'a RadiationMoments_ASC'
    RMA % RadiationMomentsType = RadiationMomentsType

    RMA % RiemannSolverType = 'HLL'
    if ( present ( RiemannSolverTypeOption ) ) &
      RMA % RiemannSolverType = RiemannSolverTypeOption
    call PROGRAM_HEADER % GetParameter &
           ( RMA % RiemannSolverType, 'RiemannSolverType' )

    RMA % UseLimiter = .false.
    if ( present ( UseLimiterOption ) ) &
      RMA % UseLimiter = UseLimiterOption
    call PROGRAM_HEADER % GetParameter &
           ( RMA % UseLimiter, 'UseLimiter' )

    RMA % LimiterParameter = 2.0_KDR
    if ( present ( LimiterParameterOption ) ) &
      RMA % LimiterParameter = LimiterParameterOption
    call PROGRAM_HEADER % GetParameter &
           ( RMA % LimiterParameter, 'LimiterParameter' )

    if ( present ( EnergyDensityUnitOption ) ) &
      RMA % EnergyDensityUnit = EnergyDensityUnitOption
     if ( present ( TemperatureUnitOption ) ) &
      RMA % TemperatureUnit = TemperatureUnitOption
   if ( present ( Velocity_U_UnitOption ) ) &
      RMA % Velocity_U_Unit = Velocity_U_UnitOption
    if ( present ( MomentumDensity_U_UnitOption ) ) &
      RMA % MomentumDensity_U_Unit = MomentumDensity_U_UnitOption
    if ( present ( MomentumDensity_D_UnitOption ) ) &
      RMA % MomentumDensity_D_Unit = MomentumDensity_D_UnitOption

    ! if ( .not. allocated ( RMA % TallyInterior ) ) then
    !   select case ( trim ( RadiationMomentsType ) )
    !   case ( 'DUST' )
    !     allocate ( Tally_F_D_Form :: RMA % TallyInterior )
    !     allocate ( Tally_F_D_Form :: RMA % TallyTotal )
    !     allocate ( Tally_F_D_Form :: RMA % TallyChange )
    !     allocate ( RMA % TallyBoundaryLocal  ( A % nBoundaries ) )
    !     allocate ( RMA % TallyBoundaryGlobal ( A % nBoundaries ) )
    !     do iB = 1, A % nBoundaries 
    !       allocate &
    !         ( Tally_F_D_Form :: RMA % TallyBoundaryLocal  ( iB ) % Element )
    !       allocate &
    !         ( Tally_F_D_Form :: RMA % TallyBoundaryGlobal ( iB ) % Element )
    !     end do !-- iB
    !   case default
    !     call Show ( 'RadiationMomentsType not recognized', CONSOLE % ERROR )
    !     call Show ( 'RadiationMoments_ASC__Form', 'module', CONSOLE % ERROR )
    !     call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
    !     call PROGRAM_HEADER % Abort ( )
    !   end select !-- RadiationMomentsType   
    ! end if !-- allocated TallyInterior
      
    ! select type ( TI => RMA % TallyInterior )
    ! class is ( Tally_F_D_Form )
    !   call TI % Initialize &
    !          ( A, MassUnitOption = MassUnitOption, &
    !            EnergyUnitOption = EnergyUnitOption, &
    !            MomentumUnitOption = MomentumUnitOption, &
    !            AngularMomentumUnitOption = AngularMomentumUnitOption )
    ! end select !-- TI

    ! select type ( TT => RMA % TallyTotal )
    ! class is ( Tally_F_D_Form )
    !   call TT % Initialize &
    !          ( A, MassUnitOption = MassUnitOption, &
    !            EnergyUnitOption = EnergyUnitOption, &
    !            MomentumUnitOption = MomentumUnitOption, &
    !            AngularMomentumUnitOption = AngularMomentumUnitOption )
    ! end select !-- TT

    ! select type ( TC => RMA % TallyChange )
    ! class is ( Tally_F_D_Form )
    !   call TC % Initialize &
    !          ( A, MassUnitOption = MassUnitOption, &
    !            EnergyUnitOption = EnergyUnitOption, &
    !            MomentumUnitOption = MomentumUnitOption, &
    !            AngularMomentumUnitOption = AngularMomentumUnitOption )
    ! end select !-- TC

    ! do iB = 1, A % nBoundaries
    !   select type ( TB => RMA % TallyBoundaryLocal ( iB ) % Element )
    !   class is ( Tally_F_D_Form )
    !     call TB % Initialize &
    !            ( A, MassUnitOption = MassUnitOption, &
    !              EnergyUnitOption = EnergyUnitOption, &
    !              MomentumUnitOption = MomentumUnitOption, &
    !              AngularMomentumUnitOption = AngularMomentumUnitOption )
    !   end select !-- TB
    !   select type ( TB => RMA % TallyBoundaryGlobal ( iB ) % Element )
    !   class is ( Tally_F_D_Form )
    !     call TB % Initialize &
    !            ( A, MassUnitOption = MassUnitOption, &
    !              EnergyUnitOption = EnergyUnitOption, &
    !              MomentumUnitOption = MomentumUnitOption, &
    !              AngularMomentumUnitOption = AngularMomentumUnitOption )
    !   end select !-- TB
    ! end do !-- iB

    NameShort = 'Radiation'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call RMA % InitializeTemplate_ASC_C ( A, NameShort, IgnorabilityOption )

    call Show ( RMA % RadiationMomentsType, 'RadiationMomentsType', &
                RMA % IGNORABILITY )
    call Show ( RMA % RiemannSolverType, 'RiemannSolverType', &
                RMA % IGNORABILITY )
    call Show ( RMA % UseLimiter, 'UseLimiter', RMA % IGNORABILITY )
    call Show ( RMA % LimiterParameter, 'LimiterParameter', &
                RMA % IGNORABILITY )

    !-- Sources

    AllocateSources = .true.
    if ( present ( AllocateSourcesOption ) ) &
      AllocateSources = AllocateSourcesOption

    if ( AllocateSources ) then
      allocate ( RMA % Sources_ASC )
      associate ( SRMA => RMA % Sources_ASC )
      call SRMA % Initialize &
             ( RMA, NameShortOption = trim ( NameShort ) // '_Sources', &
               TimeUnitOption = TimeUnitOption, &
               IgnorabilityOption = IgnorabilityOption )
      select type ( SRMC => SRMA % Chart )
      class is ( Sources_RM_CSL_Form )
        select type ( RMC => RMA % Chart )
        class is ( RadiationMoments_CSL_Form )
          call RMC % SetSources ( SRMC )
        end select !-- RMC
      end select !-- SRMC
      end associate !-- SRMA
    end if

  end subroutine Initialize


  function RadiationMoments_CSL ( RMA ) result ( RM )

    class ( RadiationMoments_ASC_Form ), intent ( in ) :: &
      RMA
    class ( RadiationMomentsForm ), pointer :: &
      RM

    select type ( RMC => RMA % Chart )
    class is ( RadiationMoments_CSL_Form )
      RM => RMC % RadiationMoments ( )
    class default
      call Show ( 'RadiationMoments Chart type not recognized', &
                  CONSOLE % ERROR )
      call Show ( 'RadiationMoments_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'RadiationMoments_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function RadiationMoments_CSL


  function PhotonMoments_G_CSL ( RMA ) result ( PM )

    class ( RadiationMoments_ASC_Form ), intent ( in ) :: &
      RMA
    class ( PhotonMoments_G_Form ), pointer :: &
      PM

    select type ( RMC => RMA % Chart )
    class is ( RadiationMoments_CSL_Form )
      PM => RMC % PhotonMoments_G ( )
    class default
      call Show ( 'RadiationMoments Chart type not recognized', &
                  CONSOLE % ERROR )
      call Show ( 'RadiationMoments_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'PhotonMoments_G_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function PhotonMoments_G_CSL


  function NeutrinoMoments_G_CSL ( RMA ) result ( NM )

    class ( RadiationMoments_ASC_Form ), intent ( in ) :: &
      RMA
    class ( NeutrinoMoments_G_Form ), pointer :: &
      NM

    select type ( RMC => RMA % Chart )
    class is ( RadiationMoments_CSL_Form )
      NM => RMC % NeutrinoMoments_G ( )
    class default
      call Show ( 'RadiationMoments Chart type not recognized', &
                  CONSOLE % ERROR )
      call Show ( 'RadiationMoments_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'NeutrinoMoments_G_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function NeutrinoMoments_G_CSL


  subroutine SetInteractions ( RMA, IA )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      RMA
    class ( Field_ASC_Template ), intent ( in ), target :: &
      IA

    RMA % Interactions_ASC => IA

    select type ( RMC => RMA % Chart )
    class is ( RadiationMoments_CSL_Form )

    select type ( IC => IA % Chart )
    class is ( Field_CSL_Template )
    call RMC % SetInteractions ( IC )
    end select !-- I

    end select !-- RMC

  end subroutine SetInteractions


  impure elemental subroutine Finalize ( RMA )

    type ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      RMA

    nullify ( RMA % Interactions_ASC )

    if ( allocated ( RMA % Sources_ASC ) ) &
      deallocate ( RMA % Sources_ASC )

    call RMA % FinalizeTemplate_ASC_C ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( RadiationMoments_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( RadiationMoments_CSL_Form )
      call FC % Initialize &
             ( C, FA % NameShort, FA % RadiationMomentsType, &
               FA % RiemannSolverType, FA % UseLimiter, FA % Velocity_U_Unit, &
               FA % MomentumDensity_U_Unit, FA % MomentumDensity_D_Unit, &
               FA % EnergyDensityUnit, FA % TemperatureUnit, &
               FA % LimiterParameter, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module RadiationMoments_ASC__Form
