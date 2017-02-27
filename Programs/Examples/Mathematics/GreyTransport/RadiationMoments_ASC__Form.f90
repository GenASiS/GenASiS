module RadiationMoments_ASC__Form

  !-- RadiationMoments_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use PhotonMoments_Form
  use RadiationMoments_CSL__Form

  implicit none
  private
  
  type, public, extends ( Current_ASC_Template ) :: RadiationMoments_ASC_Form
    type ( MeasuredValueForm ) :: &
      EnergyDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    character ( LDF ) :: &
      RadiationMomentsType = ''
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
      PhotonMoments_CSL
    generic, public :: &
      PhotonMoments => PhotonMoments_CSL
    procedure, public, pass :: &
      SetInteractions
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type RadiationMoments_ASC_Form

contains


  subroutine Initialize &
               ( RMA, A, RadiationMomentsType, NameOutputOption, &
                 Velocity_U_UnitOption, MomentumDensity_U_UnitOption, &
                 MomentumDensity_D_UnitOption, EnergyDensityUnitOption, &
                 TemperatureUnitOption, EnergyUnitOption, MomentumUnitOption, &
                 AngularMomentumUnitOption )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      RMA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      RadiationMomentsType
    character ( * ), intent ( in ), optional :: &
      NameOutputOption
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ), optional :: &
      Velocity_U_UnitOption, &
      MomentumDensity_U_UnitOption, &
      MomentumDensity_D_UnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      EnergyDensityUnitOption, &
      TemperatureUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption

!    integer ( KDI ) :: &
!      iB  !-- iBoundary

    if ( RMA % Type == '' ) &
      RMA % Type = 'a RadiationMoments_ASC'
    RMA % RadiationMomentsType = RadiationMomentsType

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

!-- If there is a tally for RM
!    call RMA % InitializeTemplate_ASC &
!           ( A, NameOutputOption = NameOutputOption )

    call RMA % InitializeTemplate_ASC_C &
           ( A, NameOutputOption = NameOutputOption )

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


  function PhotonMoments_CSL ( RMA ) result ( PM )

    class ( RadiationMoments_ASC_Form ), intent ( in ) :: &
      RMA
    class ( PhotonMomentsForm ), pointer :: &
      PM

    select type ( RMC => RMA % Chart )
    class is ( RadiationMoments_CSL_Form )
      PM => RMC % PhotonMoments ( )
    class default
      call Show ( 'RadiationMoments Chart type not recognized', &
                  CONSOLE % ERROR )
      call Show ( 'RadiationMoments_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'RadiationMoments_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function PhotonMoments_CSL


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

    call RMA % FinalizeTemplate_ASC_C ( )

  end subroutine Finalize


  subroutine SetField ( FA, NameOutputOption )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      FA
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( RadiationMoments_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( RadiationMoments_CSL_Form )
      call FC % Initialize &
             ( C, FA % RadiationMomentsType, FA % Velocity_U_Unit, &
               FA % MomentumDensity_U_Unit, FA % MomentumDensity_D_Unit, &
               FA % EnergyDensityUnit, FA % TemperatureUnit, nValues, &
               NameOutputOption = NameOutputOption )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module RadiationMoments_ASC__Form
