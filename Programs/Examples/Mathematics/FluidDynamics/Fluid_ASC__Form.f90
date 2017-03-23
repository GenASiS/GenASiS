module Fluid_ASC__Form

  !-- Fluid_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use Fluid_D__Form
  use Fluid_P_P__Form
  use Fluid_P_NR__Form
  use Fluid_P_MHN__Form
  use Tally_F_D__Form
  use Tally_F_P__Form
  use FluidFeatures_CSL__Form
  use FluidFeatures_ASC__Form
  use Fluid_CSL__Form
  
  implicit none
  private
  
  type, public, extends ( Current_ASC_Template ) :: Fluid_ASC_Form
    real ( KDR ) :: &
      LimiterParameter
    type ( MeasuredValueForm ) :: &
      MassDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      VelocityUnit
    logical ( KDL ) :: &
      UseLimiter
    character ( LDF ) :: &
      FluidType         = '', &
      RiemannSolverType = ''
    type ( FluidFeatures_ASC_Form ), allocatable :: &
      Features_ASC
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Fluid_D_CSL
    generic, public :: &
      Fluid_D => Fluid_D_CSL
    procedure, private, pass :: &
      Fluid_P_P_CSL
    generic, public :: &
      Fluid_P_P => Fluid_P_P_CSL
    procedure, private, pass :: &
      Fluid_P_NR_CSL
    generic, public :: &
      Fluid_P_NR => Fluid_P_NR_CSL
    procedure, private, pass :: &
      Fluid_P_MHN_CSL
    generic, public :: &
      Fluid_P_MHN => Fluid_P_MHN_CSL
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Fluid_ASC_Form

contains


  subroutine Initialize &
               ( FA, A, FluidType, NameShortOption, RiemannSolverTypeOption, &
                 UseLimiterOption, VelocityUnitOption, MassDensityUnitOption, &
                 EnergyDensityUnitOption, TemperatureUnitOption, &
                 MassUnitOption, EnergyUnitOption, MomentumUnitOption, &
                 AngularMomentumUnitOption, LimiterParameterOption, &
                 ShockThresholdOption, IgnorabilityOption )

    class ( Fluid_ASC_Form ), intent ( inout ) :: &
      FA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      FluidType
    character ( * ), intent ( in ), optional :: &
      NameShortOption, &
      RiemannSolverTypeOption
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterOption
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ), optional :: &
      VelocityUnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      MassDensityUnitOption, &
      EnergyDensityUnitOption, &
      TemperatureUnitOption, &
      MassUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption
    real ( KDR ), intent ( in ), optional :: &
      LimiterParameterOption, &
      ShockThresholdOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iB  !-- iBoundary
    character ( LDL ) :: &
      NameShort

    if ( FA % Type == '' ) &
      FA % Type = 'a Fluid_ASC'
    FA % FluidType = FluidType

    FA % RiemannSolverType = 'HLLC'
    if ( present ( RiemannSolverTypeOption ) ) &
      FA % RiemannSolverType = RiemannSolverTypeOption
    call PROGRAM_HEADER % GetParameter &
           ( FA % RiemannSolverType, 'RiemannSolverType' )

    FA % UseLimiter = .true.
    if ( present ( UseLimiterOption ) ) &
      FA % UseLimiter = UseLimiterOption
    call PROGRAM_HEADER % GetParameter &
           ( FA % UseLimiter, 'UseLimiter' )

    FA % LimiterParameter = 1.4_KDR
    if ( present ( LimiterParameterOption ) ) &
      FA % LimiterParameter = LimiterParameterOption
    call PROGRAM_HEADER % GetParameter &
           ( FA % LimiterParameter, 'LimiterParameter' )

    if ( present ( MassDensityUnitOption ) ) &
      FA % MassDensityUnit = MassDensityUnitOption
    if ( present ( EnergyDensityUnitOption ) ) &
      FA % EnergyDensityUnit = EnergyDensityUnitOption
    if ( present ( TemperatureUnitOption ) ) &
      FA % TemperatureUnit = TemperatureUnitOption
    if ( present ( VelocityUnitOption ) ) &
      FA % VelocityUnit = VelocityUnitOption

    if ( .not. allocated ( FA % TallyInterior ) ) then
      select case ( trim ( FluidType ) )
      case ( 'DUST' )
        allocate ( Tally_F_D_Form :: FA % TallyInterior )
        allocate ( Tally_F_D_Form :: FA % TallyTotal )
        allocate ( Tally_F_D_Form :: FA % TallyChange )
        allocate ( FA % TallyBoundaryLocal  ( A % nBoundaries ) )
        allocate ( FA % TallyBoundaryGlobal ( A % nBoundaries ) )
        do iB = 1, A % nBoundaries 
          allocate &
            ( Tally_F_D_Form :: FA % TallyBoundaryLocal  ( iB ) % Element )
          allocate &
            ( Tally_F_D_Form :: FA % TallyBoundaryGlobal ( iB ) % Element )
        end do !-- iB
      case ( 'POLYTROPIC', 'NON_RELATIVISTIC', 'MEAN_HEAVY_NUCLEUS' )
        allocate ( Tally_F_P_Form :: FA % TallyInterior )
        allocate ( Tally_F_P_Form :: FA % TallyTotal )
        allocate ( Tally_F_P_Form :: FA % TallyChange )
        allocate ( FA % TallyBoundaryLocal  ( A % nBoundaries ) )
        allocate ( FA % TallyBoundaryGlobal ( A % nBoundaries ) )
        do iB = 1, A % nBoundaries
          allocate &
            ( Tally_F_P_Form :: FA % TallyBoundaryLocal  ( iB ) % Element )
          allocate &
            ( Tally_F_P_Form :: FA % TallyBoundaryGlobal ( iB ) % Element )
        end do !-- iB
      case default
        call Show ( 'FluidType not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- FluidType   
    end if !-- allocated TallyInterior
      
    select type ( TI => FA % TallyInterior )
    class is ( Tally_F_D_Form )
      call TI % Initialize &
             ( A, MassUnitOption = MassUnitOption, &
               EnergyUnitOption = EnergyUnitOption, &
               MomentumUnitOption = MomentumUnitOption, &
               AngularMomentumUnitOption = AngularMomentumUnitOption )
    end select !-- TI

    select type ( TT => FA % TallyTotal )
    class is ( Tally_F_D_Form )
      call TT % Initialize &
             ( A, MassUnitOption = MassUnitOption, &
               EnergyUnitOption = EnergyUnitOption, &
               MomentumUnitOption = MomentumUnitOption, &
               AngularMomentumUnitOption = AngularMomentumUnitOption )
    end select !-- TT

    select type ( TC => FA % TallyChange )
    class is ( Tally_F_D_Form )
      call TC % Initialize &
             ( A, MassUnitOption = MassUnitOption, &
               EnergyUnitOption = EnergyUnitOption, &
               MomentumUnitOption = MomentumUnitOption, &
               AngularMomentumUnitOption = AngularMomentumUnitOption )
    end select !-- TC

    do iB = 1, A % nBoundaries
      select type ( TB => FA % TallyBoundaryLocal ( iB ) % Element )
      class is ( Tally_F_D_Form )
        call TB % Initialize &
               ( A, MassUnitOption = MassUnitOption, &
                 EnergyUnitOption = EnergyUnitOption, &
                 MomentumUnitOption = MomentumUnitOption, &
                 AngularMomentumUnitOption = AngularMomentumUnitOption )
      end select !-- TB
      select type ( TB => FA % TallyBoundaryGlobal ( iB ) % Element )
      class is ( Tally_F_D_Form )
        call TB % Initialize &
               ( A, MassUnitOption = MassUnitOption, &
                 EnergyUnitOption = EnergyUnitOption, &
                 MomentumUnitOption = MomentumUnitOption, &
                 AngularMomentumUnitOption = AngularMomentumUnitOption )
      end select !-- TB
    end do !-- iB

    NameShort = 'Fluid'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call FA % InitializeTemplate_ASC_C ( A, NameShort, IgnorabilityOption )

    call Show ( FA % FluidType, 'FluidType', FA % IGNORABILITY )
    call Show ( FA % RiemannSolverType, 'RiemannSolverType', &
                FA % IGNORABILITY )
    call Show ( FA % UseLimiter, 'UseLimiter', FA % IGNORABILITY )
    call Show ( FA % LimiterParameter, 'LimiterParameter', FA % IGNORABILITY )

    allocate ( FA % Features_ASC )
    associate ( FFA => FA % Features_ASC )
    call FFA % Initialize &
           ( FA, FA % FluidType, &
             NameShortOption = trim ( NameShort ) // '_Features', &
             ShockThresholdOption = ShockThresholdOption, &
             IgnorabilityOption = IgnorabilityOption )
    select type ( FFC => FFA % Chart )
    class is ( FluidFeatures_CSL_Form )
      select type ( FC => FA % Chart )
      class is ( Fluid_CSL_Form )
        call FC % SetFeatures ( FFC )
      end select !-- FF
    end select !-- FFC
    end associate !-- FFA

  end subroutine Initialize


  function Fluid_D_CSL ( FA ) result ( F )

    class ( Fluid_ASC_Form ), intent ( in ) :: &
      FA
    class ( Fluid_D_Form ), pointer :: &
      F

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      F => FC % Fluid_D ( )
    class default
      call Show ( 'Fluid Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Fluid_D_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function Fluid_D_CSL


  function Fluid_P_P_CSL ( FA ) result ( F )

    class ( Fluid_ASC_Form ), intent ( in ) :: &
      FA
    class ( Fluid_P_P_Form ), pointer :: &
      F

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      F => FC % Fluid_P_P ( )
    class default
      call Show ( 'Fluid Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Fluid_P_P_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function Fluid_P_P_CSL


  function Fluid_P_NR_CSL ( FA ) result ( F )

    class ( Fluid_ASC_Form ), intent ( in ) :: &
      FA
    class ( Fluid_P_NR_Form ), pointer :: &
      F

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      F => FC % Fluid_P_NR ( )
    class default
      call Show ( 'Fluid Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Fluid_P_NR_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function Fluid_P_NR_CSL


  function Fluid_P_MHN_CSL ( FA ) result ( F )

    class ( Fluid_ASC_Form ), intent ( in ) :: &
      FA
    class ( Fluid_P_MHN_Form ), pointer :: &
      F

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      F => FC % Fluid_P_MHN ( )
    class default
      call Show ( 'Fluid Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Fluid_P_MHN_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function Fluid_P_MHN_CSL


  impure elemental subroutine Finalize ( FA )

    type ( Fluid_ASC_Form ), intent ( inout ) :: &
      FA

    if ( allocated ( FA % Features_ASC ) ) &
      deallocate ( FA % Features_ASC )

    call FA % FinalizeTemplate_ASC_C ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Fluid_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Fluid_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      call FC % Initialize &
             ( C, FA % NameShort, FA % FluidType, FA % RiemannSolverType, &
               FA % UseLimiter, FA % VelocityUnit, FA % MassDensityUnit, &
               FA % EnergyDensityUnit, FA % TemperatureUnit, &
               FA % LimiterParameter, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Fluid_ASC__Form
