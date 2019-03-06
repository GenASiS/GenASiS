module Fluid_ASC__Form

  !-- Fluid_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergyBasics
  use Fluid_D__Form
  use Fluid_P_I__Form
  use Fluid_P_HN__Form
  use Tally_F_D__Form
  use Tally_F_P__Form
  use Tally_F_P_HN__Form
  use Sources_F_CSL__Form
  use Sources_F_ASC__Form
  use FluidFeatures_CSL__Form
  use FluidFeatures_ASC__Form
  use Fluid_CSL__Form
  
  implicit none
  private
  
  type, public, extends ( Current_ASC_Template ) :: Fluid_ASC_Form
    real ( KDR ) :: &
      LimiterParameter, &
      BaryonMassReference
    class ( StressEnergyUnitsForm ), pointer :: &
      Units => null ( )
    type ( MeasuredValueForm ) :: &
      BaryonMassUnit, &
      NumberDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U_Unit, &
      MomentumDensity_D_Unit
    logical ( KDL ) :: &
      UseLimiter, &
      UseEntropy
    character ( LDF ) :: &
      FluidType         = '', &
      RiemannSolverType = ''
    type ( Sources_F_ASC_Form ), allocatable :: &
      Sources_ASC
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
      Fluid_P_I_CSL
    generic, public :: &
      Fluid_P_I => Fluid_P_I_CSL
    procedure, private, pass :: &
      Fluid_P_HN_CSL
    generic, public :: &
      Fluid_P_HN => Fluid_P_HN_CSL
    procedure, public, pass :: &
      ComputeTally
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Fluid_ASC_Form

contains


  subroutine Initialize &
               ( FA, A, FluidType, Units, NameShortOption, &
                 RiemannSolverTypeOption, UseEntropyOption, UseLimiterOption, &
                 AllocateTallyOption, AllocateSourcesOption, &
                 BaryonMassReferenceOption, LimiterParameterOption, &
                 ShockThresholdOption, IgnorabilityOption )

    class ( Fluid_ASC_Form ), intent ( inout ) :: &
      FA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      FluidType
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    character ( * ), intent ( in ), optional :: &
      NameShortOption, &
      RiemannSolverTypeOption
    logical ( KDL ), intent ( in ), optional :: &
      UseEntropyOption, &
      UseLimiterOption, &
      AllocateTallyOption, &
      AllocateSourcesOption
    real ( KDR ), intent ( in ), optional :: &
      BaryonMassReferenceOption, &
      LimiterParameterOption, &
      ShockThresholdOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iB  !-- iBoundary
    character ( LDL ) :: &
      NameShort
    logical ( KDL ) :: &
      AllocateTally, &
      AllocateSources

    if ( FA % Type == '' ) &
      FA % Type = 'a Fluid_ASC'
    FA % FluidType = FluidType

    FA % Units => Units

    select case ( trim ( FA % FluidType ) )
    case ( 'DUST' )
      FA % RiemannSolverType = 'HLL'
    case default
      FA % RiemannSolverType = 'HLLC'
    end select
    if ( present ( RiemannSolverTypeOption ) ) &
      FA % RiemannSolverType = RiemannSolverTypeOption
    call PROGRAM_HEADER % GetParameter &
           ( FA % RiemannSolverType, 'RiemannSolverType' )

    select case ( trim ( FA % FluidType ) )
    case ( 'HEAVY_NUCLEUS' )
      FA % UseEntropy = .true.
    case default
      FA % UseEntropy = .false.
    end select
    if ( present ( UseEntropyOption ) ) &
      FA % UseEntropy = UseEntropyOption
    call PROGRAM_HEADER % GetParameter &
           ( FA % UseEntropy, 'UseEntropy' )

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

    FA % BaryonMassReference = 1.0_KDR
    if ( present ( BaryonMassReferenceOption ) ) &
      FA % BaryonMassReference = BaryonMassReferenceOption

    AllocateTally = .true.
    if ( present ( AllocateTallyOption ) ) &
      AllocateTally = AllocateTallyOption

    if ( AllocateTally ) then

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
        case ( 'IDEAL' )
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
        case ( 'HEAVY_NUCLEUS' )
          allocate ( Tally_F_P_HN_Form :: FA % TallyInterior )
          allocate ( Tally_F_P_HN_Form :: FA % TallyTotal )
          allocate ( Tally_F_P_HN_Form :: FA % TallyChange )
          allocate ( FA % TallyBoundaryLocal  ( A % nBoundaries ) )
          allocate ( FA % TallyBoundaryGlobal ( A % nBoundaries ) )
          do iB = 1, A % nBoundaries
            allocate &
              ( Tally_F_P_HN_Form :: FA % TallyBoundaryLocal ( iB ) &
                                        % Element )
            allocate &
              ( Tally_F_P_HN_Form :: FA % TallyBoundaryGlobal ( iB ) &
                                        % Element )
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
        call TI % Initialize ( A, Units )
      end select !-- TI

      select type ( TT => FA % TallyTotal )
      class is ( Tally_F_D_Form )
        call TT % Initialize ( A, Units )
      end select !-- TT

      select type ( TC => FA % TallyChange )
      class is ( Tally_F_D_Form )
        call TC % Initialize ( A, Units )
      end select !-- TC

      do iB = 1, A % nBoundaries
        select type ( TB => FA % TallyBoundaryLocal ( iB ) % Element )
        class is ( Tally_F_D_Form )
          call TB % Initialize ( A, Units )
        end select !-- TB
        select type ( TB => FA % TallyBoundaryGlobal ( iB ) % Element )
        class is ( Tally_F_D_Form )
          call TB % Initialize ( A, Units )
        end select !-- TB
      end do !-- iB

    end if !-- AllocateTally 

    NameShort = 'Fluid'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call FA % InitializeTemplate_ASC_C &
           ( A, NameShort, AllocateTallyOption = AllocateTallyOption, &
             IgnorabilityOption = IgnorabilityOption )

    call Show ( FA % FluidType, 'FluidType', FA % IGNORABILITY )
    call Show ( FA % RiemannSolverType, 'RiemannSolverType', &
                FA % IGNORABILITY )
    call Show ( FA % UseEntropy, 'UseEntropy', FA % IGNORABILITY )
    call Show ( FA % UseLimiter, 'UseLimiter', FA % IGNORABILITY )
    call Show ( FA % LimiterParameter, 'LimiterParameter', FA % IGNORABILITY )

    !-- Sources

    AllocateSources = .true.
    if ( present ( AllocateSourcesOption ) ) &
      AllocateSources = AllocateSourcesOption

    if ( AllocateSources ) then
      allocate ( FA % Sources_ASC )
      associate ( SFA => FA % Sources_ASC )
      call SFA % Initialize &
             ( FA, NameShortOption = trim ( NameShort ) // '_Sources', &
               TimeUnitOption = Units % Time, &
               IgnorabilityOption = IgnorabilityOption )
      select type ( SFC => SFA % Chart )
      class is ( Sources_F_CSL_Form )
        select type ( FC => FA % Chart )
        class is ( Fluid_CSL_Form )
          call FC % SetSources ( SFC )
        end select !-- FC
      end select !-- SFC
      end associate !-- SFA
    end if

    !-- Features

    if ( trim ( FA % FluidType ) /= 'DUST' ) then
      allocate ( FA % Features_ASC )
      associate ( FFA => FA % Features_ASC )
      call FFA % Initialize &
             ( FA, FA % FluidType, FA % RiemannSolverType, &
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
    end if !-- not DUST

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


  function Fluid_P_I_CSL ( FA ) result ( F )

    class ( Fluid_ASC_Form ), intent ( in ) :: &
      FA
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      F => FC % Fluid_P_I ( )
    class default
      call Show ( 'Fluid Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Fluid_P_I_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function Fluid_P_I_CSL


  function Fluid_P_HN_CSL ( FA ) result ( F )

    class ( Fluid_ASC_Form ), intent ( in ) :: &
      FA
    class ( Fluid_P_HN_Form ), pointer :: &
      F

    select type ( FC => FA % Chart )
    class is ( Fluid_CSL_Form )
      F => FC % Fluid_P_HN ( )
    class default
      call Show ( 'Fluid Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Fluid_P_HN_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FC

  end function Fluid_P_HN_CSL


  subroutine ComputeTally ( CA, ComputeChangeOption, IgnorabilityOption )
    
    class ( Fluid_ASC_Form ), intent ( inout ) :: &
      CA
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call CA % ComputeTallyTemplate ( ComputeChangeOption, IgnorabilityOption )

  end subroutine ComputeTally


  impure elemental subroutine Finalize ( FA )

    type ( Fluid_ASC_Form ), intent ( inout ) :: &
      FA

    if ( allocated ( FA % Features_ASC ) ) &
      deallocate ( FA % Features_ASC )
    if ( allocated ( FA % Sources_ASC ) ) &
      deallocate ( FA % Sources_ASC )

    nullify ( FA % Units )

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
               FA % UseEntropy, FA % UseLimiter, FA % Units, &
               FA % BaryonMassReference, FA % LimiterParameter, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Fluid_ASC__Form
