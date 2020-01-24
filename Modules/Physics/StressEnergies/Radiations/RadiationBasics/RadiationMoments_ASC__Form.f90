module RadiationMoments_ASC__Form

  !-- RadiationMoments_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationMoments_Form
  use Sources_RM_CSL__Form
  use Sources_RM_ASC__Form
  use RadiationMoments_CSL__Form

  implicit none
  private
  
  type, public, extends ( Current_ASC_Template ) :: RadiationMoments_ASC_Form
    real ( KDR ) :: &
      LimiterParameter
    class ( StressEnergyUnitsForm ), pointer :: &
      Units => null ( )
    logical ( KDL ) :: &
      UseLimiter, &
      SuppressWrite
    character ( LDF ) :: &
      RadiationType = '', &
      MomentsType = '', &
      ReconstructedType = '', &
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
    procedure, public, pass :: &
      SetInteractions
    final :: &
      Finalize
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      SetField
    procedure, private, pass :: &
      AllocateField
  end type RadiationMoments_ASC_Form

contains


  subroutine Initialize &
               ( RMA, A, RadiationType, MomentsType, Units, NameShortOption, &
                 RiemannSolverTypeOption, ReconstructedTypeOption, &
                 UseLimiterOption, AllocateTallyOption, AllocateSourcesOption, &
                 SuppressWriteOption, SuppressWriteSourcesOption, &
                 LimiterParameterOption, IgnorabilityOption )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      RMA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      RadiationType, &
      MomentsType
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    character ( * ), intent ( in ), optional :: &
      NameShortOption, &
      RiemannSolverTypeOption, &
      ReconstructedTypeOption
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterOption, &
      AllocateTallyOption, &
      AllocateSourcesOption, &
      SuppressWriteOption, &
      SuppressWriteSourcesOption
    real ( KDR ), intent ( in ), optional :: &
      LimiterParameterOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

!    integer ( KDI ) :: &
!      iB  !-- iBoundary
    character ( LDL ) :: &
      NameShort
    character ( LDL ), dimension ( 5 ) :: &
      TallyVariable
    logical ( KDL ) :: &
      AllocateSources
    type ( MeasuredValueForm ), dimension ( 5 ) :: &
      TallyUnit

    call RMA % SetType ( )

    RMA % RadiationType = RadiationType
    RMA % MomentsType   = MomentsType

    RMA % Units => Units

    RMA % RiemannSolverType = 'HLL'
    if ( present ( RiemannSolverTypeOption ) ) &
      RMA % RiemannSolverType = RiemannSolverTypeOption
    call PROGRAM_HEADER % GetParameter &
           ( RMA % RiemannSolverType, 'RiemannSolverType' )

    RMA % ReconstructedType = 'PRIMITIVE'
    if ( present ( ReconstructedTypeOption ) ) &
      RMA % ReconstructedType = ReconstructedTypeOption
    call PROGRAM_HEADER % GetParameter &
           ( RMA % ReconstructedType, 'ReconstructedType' )

    RMA % UseLimiter = .false.
    if ( present ( UseLimiterOption ) ) &
      RMA % UseLimiter = UseLimiterOption
    call PROGRAM_HEADER % GetParameter &
           ( RMA % UseLimiter, 'UseLimiter' )

    RMA % LimiterParameter = 1.4_KDR
    if ( present ( LimiterParameterOption ) ) &
      RMA % LimiterParameter = LimiterParameterOption
    call PROGRAM_HEADER % GetParameter &
           ( RMA % LimiterParameter, 'LimiterParameter' )

    NameShort = 'Radiation'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    RMA % SuppressWrite = .false.
    if ( present ( SuppressWriteOption ) ) &
      RMA % SuppressWrite = SuppressWriteOption

    TallyVariable ( 1 )  =  'Energy'
    TallyVariable ( 2 )  =  'Momentum_1'
    TallyVariable ( 3 )  =  'Momentum_2'
    TallyVariable ( 4 )  =  'Momentum_3'
    TallyVariable ( 5 )  =  'Number'

    TallyUnit ( 1 )  =  Units % Energy
    TallyUnit ( 2 )  =  Units % Momentum
    TallyUnit ( 3 )  =  Units % Momentum
    TallyUnit ( 4 )  =  Units % Momentum
    TallyUnit ( 5 )  =  Units % Number 

    call RMA % InitializeTemplate_ASC_C &
           ( A, NameShort, TallyVariableOption = TallyVariable, &
             AllocateTallyOption = AllocateTallyOption, &
             TallyUnitOption = TallyUnit, &
             IgnorabilityOption = IgnorabilityOption )

    call Show ( RMA % RadiationType, 'RadiationType', &
                RMA % IGNORABILITY )
    call Show ( RMA % MomentsType, 'MomentsType', &
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
               TimeUnitOption = Units % Time, &
               EnergyUnitOption = Units % Coordinate_MS ( 1 ), &
               IgnorabilityOption = IgnorabilityOption, &
               SuppressWriteOption = SuppressWriteSourcesOption )
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

    nullify ( RMA % Units )

    call RMA % FinalizeTemplate_ASC_C ( )

  end subroutine Finalize


  subroutine SetType ( RMA )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      RMA

    RMA % Type = 'a RadiationMoments_ASC'

  end subroutine SetType


  subroutine SetField ( FA )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    call FA % AllocateField ( )

    select type ( FC => FA % Chart )
    class is ( RadiationMoments_CSL_Form )
      call FC % Initialize &
             ( C, FA % NameShort, FA % RadiationType, FA % MomentsType, &
               FA % RiemannSolverType, FA % ReconstructedType, &
               FA % UseLimiter, FA % Units, FA % LimiterParameter, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    if ( .not. FA % SuppressWrite ) &
      call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


  subroutine AllocateField ( RMA )

    class ( RadiationMoments_ASC_Form ), intent ( inout ) :: &
      RMA

    allocate ( RadiationMoments_CSL_Form :: RMA % Chart )

  end subroutine AllocateField


end module RadiationMoments_ASC__Form
