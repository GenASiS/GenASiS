module RadiationMoments_CSL__Form

  !-- RadiationMoments_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use Interactions_Template
  use RadiationMoments_Form
  use Sources_RM__Form
  use Sources_RM_CSL__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: RadiationMoments_CSL_Form
    real ( KDR ) :: &
      LimiterParameter
    class ( StressEnergyUnitsForm ), pointer :: &
      Units => null ( )
    logical ( KDL ) :: &
      UseLimiter
    character ( LDF ) :: &
      RadiationType = '', &
      MomentsType = '', &
      ReconstructedType = '', &
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
      SetSources
    procedure, public, pass :: &
      SetInteractions
    final :: &
      Finalize
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      SetField
  end type RadiationMoments_CSL_Form

contains


  subroutine Initialize &
               ( RMC, C, NameShort, RadiationType, MomentsType, &
                 RiemannSolverType, ReconstructedType, UseLimiter, &
                 UsePinnedMemory, Units, LimiterParameter, nValues, &
                 IgnorabilityOption )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      RMC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      RadiationType, &
      MomentsType, &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter, &
      UsePinnedMemory
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    real ( KDR ), intent ( in ) :: &
      LimiterParameter
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call RMC % SetType ( )

    RMC % RadiationType        = RadiationType
    RMC % MomentsType          = MomentsType
    RMC % ReconstructedType    = ReconstructedType
    RMC % RiemannSolverType    = RiemannSolverType
    RMC % UseLimiter           = UseLimiter
    RMC % LimiterParameter     = LimiterParameter

    RMC % Units => Units

    call RMC % InitializeTemplate_CSL &
           ( C, NameShort, UsePinnedMemory, nValues, IgnorabilityOption )

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
    nullify ( RMC % Units )

    call RMC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetType ( RMC )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      RMC

    RMC % Type = 'a RadiationMoments_CSL'

  end subroutine SetType


  subroutine SetField ( FC )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % RadiationType ) )
    case ( 'GENERIC' )
      allocate ( RadiationMomentsForm :: FC % Field )
      select type ( RM => FC % Field )
      type is ( RadiationMomentsForm )
        call RM % Initialize &
               ( FC % RadiationType, FC % MomentsType, FC % RiemannSolverType, &
                 FC % ReconstructedType, FC % UseLimiter, FC % Units, &
                 FC % LimiterParameter, FC % nValues, &
                 PinnedOption = FC % UsePinnedMemory, &
                 NameOption = FC % NameShort )
        call RM % SetPrimitiveConserved ( )
        call RM % SetReconstructed ( )
        call RM % SetOutput ( FC % FieldOutput )
      end select !-- RM
    case default
      call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
      call Show ( FC % RadiationType, 'RadiationType', &
                  CONSOLE % ERROR )
      call Show ( 'RadiationMoments_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- RadiationMomentsType

  end subroutine SetField


end module RadiationMoments_CSL__Form
