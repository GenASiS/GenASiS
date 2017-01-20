module RadiationMoments_CSL__Form

  !-- RadiationMoments_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use Interactions_Template
  use RadiationMoments_Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: RadiationMoments_CSL_Form
    type ( MeasuredValueForm ) :: &
      EnergyDensityUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      VelocityUnit
    character ( LDF ) :: &
      RadiationMomentsType = ''
    class ( Field_CSL_Template ), pointer :: &
      Interactions_CSL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      RadiationMoments
    procedure, public, pass :: &
      SetInteractions
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type RadiationMoments_CSL_Form

contains


  subroutine Initialize &
               ( RMC, C, RadiationMomentsType, VelocityUnit, &
                 EnergyDensityUnit, nValues, NameOutputOption )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      RMC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    character ( * ), intent ( in ) :: &
      RadiationMomentsType
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    if ( RMC % Type == '' ) &
      RMC % Type = 'a Fluid_CSL'
    RMC % RadiationMomentsType = RadiationMomentsType

    RMC % EnergyDensityUnit = EnergyDensityUnit
    RMC % VelocityUnit      = VelocityUnit

    call RMC % InitializeTemplate_CSL &
           ( C, nValues, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  function RadiationMoments ( RMC ) result ( RM )

    class ( RadiationMoments_CSL_Form ), intent ( in ), target :: &
      RMC
    class ( RadiationMomentsForm ), pointer :: &
      RM
      
    class ( VariableGroupForm ), pointer :: &
      Field

    RM => null ( )

    Field => RMC % Field
    select type ( Field )
    class is ( RadiationMomentsForm )
    RM => Field
    end select !-- Field

  end function RadiationMoments


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

    call RMC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC, NameOption )

    class ( RadiationMoments_CSL_Form ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in ), optional :: &
      NameOption

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % RadiationMomentsType ) )
    case ( 'GENERIC' )
      allocate ( RadiationMomentsForm :: FC % Field )
      select type ( RM => FC % Field )
      type is ( RadiationMomentsForm )
        call RM % Initialize &
               ( FC % VelocityUnit, FC % EnergyDensityUnit, FC % nValues, &
                 NameOption = NameOption )
        call RM % SetOutput ( FC % FieldOutput )
      end select !-- RM
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
