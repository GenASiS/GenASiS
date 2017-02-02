module Interactions_Template

  use Basics
  use Mathematics

  implicit none
  private

  type, public, extends ( VariableGroupForm ), abstract :: InteractionsTemplate
    integer ( KDI ) :: &
      IGNORABILITY        = 0, &
      N_FIELDS            = 3, &
      EQUILIBRIUM_DENSITY = 1, &
      EFFECTIVE_OPACITY   = 2, &
      TRANSPORT_OPACITY   = 3
    character ( LDL ) :: &
      Type = ''
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure ( C ), public, pass, deferred :: &
      Compute
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass :: &
      FinalizeTemplate
  end type InteractionsTemplate

  abstract interface
    subroutine C ( I )
      import InteractionsTemplate
      class ( InteractionsTemplate ), intent ( inout ) :: &
        I
    end subroutine C
  end interface

contains


  subroutine InitializeTemplate &
               ( I, LengthUnit, EnergyDensityUnit, nValues, NameOption, &
                 ClearOption, UnitOption )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    character ( LDF ) :: &
      Name 
    logical ( KDL ) :: &
      Clear
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit

    if ( I % Type == '' ) &
      I % Type = 'an Interactions'

    Name = 'Interactions'
    if ( present ( NameOption ) ) &
      Name = NameOption

    I % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( Name, 'Name', I % IGNORABILITY )

    Clear = .true.
    if ( present ( ClearOption ) ) Clear = ClearOption

    allocate ( VariableUnit ( I % N_FIELDS ) )
    VariableUnit ( I % EQUILIBRIUM_DENSITY ) = EnergyDensityUnit
    VariableUnit ( I % EFFECTIVE_OPACITY )   = LengthUnit ** (-1)
    VariableUnit ( I % TRANSPORT_OPACITY )   = LengthUnit ** (-1)

    call I % VariableGroupForm % Initialize &
           ( [ nValues, I % N_FIELDS ], &
             VariableOption = [ 'EquilibriumDensity', &
                                'EffectiveOpacity  ', &
                                'TransportOpacity  ' ], &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = VariableUnit )

  end subroutine InitializeTemplate


  subroutine SetOutput ( I, Output )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    call Output % Initialize &
           ( I, iaSelectedOption = [ I % EQUILIBRIUM_DENSITY, &
                                     I % EFFECTIVE_OPACITY, &
                                     I % TRANSPORT_OPACITY ] )

  end subroutine SetOutput


  impure elemental subroutine FinalizeTemplate ( I )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I

    call Show ( 'Finalizing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )
   
  end subroutine FinalizeTemplate


end module Interactions_Template
