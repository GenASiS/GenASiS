module Interactions_Template

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_TEMPLATE = 3

  type, public, extends ( VariableGroupForm ), abstract :: InteractionsTemplate
    integer ( KDI ) :: &
      IGNORABILITY             = 0, &
      N_FIELDS                 = 0, &
      N_FIELDS_TEMPLATE        = N_FIELDS_TEMPLATE, &
      EMISSIVITY               = 0, &
      EFFECTIVE_OPACITY        = 0, &
      TRANSPORT_OPACITY        = 0!, &
!       EMISSIVITY_NUMBER        = 0, &
!       EFFECTIVE_OPACITY_NUMBER = 0
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
    procedure, private, pass ( I ) :: &
      ComputeEquilibrium_T
    generic, public :: &
      ComputeEquilibriumParameters => ComputeEquilibrium_T
  end type InteractionsTemplate

  abstract interface
    subroutine C ( I, Current )
      use Mathematics
      import InteractionsTemplate
      class ( InteractionsTemplate ), intent ( inout ) :: &
        I
      class ( CurrentTemplate ), intent ( in ) :: &
        Current
    end subroutine C
  end interface

    private :: &
      InitializeBasics, &
      SetUnits

contains


  subroutine InitializeTemplate &
               ( I, LengthUnit, EnergyDensityUnit, nValues, VariableOption, &
                 NameOption, ClearOption, UnitOption )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDF ) :: &
      Name 
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    logical ( KDL ) :: &
      Clear

    call InitializeBasics &
           ( I, Variable, Name, VariableUnit, VariableOption, NameOption, &
             UnitOption )

    call SetUnits ( VariableUnit, I, LengthUnit, EnergyDensityUnit )

    Clear = .true.
    if ( present ( ClearOption ) ) &
      Clear = ClearOption

    call I % VariableGroupForm % Initialize &
           ( [ nValues, I % N_FIELDS ], VariableOption = Variable, &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = VariableUnit )

  end subroutine InitializeTemplate


  subroutine SetOutput ( I, Output )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    call Output % Initialize &
           ( I, iaSelectedOption = [ I % EMISSIVITY, &
                                     I % EFFECTIVE_OPACITY, &
                                     I % TRANSPORT_OPACITY ] )

  end subroutine SetOutput


  impure elemental subroutine FinalizeTemplate ( I )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I

    call Show ( 'Finalizing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )
   
  end subroutine FinalizeTemplate


  subroutine ComputeEquilibrium_T ( T_EQ, I, C )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T_EQ
    class ( InteractionsTemplate ), intent ( in ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      C

    !-- Empty interface to be overridden later as needed

  end subroutine ComputeEquilibrium_T


  subroutine InitializeBasics &
               ( I, Variable, Name, VariableUnit, VariableOption, &
                 NameOption, VariableUnitOption )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable
    character ( LDF ), intent ( out ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption

    if ( I % Type == '' ) &
      I % Type = 'an Interactions'

    Name = 'Interactions'
    if ( present ( NameOption ) ) &
      Name = NameOption

    I % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( Name, 'Name', I % IGNORABILITY )

    !-- variable indices

    if ( I % N_FIELDS == 0 ) &
      I % N_FIELDS = I % N_FIELDS_TEMPLATE

    I % EMISSIVITY               =  1
    I % EFFECTIVE_OPACITY        =  2
    I % TRANSPORT_OPACITY        =  3
!     I % EMISSIVITY_NUMBER        =  4
!     I % EFFECTIVE_OPACITY_NUMBER =  5

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( I % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : I % N_FIELDS_TEMPLATE ) &
      = [ 'Emissivity      ', &
          'EffectiveOpacity', &
          'TransportOpacity' ]
!           'EmissivityNumber      ', &
!           'EffectiveOpacityNumber' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( I % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, I, LengthUnit, EnergyDensityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( InteractionsTemplate ), intent ( in ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit

    VariableUnit ( I % EMISSIVITY ) &
      =  EnergyDensityUnit / LengthUnit ** (-1)
    VariableUnit ( I % EFFECTIVE_OPACITY ) &
      =  LengthUnit ** (-1)
    VariableUnit ( I % TRANSPORT_OPACITY ) &
      =  LengthUnit ** (-1)

  end subroutine SetUnits


end module Interactions_Template
