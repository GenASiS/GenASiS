module Interactions_Template

  use Basics
  use Mathematics
  use Fluids

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_TEMPLATE = 9

  type, public, extends ( StorageForm ), abstract :: InteractionsTemplate
    integer ( KDI ) :: &
      IGNORABILITY      = 0, &
      N_FIELDS          = 0, &
      N_FIELDS_TEMPLATE = N_FIELDS_TEMPLATE, &
      EMISSIVITY_J      = 0, &
      EMISSIVITY_H      = 0, &
      EMISSIVITY_N      = 0, &
      OPACITY_J         = 0, &
      OPACITY_H         = 0, &
      OPACITY_N         = 0, &
      EQUILIBRIUM_J     = 0, &
      EQUILIBRIUM_N     = 0, &
      TIME_SCALE        = 0
    integer ( KDI ) :: &
      iBaseCell = 0
    real ( KDR ), dimension ( : ), pointer :: &
      Energy => null ( )
    character ( LDL ) :: &
      Type = '', &
      MomentsType = ''
    class ( Fluid_P_Template ), pointer :: &
      Fluid => null ( )
  contains
    procedure ( IAI ), private, pass, deferred :: &
      InitializeAllocate_I
    generic, public :: &
      Initialize => InitializeAllocate_I
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      SetOutput
    procedure ( C ), public, pass, deferred :: &
      Compute
    procedure, private, pass ( I ) :: &
      ComputeEquilibrium_T
    procedure, private, pass ( I ) :: &
      ComputeEquilibrium_T_Eta
    generic, public :: &
      ComputeEquilibriumParameters &
        => ComputeEquilibrium_T, ComputeEquilibrium_T_Eta
    procedure ( CTS ), public, pass, deferred :: &
      ComputeTimeScale
    procedure, public, pass :: &
      FinalizeTemplate
  end type InteractionsTemplate

  abstract interface

    subroutine IAI ( I, MomentsType, LengthUnit, EnergyDensityUnit, &
                     TemperatureUnit, nValues, VariableOption, NameOption, &
                     ClearOption, UnitOption )
      use Basics
      import InteractionsTemplate
      class ( InteractionsTemplate ), intent ( inout ) :: &
        I
      character ( * ), intent ( in ) :: &
        MomentsType
      type ( MeasuredValueForm ), intent ( in ) :: &
        LengthUnit, &
        EnergyDensityUnit, &
        TemperatureUnit
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
    end subroutine IAI

    subroutine C ( I, R )
      use Mathematics
      import InteractionsTemplate
      class ( InteractionsTemplate ), intent ( inout ) :: &
        I
      class ( CurrentTemplate ), intent ( in ) :: &
        R
    end subroutine C

    subroutine CTS ( I, R )
      use Basics
      use Mathematics
      import InteractionsTemplate
      class ( InteractionsTemplate ), intent ( inout ) :: &
        I
      class ( CurrentTemplate ), intent ( in ) :: &
        R
    end subroutine CTS

  end interface

    private :: &
      InitializeBasics, &
      SetUnits


contains


  subroutine InitializeTemplate &
               ( I, MomentsType, LengthUnit, EnergyDensityUnit, &
                 TemperatureUnit, nValues, VariableOption, NameOption, &
                 ClearOption, UnitOption )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I
    character ( * ), intent ( in ) :: &
      MomentsType
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
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

    I % MomentsType = MomentsType

    call SetUnits &
           ( Variableunit, I, LengthUnit, EnergyDensityUnit, TemperatureUnit )

    Clear = .true.
    if ( present ( ClearOption ) ) &
      Clear = ClearOption

    call I % StorageForm % Initialize &
           ( [ nValues, I % N_FIELDS ], VariableOption = Variable, &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = VariableUnit )

  end subroutine InitializeTemplate


  subroutine SetOutput ( I, Output )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I
    type ( StorageForm ), intent ( inout ) :: &
      Output

    call Output % Initialize &
           ( I, iaSelectedOption = [ I % EMISSIVITY_J, &
                                     I % EMISSIVITY_H, &
                                     I % EMISSIVITY_N, &
                                     I % OPACITY_J, &
                                     I % OPACITY_H, &
                                     I % OPACITY_N, &
                                     I % EQUILIBRIUM_J, &
                                     I % EQUILIBRIUM_N, &
                                     I % TIME_SCALE ] )

  end subroutine SetOutput


  subroutine ComputeEquilibrium_T ( T_EQ, I )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T_EQ
    class ( InteractionsTemplate ), intent ( in ) :: &
      I

    associate ( F => I % Fluid )
    call Copy ( F % Value ( :, F % TEMPERATURE ), T_EQ )
    end associate !-- F

  end subroutine ComputeEquilibrium_T


  subroutine ComputeEquilibrium_T_Eta ( T_EQ, Eta_EQ, I )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T_EQ, &
      Eta_EQ
    class ( InteractionsTemplate ), intent ( in ) :: &
      I

    !-- Empty interface to be overridden later as needed

  end subroutine ComputeEquilibrium_T_Eta


  impure elemental subroutine FinalizeTemplate ( I )

    class ( InteractionsTemplate ), intent ( inout ) :: &
      I

    nullify ( I % Fluid )

    call Show ( 'Finalizing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )
   
  end subroutine FinalizeTemplate


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

    I % EMISSIVITY_J   =  1
    I % EMISSIVITY_H   =  2
    I % EMISSIVITY_N   =  3
    I % OPACITY_J      =  4
    I % OPACITY_H      =  5
    I % OPACITY_N      =  6
    I % EQUILIBRIUM_J  =  7
    I % EQUILIBRIUM_N  =  8
    I % TIME_SCALE     =  9

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( I % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : I % N_FIELDS_TEMPLATE ) &
      = [ 'Emissivity_J ', &
          'Emissivity_H ', &
          'Emissivity_N ', &
          'Opacity_J    ', &
          'Opacity_H    ', &
          'Opacity_N    ', &
          'Equilibrium_J', &
          'Equilibrium_N', &
          'TimeScale    ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( I % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
  end subroutine InitializeBasics


  subroutine SetUnits &
               ( VariableUnit, I, LengthUnit, EnergyDensityUnit, &
                 TemperatureUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( InteractionsTemplate ), intent ( in ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit, &
      TemperatureUnit

    VariableUnit ( I % EMISSIVITY_J ) &
      =  EnergyDensityUnit  *  LengthUnit ** (-1)
    VariableUnit ( I % EMISSIVITY_H ) &
      =  EnergyDensityUnit  *  LengthUnit ** (-1)
    VariableUnit ( I % EMISSIVITY_N ) &
      =  EnergyDensityUnit  *  TemperatureUnit ** (-1)  *  LengthUnit ** (-1)
    VariableUnit ( I % OPACITY_J ) &
      =  LengthUnit ** (-1)
    VariableUnit ( I % OPACITY_H ) &
      =  LengthUnit ** (-1)
    VariableUnit ( I % OPACITY_N ) &
      =  LengthUnit ** (-1)
    VariableUnit ( I % EQUILIBRIUM_J ) &
      =  EnergyDensityUnit
    VariableUnit ( I % EQUILIBRIUM_N ) &
      =  EnergyDensityUnit  *  TemperatureUnit ** (-1)

  end subroutine SetUnits


end module Interactions_Template
