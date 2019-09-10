module FluidFeatures_Template

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_TEMPLATE  = 3, &
      N_VECTORS_TEMPLATE = 0

  type, public, extends ( StorageForm ), abstract :: &
    FluidFeaturesTemplate
      integer ( KDI ) :: &
        IGNORABILITY       = 0, &
        N_FIELDS_TEMPLATE  = N_FIELDS_TEMPLATE, &
        N_VECTORS_TEMPLATE = N_VECTORS_TEMPLATE, &
        N_FIELDS          = 0, &
        N_VECTORS         = 0
      integer ( KDI ), dimension ( 3 ) :: &
        DIFFUSIVE_FLUX_I = 0
      character ( LDL ) :: &
        Type = ''
      class ( * ), pointer :: &
        Grid => null ( )
      class ( StorageForm ), pointer :: &
        Fluid => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure ( D ), public, pass, deferred :: &
      Detect
    procedure, public, pass :: &
      FinalizeTemplate
  end type FluidFeaturesTemplate

  abstract interface
    subroutine D ( FF )
      import FluidFeaturesTemplate
      class ( FluidFeaturesTemplate ), intent ( inout ) :: &
        FF
    end subroutine D
  end interface

    private :: &
      InitializeBasics

contains


  subroutine InitializeTemplate &
               ( FF, Fluid, Grid, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, PinnedOption, UnitOption, &
                 VectorIndicesOption )
                 
    class ( FluidFeaturesTemplate ), intent ( inout ) :: &
      FF
    class ( StorageForm ), intent ( in ), target :: &
      Fluid
    class ( * ), intent ( in ), target :: &
      Grid
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption, &
      PinnedOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), &
      optional :: &
        VectorIndicesOption

    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDF ) :: &
      Name 
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    logical ( KDL ) :: &
      Clear

    call InitializeBasics &
           ( FF, Variable, Name, VariableUnit, VariableOption, &
             NameOption, UnitOption )

    Clear = .true.
    if ( present ( ClearOption ) ) Clear = ClearOption

    call FF % StorageForm % Initialize &
           ( [ nValues, FF % N_FIELDS ], &
             VariableOption = Variable, VectorOption = VectorOption, &
             NameOption = Name, ClearOption = Clear, &
             PinnedOption = PinnedOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    FF % Grid  => Grid
    FF % Fluid => Fluid

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( FF )
  
    class ( FluidFeaturesTemplate ), intent ( inout ) :: &
      FF

    nullify ( FF % Fluid )
    nullify ( FF % Grid )

    call Show ( 'Finalizing ' // trim ( FF % Type ), FF % IGNORABILITY )
    call Show ( FF % Name, 'Name', FF % IGNORABILITY )
   
  end subroutine FinalizeTemplate


  subroutine InitializeBasics &
               ( FF, Variable, Name, VariableUnit, VariableOption, &
                 NameOption, VariableUnitOption )

    class ( FluidFeaturesTemplate ), intent ( inout ) :: &
      FF
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

    if ( FF % Type == '' ) &
      FF % Type = 'a FluidFeatures'

    Name = 'FluidFeatures'
    if ( present ( NameOption ) ) &
      Name = NameOption

    FF % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( FF % Type ), FF % IGNORABILITY )
    call Show ( Name, 'Name', FF % IGNORABILITY )

    !-- variable indices

    if ( FF % N_FIELDS == 0 ) &
      FF % N_FIELDS = FF % N_FIELDS_TEMPLATE

    FF % DIFFUSIVE_FLUX_I  =  [ 1, 2, 3 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( FF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : FF % N_FIELDS_TEMPLATE ) &
      = [ 'DiffusiveFlux_I_1', &
          'DiffusiveFlux_I_2', &
          'DiffusiveFlux_I_3' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( FF % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    if ( FF % N_VECTORS == 0 ) &
      FF % N_VECTORS = FF % N_VECTORS_TEMPLATE

  end subroutine InitializeBasics


end module FluidFeatures_Template
