module FluidFeatures_Template

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_TEMPLATE  = 0, &
      N_VECTORS_TEMPLATE = 0

  type, public, extends ( VariableGroupForm ), abstract :: &
    FluidFeaturesTemplate
      integer ( KDI ) :: &
        IGNORABILITY       = 0, &
        N_FIELDS_TEMPLATE  = N_FIELDS_TEMPLATE, &
        N_VECTORS_TEMPLATE = N_VECTORS_TEMPLATE, &
        N_FIELDS          = 0, &
        N_VECTORS         = 0
      character ( LDL ) :: &
        Type = ''
      class ( * ), pointer :: &
        Grid => null ( )
      class ( VariableGroupForm ), pointer :: &
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
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( FluidFeaturesTemplate ), intent ( inout ) :: &
      FF
    class ( VariableGroupForm ), intent ( in ), target :: &
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
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), &
      optional :: &
        VectorIndicesOption

    character ( LDF ) :: &
      Name 
    logical ( KDL ) :: &
      Clear

    call InitializeBasics ( FF, Name, NameOption )

    Clear = .true.
    if ( present ( ClearOption ) ) Clear = ClearOption

    call FF % VariableGroupForm % Initialize &
           ( [ nValues, FF % N_FIELDS ], &
             VariableOption = VariableOption, VectorOption = VectorOption, &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = UnitOption, &
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


  subroutine InitializeBasics ( FF, Name, NameOption )

    class ( FluidFeaturesTemplate ), intent ( inout ) :: &
      FF
    character ( LDF ), intent ( out ) :: &
      Name
    character ( * ), intent ( in ), optional :: &
      NameOption

    if ( FF % Type == '' ) &
      FF % Type = 'a FluidFeatures'

    Name = 'FluidFeatures'
    if ( present ( NameOption ) ) &
      Name = NameOption

    FF % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( FF % Type ), FF % IGNORABILITY )
    call Show ( Name, 'Name', FF % IGNORABILITY )

  end subroutine InitializeBasics


end module FluidFeatures_Template
