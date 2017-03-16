module FluidFeatures_P__Form

  !-- FluidFeatures_Perfect__Form

  use Basics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_PERFECT  = 6, &
      N_VECTORS_PERFECT = 2

  type, public, extends ( VariableGroupForm ) :: FluidFeatures_P_Form
    integer ( KDI ) :: &
      IGNORABILITY      = 0, &
      N_FIELDS_PERFECT  = N_FIELDS_PERFECT, &
      N_VECTORS_PERFECT = N_VECTORS_PERFECT, &
      N_FIELDS          = 0, &
      N_VECTORS         = 0
    integer ( KDI ), dimension ( 3 ) :: &
      SHOCK_I   = 0, &
      CONTACT_I = 0
    real ( KDR ) :: &
      ShockThreshold
    character ( LDL ) :: &
      Type = ''
    class ( * ), pointer :: &
      Grid => null ( )
  contains
    procedure, public, pass :: &
      InitializeAllocate_P
    generic, public :: &
      Initialize => InitializeAllocate_P
    final :: &
      Finalize
    procedure, public, pass :: &
      SetOutput
  end type FluidFeatures_P_Form

    private :: &
      InitializeBasics

contains


  subroutine InitializeAllocate_P &
               ( FF, Grid, ShockThreshold, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    class ( * ), intent ( in ), target :: &
      Grid
    real ( KDR ), intent ( in ) :: &
      ShockThreshold
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

    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDF ) :: &
      Name 
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector
    logical ( KDL ) :: &
      Clear

    Clear = .true.
    if ( present ( ClearOption ) ) Clear = ClearOption

    call FF % VariableGroupForm % Initialize &
           ( [ nValues, FF % N_FIELDS ], &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

    FF % ShockThreshold = ShockThreshold
    FF % Grid => Grid

  end subroutine InitializeAllocate_P


  impure elemental subroutine Finalize ( FF )

    type ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF

    call Show ( 'Finalizing ' // trim ( FF % Type ), FF % IGNORABILITY )
    call Show ( FF % Name, 'Name', FF % IGNORABILITY )
   
  end subroutine Finalize


  subroutine SetOutput ( FF, Output )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    call Output % Initialize ( FF )

  end subroutine SetOutput


  subroutine InitializeBasics &
               ( FF, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable, &
      Vector
    character ( LDF ), intent ( out ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    !-- FIXME: intent(out) here caused ICE with Intel Compiler 15
    !          Temporarily set to intent(inout)
    !type ( Integer_1D_Form ), dimension ( : ), allocatable, &
    !  intent ( out ) :: &
    type ( Integer_1D_Form ), dimension ( : ), allocatable, &
      intent ( inout ) :: &
        VectorIndices
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      iV  !-- iVector

    if ( FF % Type == '' ) &
      FF % Type = 'a FluidFeatures_P'

    Name = 'FluidFeatures'
    if ( present ( NameOption ) ) &
      Name = NameOption

    FF % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( FF % Type ), FF % IGNORABILITY )
    call Show ( Name, 'Name', FF % IGNORABILITY )

    !-- variable indices

    if ( FF % N_FIELDS == 0 ) &
      FF % N_FIELDS = FF % N_FIELDS_PERFECT

    FF % SHOCK_I   = [ 1, 2, 3 ]
    FF % CONTACT_I = [ 4, 5, 6 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( FF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : FF % N_FIELDS_PERFECT ) &
      = [ 'Shock_I_1  ', &
          'Shock_I_2  ', &
          'Shock_I_3  ', &
          'Contact_I_1', &
          'Contact_I_2', &
          'Contact_I_3' ]
          
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
      FF % N_VECTORS = FF % N_VECTORS_PERFECT

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( FF % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( 1 : FF % N_VECTORS_PERFECT ) &
      = [ 'Shock_I  ', &
          'Contact_I' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = FF % N_VECTORS_PERFECT + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( FF % N_VECTORS ) )
    end if

    call VectorIndices ( 1 ) % Initialize ( FF % SHOCK_I )
    call VectorIndices ( 2 ) % Initialize ( FF % CONTACT_I )

  end subroutine InitializeBasics


end module FluidFeatures_P__Form
