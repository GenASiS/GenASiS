!-- CurrentSources contains sources for Currents.

module CurrentSources_Form

  use Basics

  implicit none
  private

  type, public, extends ( VariableGroupForm ) :: CurrentSourcesForm
    integer ( KDI ) :: &
      IGNORABILITY        = 0, &
      N_FIELDS_CONSERVED  = 0, &
      N_VECTORS_CONSERVED = 0, &
      N_FIELDS            = 0, &
      N_VECTORS           = 0
    character ( LDL ) :: &
      Type = ''
  contains
    procedure, private, pass :: &
      InitializeConserved
    generic, public :: &
      Initialize => InitializeConserved
    final :: &
      Finalize
    procedure, public, pass :: &
      SetOutput
  end type CurrentSourcesForm

    private :: &
      InitializeBasics, &
      SetUnits

contains


  subroutine InitializeConserved &
               ( CS, Current, TimeUnit, iaConserved, VariableOption, &
                 VectorOption, NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( CurrentSourcesForm ), intent ( inout ) :: &
      CS
    class ( VariableGroupForm ), intent ( in ) :: &
      Current
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaConserved
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

    CS % N_FIELDS_CONSERVED = size ( iaConserved )

    call InitializeBasics &
           ( CS, Current % Variable, iaConserved, Variable, Vector, Name, &
             VariableUnit, VectorIndices, VariableOption, VectorOption, &
             NameOption, UnitOption, VectorIndicesOption )

    call SetUnits &
           ( VariableUnit, CS, Current % Unit, TimeUnit, iaConserved )

    Clear = .true.
    if ( present ( ClearOption ) ) Clear = ClearOption

    call CS % VariableGroupForm % Initialize &
           ( [ Current % nValues, CS % N_FIELDS ], &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeConserved


  impure elemental subroutine Finalize ( CS )

    type ( CurrentSourcesForm ), intent ( inout ) :: &
      CS

    call Show ( 'Finalizing ' // trim ( CS % Type ), CS % IGNORABILITY )
    call Show ( CS % Name, 'Name', CS % IGNORABILITY )
   
  end subroutine Finalize


  subroutine SetOutput ( CS, Output )

    class ( CurrentSourcesForm ), intent ( inout ) :: &
      CS
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    call Output % Initialize ( CS )

  end subroutine SetOutput


  subroutine InitializeBasics &
               ( CS, VariableCurrent, iaConserved, Variable, Vector, Name, &
                 VariableUnit, VectorIndices, VariableOption, VectorOption, &
                 NameOption, VariableUnitOption, VectorIndicesOption )

    class ( CurrentSourcesForm ), intent ( inout ) :: &
      CS
    character ( LDL ), dimension ( : ), intent ( in ) :: &
      VariableCurrent
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaConserved
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
      iC, &  !-- iConserved
      iV     !-- iVector

    if ( CS % Type == '' ) &
      CS % Type = 'a CurrentSources'

    Name = 'CurrentSources'
    if ( present ( NameOption ) ) &
      Name = NameOption

    CS % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( CS % Type ), CS % IGNORABILITY )
    call Show ( Name, 'Name', CS % IGNORABILITY )

    !-- variable indices

    if ( CS % N_FIELDS == 0 ) &
      CS % N_FIELDS = CS % N_FIELDS_CONSERVED

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( CS % N_FIELDS ) )
      Variable = ''
    end if

    do iC = 1, CS % N_FIELDS_CONSERVED
      Variable ( iC ) = 'Div_Flux_' // VariableCurrent ( iaConserved ( iC ) )
    end do !-- iC
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( CS % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    if ( CS % N_VECTORS == 0 ) &
      CS % N_VECTORS = CS % N_VECTORS_CONSERVED

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( CS % N_VECTORS ) )
      Vector = ''
    end if

!    Vector ( 1 : C % N_VECTORS_CONSERVED ) &
!      = [ '' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = CS % N_VECTORS_CONSERVED + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( CS % N_VECTORS ) )
    end if

!    call VectorIndices ( 1 ) % Initialize ( CS %  )

  end subroutine InitializeBasics


  subroutine SetUnits &
               ( VariableUnit, CS, VariableUnitCurrent, TimeUnit, &
                 iaConserved )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( CurrentSourcesForm ), intent ( in ) :: &
      CS
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ) :: &
      VariableUnitCurrent
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaConserved

    integer ( KDI ) :: &
      iC  !-- iConserved

    do iC = 1, CS % N_FIELDS_CONSERVED
      VariableUnit ( iC )  &
        =  VariableUnitCurrent ( iaConserved ( iC ) )  /  TimeUnit
    end do

  end subroutine SetUnits


end module CurrentSources_Form
