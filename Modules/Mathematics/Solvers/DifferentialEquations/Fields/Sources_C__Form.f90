!-- Sources_C contains sources for Currents.

module Sources_C__Form

  !-- Sources_Current_Form

  use Basics

  implicit none
  private

  type, public, extends ( VariableGroupForm ) :: Sources_C_Form
    integer ( KDI ) :: &
      IGNORABILITY        = 0, &
      N_FIELDS_C  = 0, &
      N_VECTORS_C = 0, &
      N_FIELDS            = 0, &
      N_VECTORS           = 0
    character ( LDL ) :: &
      Type = ''
  contains
    procedure, private, pass :: &
      InitializeAllocate_SC
    generic, public :: &
      Initialize => InitializeAllocate_SC
    final :: &
      Finalize
    procedure, public, pass :: &
      SetOutput
  end type Sources_C_Form

    private :: &
      InitializeBasics, &
      SetUnits

contains


  subroutine InitializeAllocate_SC &
               ( SC, Current, TimeUnit, iaConserved, VariableOption, &
                 VectorOption, NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( Sources_C_Form ), intent ( inout ) :: &
      SC
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

    SC % N_FIELDS_C = size ( iaConserved )

    call InitializeBasics &
           ( SC, Current % Variable, iaConserved, Variable, Vector, Name, &
             VariableUnit, VectorIndices, VariableOption, VectorOption, &
             NameOption, UnitOption, VectorIndicesOption )

    call SetUnits &
           ( VariableUnit, SC, Current % Unit, TimeUnit, iaConserved )

    Clear = .true.
    if ( present ( ClearOption ) ) Clear = ClearOption

    call SC % VariableGroupForm % Initialize &
           ( [ Current % nValues, SC % N_FIELDS ], &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_SC


  impure elemental subroutine Finalize ( SC )

    type ( Sources_C_Form ), intent ( inout ) :: &
      SC

    call Show ( 'Finalizing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )
   
  end subroutine Finalize


  subroutine SetOutput ( SC, Output )

    class ( Sources_C_Form ), intent ( inout ) :: &
      SC
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    call Output % Initialize ( SC )

  end subroutine SetOutput


  subroutine InitializeBasics &
               ( SC, VariableCurrent, iaConserved, Variable, Vector, Name, &
                 VariableUnit, VectorIndices, VariableOption, VectorOption, &
                 NameOption, VariableUnitOption, VectorIndicesOption )

    class ( Sources_C_Form ), intent ( inout ) :: &
      SC
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

    if ( SC % Type == '' ) &
      SC % Type = 'a Sources_C'

    Name = 'Sources_C_'
    if ( present ( NameOption ) ) &
      Name = NameOption

    SC % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( Name, 'Name', SC % IGNORABILITY )

    !-- variable indices

    if ( SC % N_FIELDS == 0 ) &
      SC % N_FIELDS = SC % N_FIELDS_C

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( SC % N_FIELDS ) )
      Variable = ''
    end if

    do iC = 1, SC % N_FIELDS_C
      Variable ( iC ) = 'Div_Flux_' // VariableCurrent ( iaConserved ( iC ) )
    end do !-- iC
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( SC % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    if ( SC % N_VECTORS == 0 ) &
      SC % N_VECTORS = SC % N_VECTORS_C

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( SC % N_VECTORS ) )
      Vector = ''
    end if

!    Vector ( 1 : C % N_VECTORS_C ) &
!      = [ '' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = SC % N_VECTORS_C + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( SC % N_VECTORS ) )
    end if

!    call VectorIndices ( 1 ) % Initialize ( S_C %  )

  end subroutine InitializeBasics


  subroutine SetUnits &
               ( VariableUnit, SC, VariableUnitCurrent, TimeUnit, &
                 iaConserved )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Sources_C_Form ), intent ( in ) :: &
      SC
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ) :: &
      VariableUnitCurrent
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaConserved

    integer ( KDI ) :: &
      iC  !-- iConserved

    do iC = 1, SC % N_FIELDS_C
      VariableUnit ( iC )  &
        =  VariableUnitCurrent ( iaConserved ( iC ) )  /  TimeUnit
    end do

  end subroutine SetUnits


end module Sources_C__Form
