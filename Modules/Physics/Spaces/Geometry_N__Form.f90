!-- Geometry_N represents Galilean geometry.

module Geometry_N__Form

  !-- Geometry_Newtonian__Form

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_NEWTONIAN  = 4, &
      N_VECTORS_NEWTONIAN = 1

  type, public, extends ( GeometryFlatForm ) :: Geometry_N_Form
    integer ( KDI ) :: &
      N_FIELDS_NEWTONIAN  = N_FIELDS_NEWTONIAN, &
      N_VECTORS_NEWTONIAN = N_VECTORS_NEWTONIAN
    integer ( KDI ) :: &
      POTENTIAL = 0
    integer ( KDI ), dimension ( 3 ) :: &
      GRAD_POTENTIAL_D = 0
    contains
      procedure, public, pass :: &
        InitializeAllocate_G
      final :: &
        Finalize
  end type Geometry_N_Form


contains


  subroutine InitializeAllocate_G &
               ( G, CoordinateSystem, CoordinateUnit, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

    class ( Geometry_N_Form ), intent ( inout ) :: &
      G
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      CoordinateUnit
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
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    integer ( KDI ) :: &
      iD
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

    call InitializeBasics &
           ( G, Variable, Vector, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, UnitOption, VectorIndicesOption )

    VariableUnit ( G % POTENTIAL )  &
      =  UNIT % SPEED_OF_LIGHT ** 2
    do iD = 1, 3
      VariableUnit ( G % GRAD_POTENTIAL_D ( iD ) ) &
        =  UNIT % SPEED_OF_LIGHT ** 2  /  CoordinateUnit ( iD )
    end do

    call G % GeometryFlatForm % Initialize &
           ( CoordinateSystem, CoordinateUnit, nValues, &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = VariableUnit, VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_G


  impure elemental subroutine Finalize ( G )

    type ( Geometry_N_Form ), intent ( inout ) :: &
      G

    !-- Empty routine to trigger finalization of parent type

  end subroutine Finalize


  subroutine InitializeBasics &
               ( G, Variable, Vector, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, VariableUnitOption, &
                 VectorIndicesOption )

    class ( Geometry_N_Form ), intent ( inout ) :: &
      G
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable, &
      Vector
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
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      iV, &  !-- iVector
      oF, &  !-- oField
      oV     !-- oVector

    if ( G % Type == '' ) &
      G % Type = 'Geometry_N'

    !-- variable indices

    oF = G % N_FIELDS_FLAT
    if ( G % N_FIELDS == 0 ) &
      G % N_FIELDS = oF + G % N_FIELDS_NEWTONIAN

    G % POTENTIAL        = oF + 1
    G % GRAD_POTENTIAL_D = oF + [ 2, 3, 4 ]

    !-- variable names

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( G % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + G % N_FIELDS_NEWTONIAN ) &
      = [ 'Potential         ', &
          'Grad_Potential_D_1', &
          'Grad_Potential_D_2', &
          'Grad_Potential_D_3' ]

    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( G % N_FIELDS ) )
    end if
    
    !-- vectors

    oV = G % N_VECTORS_FLAT
    if ( G % N_VECTORS == 0 ) &
      G % N_VECTORS = oV + G % N_VECTORS_NEWTONIAN

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( G % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( oV + 1 : oV + G % N_VECTORS_NEWTONIAN ) &
      = [ 'Grad_Potential_D' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = G % N_VECTORS_FLAT + G % N_VECTORS_NEWTONIAN + 1, &
              size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( G % N_VECTORS ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( G % GRAD_POTENTIAL_D )

  end subroutine InitializeBasics


end module Geometry_N__Form
