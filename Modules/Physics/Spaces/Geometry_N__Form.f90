!-- Geometry_N represents Newtonian geometry.

module Geometry_N__Form

  !-- Geometry_Newtonian__Form

  use Basics
  use Mathematics
  use Geometry_G__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_NEWTONIAN  = 7, &
      N_VECTORS_NEWTONIAN = 2

  type, public, extends ( Geometry_G_Form ) :: Geometry_N_Form
    integer ( KDI ) :: &
      N_FIELDS_NEWTONIAN  = N_FIELDS_NEWTONIAN, &
      N_VECTORS_NEWTONIAN = N_VECTORS_NEWTONIAN
    integer ( KDI ) :: &
      POTENTIAL = 0
    integer ( KDI ), dimension ( 3 ) :: &
      POTENTIAL_GRADIENT_D = 0, &
      PRESSURE_GRADIENT_D  = 0
  contains
    procedure, public, pass :: &
      InitializeAllocate_G
    procedure, public, pass :: &
      SetOutput
    final :: &
      Finalize
  end type Geometry_N_Form

    private :: &
      InitializeBasics

contains


  subroutine InitializeAllocate_G &
               ( G, CoordinateSystem, CoordinateUnit, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 PinnedOption, UnitOption, VectorIndicesOption )

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
      ClearOption, & 
      PinnedOption
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

    if ( CoordinateUnit ( 1 )  /=  UNIT % IDENTITY ) then
      VariableUnit ( G % POTENTIAL )  &
        =  UNIT % SPEED_OF_LIGHT ** 2
      do iD = 1, 3
        VariableUnit ( G % POTENTIAL_GRADIENT_D ( iD ) ) &
          =  UNIT % SPEED_OF_LIGHT ** 2  /  CoordinateUnit ( iD )
        VariableUnit ( G % PRESSURE_GRADIENT_D ( iD ) ) &
          =  UNIT % ATOMIC_MASS_UNIT  *  UNIT % SPEED_OF_LIGHT ** 2  &
             *  UNIT % NUMBER_DENSITY_NUCLEAR /  CoordinateUnit ( iD )
      end do
    else
      !-- Dimensionless
    end if

    call G % Geometry_G_Form % Initialize &
           ( CoordinateSystem, CoordinateUnit, nValues, &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = NameOption, ClearOption = ClearOption, &
             PinnedOption = PinnedOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_G


  subroutine SetOutput ( G, Output )

    class ( Geometry_N_Form ), intent ( inout ) :: &
      G
    class ( StorageForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 2 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( G % POTENTIAL_GRADIENT_D )
    call VectorIndices ( 2 ) % Initialize ( G % PRESSURE_GRADIENT_D )
    call Output % Initialize &
           ( G, iaSelectedOption &
                  = [ G % CENTER_U, G % METRIC_DD_22, G % METRIC_DD_33, &
                      G % POTENTIAL, &
                      G % POTENTIAL_GRADIENT_D, &
                      G % PRESSURE_GRADIENT_D ], &
             VectorOption = [ 'PotentialGradient_D', &
                              'PressureGradient_D ' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  impure elemental subroutine Finalize ( G )

    type ( Geometry_N_Form ), intent ( inout ) :: &
      G

    !-- Trigger finalization of parent

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
      G % Type = 'a Geometry_N'

    !-- variable indices

    oF = G % N_FIELDS_FLAT
    if ( G % N_FIELDS == 0 ) &
      G % N_FIELDS = oF + G % N_FIELDS_NEWTONIAN

    G % POTENTIAL             =  oF + 1
    G % POTENTIAL_GRADIENT_D  =  oF + [ 2, 3, 4 ]
    G % PRESSURE_GRADIENT_D   =  oF + [ 5, 6, 7 ]

    !-- variable names

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( G % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + G % N_FIELDS_NEWTONIAN ) &
      = [ 'Potential            ', &
          'PotentialGradient_D_1', &
          'PotentialGradient_D_2', &
          'PotentialGradient_D_3', &
          'PressureGradient_D_1 ', &
          'PressureGradient_D_2 ', &
          'PressureGradient_D_3 ' ]

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
      = [ 'PotentialGradient_D', &
          'PressureGradient_D ' ]

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

    call VectorIndices ( oV + 1 ) % Initialize &
           ( G % POTENTIAL_GRADIENT_D )
    call VectorIndices ( oV + 2 ) % Initialize &
           ( G % PRESSURE_GRADIENT_D )

  end subroutine InitializeBasics


end module Geometry_N__Form
