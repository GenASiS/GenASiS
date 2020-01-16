module Sources_F__Form

  use Basics
  use Mathematics
  use Fluid_D__Form
  use Fluid_P__Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_F  = 14, &
      N_VECTORS_F = 4

  type, public, extends ( Sources_C_Form ) :: Sources_F_Form
    integer ( KDI ) :: &
      N_FIELDS_F       = N_FIELDS_F, &
      N_VECTORS_F      = N_VECTORS_F, &
      GRAVITATIONAL_G  = 0, &
      RADIATION_G      = 0, &
      RADIATION_DS     = 0, &
      RADIATION_DE     = 0, &
      RADIATION_TIME   = 0
    integer ( KDI ), dimension ( 3 ) :: &
      CURVILINEAR_S_D   = 0, &
      GRAVITATIONAL_S_D = 0, &
      RADIATION_S_D     = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_SF
    generic, public :: &
      Initialize => InitializeAllocate_SF
    final :: &
      Finalize
  end type Sources_F_Form

    private :: &
      InitializeBasics, &
      SetUnits

contains


  subroutine InitializeAllocate_SF &
               ( SF, F, TimeUnit, VariableOption, VectorOption, NameOption, &
                 ClearOption, UnitOption, VectorIndicesOption )

    class ( Sources_F_Form ), intent ( inout ) :: &
      SF
    class ( Fluid_D_Form ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit
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

    SF % N_FIELDS_C = F % N_CONSERVED

    call InitializeBasics &
           ( SF, F, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call SetUnits ( VariableUnit, SF, F, TimeUnit )

    call SF % Sources_C_Form % Initialize &
           ( F, TimeUnit, F % iaConserved, VariableOption = Variable, &
             VectorOption = Vector, NameOption = Name, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_SF


  impure elemental subroutine Finalize ( SF )

    type ( Sources_F_Form ), intent ( inout ) :: &
      SF

    !-- Trigger finalization of parent

  end subroutine Finalize


  subroutine InitializeBasics &
               ( SF, F, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( Sources_F_Form ), intent ( inout ) :: &
      SF
    class ( Fluid_D_Form ), intent ( in ) :: &
      F
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
      iD, &  !-- iDimension
      iV, &  !-- iVector
      oF, &  !-- oField
      oV     !-- oVector
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum

    if ( SF % Type == '' ) &
      SF % Type = 'a Sources_F'

    Name = 'Sources_F'
    if ( present ( NameOption ) ) &
      Name = NameOption

    !-- variable indices

    oF = SF % N_FIELDS_C
    if ( SF % N_FIELDS == 0 ) &
      SF % N_FIELDS = oF + SF % N_FIELDS_F

    SF % GRAVITATIONAL_G    =  oF + 1
    SF % RADIATION_G        =  oF + 2
    SF % RADIATION_DS       =  oF + 3
    SF % RADIATION_DE       =  oF + 4
    SF % RADIATION_TIME     =  oF + 5
    SF % CURVILINEAR_S_D    =  oF + [  6,  7,  8 ]
    SF % GRAVITATIONAL_S_D  =  oF + [  9, 10, 11 ]
    SF % RADIATION_S_D      =  oF + [ 12, 13, 14 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( SF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + SF % N_FIELDS_F ) &
      = [ 'Gravitational_G    ', &
          'Radiation_G        ', &
          'Radiation_DS       ', &
          'Radiation_DE       ', &
          'RadiationTime      ', &
          'Curvilinear_S_D_1  ', &
          'Curvilinear_S_D_2  ', &
          'Curvilinear_S_D_3  ', &
          'Gravitational_S_D_1', &
          'Gravitational_S_D_2', &
          'Gravitational_S_D_3', &
          'Radiation_S_D_1    ', &
          'Radiation_S_D_2    ', &
          'Radiation_S_D_3    ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( SF % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = SF % N_VECTORS_C
    if ( SF % N_VECTORS == 0 ) &
      SF % N_VECTORS = oV + SF % N_VECTORS_F

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( SF % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( oV + 1 : oV + SF % N_VECTORS_F ) &
      = [ 'Div_F_S_D        ', &
          'Curvilinear_S_D  ', &
          'Gravitational_S_D', &
          'Radiation_S_D    ' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = oV + SF % N_VECTORS_F + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( SF % N_VECTORS ) )
    end if

    do iD = 1, 3
      call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( iD ), &
                    iMomentum ( iD ) )
    end do

    call VectorIndices ( oV + 1 ) % Initialize ( iMomentum )
    call VectorIndices ( oV + 2 ) % Initialize ( SF % CURVILINEAR_S_D )
    call VectorIndices ( oV + 3 ) % Initialize ( SF % GRAVITATIONAL_S_D )
    call VectorIndices ( oV + 4 ) % Initialize ( SF % RADIATION_S_D )

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, SF, F, TimeUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Sources_F_Form ), intent ( in ) :: &
      SF
    class ( Fluid_D_Form ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit

    integer ( KDI ) :: &
      iD

    select type ( F )
    class is ( Fluid_P_Template )
      VariableUnit ( SF % GRAVITATIONAL_G )  &
        =  F % Unit ( F % CONSERVED_ENERGY )  /  TimeUnit
      VariableUnit ( SF % RADIATION_G )  &
        =  F % Unit ( F % CONSERVED_ENERGY )  /  TimeUnit
      VariableUnit ( SF % RADIATION_DS )  &
        =  F % Unit ( F % CONSERVED_ENTROPY )  /  TimeUnit
      VariableUnit ( SF % RADIATION_DE )  &
        =  F % Unit ( F % CONSERVED_BARYON_DENSITY )  /  TimeUnit
      VariableUnit ( SF % RADIATION_TIME )  &
        =  TimeUnit
    end select !-- F

    do iD = 1, 3
      VariableUnit ( SF % CURVILINEAR_S_D ( iD ) )  &
        =  F % Unit ( F % MOMENTUM_DENSITY_D ( iD ) )  /  TimeUnit
      VariableUnit ( SF % GRAVITATIONAL_S_D ( iD ) )  &
        =  F % Unit ( F % MOMENTUM_DENSITY_D ( iD ) )  /  TimeUnit
      VariableUnit ( SF % RADIATION_S_D ( iD ) )  &
        =  F % Unit ( F % MOMENTUM_DENSITY_D ( iD ) )  /  TimeUnit
    end do

  end subroutine SetUnits


end module Sources_F__Form
