module FluidSources_Form

  use Basics
  use Mathematics
  use Fluid_D__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_FLUID  = 2, &
      N_VECTORS_FLUID = 0

  type, public, extends ( CurrentSourcesForm ) :: FluidSourcesForm
    integer ( KDI ) :: &
      N_FIELDS_FLUID  = N_FIELDS_FLUID, &
      N_VECTORS_FLUID = 0, &
      CURVILINEAR_S_1 = 0, &
      CURVILINEAR_S_2 = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    generic, public :: &
      Initialize => InitializeAllocate_FS
  end type FluidSourcesForm

    private :: &
      InitializeBasics, &
      SetUnits

contains


  subroutine InitializeAllocate_FS &
               ( FS, F, TimeUnit, VariableOption, VectorOption, NameOption, &
                 ClearOption, UnitOption, VectorIndicesOption )

    class ( FluidSourcesForm ), intent ( inout ) :: &
      FS
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

    FS % N_FIELDS_CONSERVED = F % N_CONSERVED

    call InitializeBasics &
           ( FS, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call SetUnits ( VariableUnit, FS, F, TimeUnit )

    call FS % CurrentSourcesForm % Initialize &
           ( F, TimeUnit, F % iaConserved, VariableOption = Variable, &
             VectorOption = Vector, NameOption = Name, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_FS


  subroutine InitializeBasics &
               ( FS, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( FluidSourcesForm ), intent ( inout ) :: &
      FS
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
      iV, &  !-- iVector
      oF, &  !-- oField
      oV     !-- oVector

    if ( FS % Type == '' ) &
      FS % Type = 'a FluidSources'

    Name = 'FluidSources'
    if ( present ( NameOption ) ) &
      Name = NameOption

    !-- variable indices

    oF = FS % N_FIELDS_CONSERVED
    if ( FS % N_FIELDS == 0 ) &
      FS % N_FIELDS = oF + FS % N_FIELDS_FLUID

    FS % CURVILINEAR_S_1  =  oF + 1
    FS % CURVILINEAR_S_2  =  oF + 2

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( FS % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + FS % N_FIELDS_FLUID ) &
      = [ 'Curvilinear_S_1', &
          'Curvilinear_S_2' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( FS % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = FS % N_VECTORS_CONSERVED
    if ( FS % N_VECTORS == 0 ) &
      FS % N_VECTORS = oV + FS % N_VECTORS_FLUID

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( FS % N_VECTORS ) )
      Vector = ''
    end if

  !   Vector ( oV + 1 : oV + FS % N_VECTORS_DUST ) &
  !     = [ 'Velocity                       ', &
  !         'MomentumDensity                ' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = oV + FS % N_VECTORS_FLUID + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( FS % N_VECTORS ) )
    end if

  !   call VectorIndices ( oV + 1 ) % Initialize ( FS % VELOCITY_U )
  !   call VectorIndices ( oV + 2 ) % Initialize ( FS % MOMENTUM_DENSITY_D )

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, FS, F, TimeUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( FluidSourcesForm ), intent ( in ) :: &
      FS
    class ( Fluid_D_Form ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit

    integer ( KDI ) :: &
      iMomentum_1, &
      iMomentum_2

    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )

    VariableUnit ( FS % CURVILINEAR_S_1 )  &
      =  F % Unit ( iMomentum_1 )  /  TimeUnit
    VariableUnit ( FS % CURVILINEAR_S_2 )  &
      =  F % Unit ( iMomentum_2 )  /  TimeUnit

  end subroutine SetUnits


end module FluidSources_Form
