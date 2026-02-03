module Features_CS__Form

  !-- Features_CurrentSet__Form

  use Basics
  use FieldSets
  use Geometries

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_CS  = 3, &
      N_VECTORS_CS = 0

  type, public, extends ( FieldSet_BM_Form ) :: Features_CS_Form
    integer ( KDI ) :: &
      N_FIELDS_CS    = N_FIELDS_CS, &
      N_VECTORS_CS   = N_VECTORS_CS
    integer ( KDI ), dimension ( 3 ) :: &
      DIFFUSIVE_FLUX_I = 0
    class ( FieldSet_BM_Form ), pointer :: &
      CurrentSet => null ( )
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
  contains
    procedure, public, pass :: &
      InitializeAllocate_CS
    generic, public :: &
      Initialize => InitializeAllocate_CS
    procedure, public, pass :: &
      Detect
    final :: &
      Finalize
  end type Features_CS_Form

contains


  subroutine InitializeAllocate_CS &
               ( F, G, CS, FieldOption, VectorOption, NameOption, UnitOption, &
                 VectorIndicesOption, nFieldsOption, IgnorabilityOption )

    class ( Features_CS_Form ), intent ( inout ) :: &
      F
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      CS
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iV, &  !-- iVector
      nFields, &
      nVectors
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      Unit
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector
    character ( LDL ) :: &
      Name

    if ( F % Type  ==  '' ) &
      F % Type  =  'a Features_CS' 
    
    Name  =  'Features' // trim ( CS % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    F % Geometry    =>  G
    F % CurrentSet  =>  CS

    !-- Field indices

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else
      nFields  =  F % N_FIELDS_CS
    end if

    F % DIFFUSIVE_FLUX_I  =  [ 1, 2, 3 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( 1 )  =  'DiffusiveFlux_I_1'
    Field ( 2 )  =  'DiffusiveFlux_I_2'
    Field ( 3 )  =  'DiffusiveFlux_I_3'

    !-- Units

    if ( present ( UnitOption ) ) then
      allocate ( Unit, source = UnitOption )
    else
      allocate ( Unit ( nFields, G % Atlas % nCharts ) )
    end if !-- UnitOption

    !-- Vector indices

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  F % N_VECTORS_CS + 1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  F % N_VECTORS_CS
      allocate ( VectorIndices ( nVectors ) )
    end if

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
    end if !-- FieldOption

    !-- FieldSet

    call F % FieldSet_BM_Form % Initialize &
           ( G % Atlas, &
             FieldOption = Field, &
             VectorOption = Vector, &
             NameOption = Name, &
             DeviceMemoryOption = G % DeviceMemory, &
             PinnedMemoryOption = G % PinnedMemory, &
             DevicesCommunicateOption = G % DevicesCommunicate, &
             UnitOption = Unit, &
             VectorIndicesOption = VectorIndices, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_CS


  subroutine Detect ( F )

    class ( Features_CS_Form ), intent ( inout ) :: &
      F

    call Show ( 'Detect should be overridden', CONSOLE % WARNING )
    call Show ( 'Features_CS_Form', 'module', CONSOLE % WARNING )
    call Show ( 'Detect', 'subroutine', CONSOLE % WARNING )

  end subroutine Detect


  impure elemental subroutine Finalize ( F )

    type ( Features_CS_Form ), intent ( inout ) :: &
      F

    nullify ( F % CurrentSet )
    nullify ( F % Geometry )

  end subroutine Finalize


end module Features_CS__Form
