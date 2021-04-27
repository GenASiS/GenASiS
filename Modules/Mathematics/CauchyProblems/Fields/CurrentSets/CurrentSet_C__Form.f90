module CurrentSet_C__Form

  !-- CurrentSet_Chart_Form

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_CS    = 6, &
      N_VECTORS_CS   = 2, &
      N_PRIMITIVE_CS = 0, &
      N_BALANCED_CS  = 0

  type, public, extends ( FieldSet_C_Form ) :: CurrentSet_C_Form
    !-- Fields and vectors
    integer ( KDI ) :: &
      N_FIELDS_CS    = N_FIELDS_CS, &
      N_VECTORS_CS   = N_VECTORS_CS
    integer ( KDI ) :: &
      FAST_EIGENSPEED_PLUS_U_1  = 0, &
      FAST_EIGENSPEED_PLUS_U_2  = 0, &
      FAST_EIGENSPEED_PLUS_U_3  = 0, &
      FAST_EIGENSPEED_MINUS_U_1 = 0, &
      FAST_EIGENSPEED_MINUS_U_2 = 0, &
      FAST_EIGENSPEED_MINUS_U_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      FAST_EIGENSPEED_PLUS_U, &
      FAST_EIGENSPEED_MINUS_U
    !-- Primtive and Balanced
    integer ( KDI ) :: &
      nPrimitive = 0, &
      nBalanced  = 0
    integer ( KDI ) :: &
      N_PRIMITIVE_CS = N_PRIMITIVE_CS, &
      N_BALANCED_CS  = N_BALANCED_CS
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaBalanced
  contains
    procedure, private, pass :: &
      InitializeAllocate_CS
    generic, public :: &
      Initialize => InitializeAllocate_CS
    procedure, public, pass ( CSC ) :: &
      SetStream
    procedure, private, pass :: &
      Show_FSC
    final :: &
      Finalize
  end type CurrentSet_C_Form


contains


  subroutine InitializeAllocate_CS &
               ( CSC, C, Velocity_U_Unit, FieldOption, VectorOption, &
                 NameOption, DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 iaPrimitiveOption, iaBalancedOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaPrimitiveOption, &
      iaBalancedOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iV, &  !-- iVector
      nFields, &
      nVectors
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector

    if ( CSC % Type  ==  '' ) &
      CSC % Type  =  'a CurrentSet_C' 
    
    Name  =  'CurrentSet'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    CSC % FAST_EIGENSPEED_PLUS_U_1   =  1
    CSC % FAST_EIGENSPEED_PLUS_U_2   =  2
    CSC % FAST_EIGENSPEED_PLUS_U_3   =  3
    CSC % FAST_EIGENSPEED_MINUS_U_1  =  4
    CSC % FAST_EIGENSPEED_MINUS_U_2  =  5
    CSC % FAST_EIGENSPEED_MINUS_U_3  =  6

!    !-- FIXME: GCC 10.1 is not initializing this correction 
!    !          in the type definition 
!    CSC % N_FIELDS_CS  =  N_FIELDS_CS

    nFields  =  CSC % N_FIELDS_CS
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

   CSC % FAST_EIGENSPEED_PLUS_U  &
     =  [ CSC % FAST_EIGENSPEED_PLUS_U_1, &
          CSC % FAST_EIGENSPEED_PLUS_U_2, &
          CSC % FAST_EIGENSPEED_PLUS_U_3 ]
   CSC % FAST_EIGENSPEED_MINUS_U  &
     =  [ CSC % FAST_EIGENSPEED_MINUS_U_1, &
          CSC % FAST_EIGENSPEED_MINUS_U_2, &
          CSC % FAST_EIGENSPEED_MINUS_U_3 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( 1 : CSC % N_FIELDS_CS ) &
      = [ 'FastEigenspeedPlus_U_1 ', &
          'FastEigenspeedPlus_U_2 ', &
          'FastEigenspeedPlus_U_3 ', &
          'FastEigenspeedMinus_U_1', &
          'FastEigenspeedMinus_U_2', &
          'FastEigenspeedMinus_U_3' ]
          
    !-- Units

    if ( present ( UnitOption ) ) then
      allocate ( Unit, source = UnitOption )
    else
      allocate ( Unit ( nFields ) )
    end if !-- UnitOption

    Unit ( CSC % FAST_EIGENSPEED_PLUS_U_1 : CSC % FAST_EIGENSPEED_PLUS_U_3 ) &
      =  Velocity_U_Unit
    Unit ( CSC % FAST_EIGENSPEED_MINUS_U_1 : CSC % FAST_EIGENSPEED_MINUS_U_3 ) &
      =  Velocity_U_Unit

    !-- Vector indices

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  CSC % N_VECTORS_CS + 1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  CSC % N_VECTORS_CS
      allocate ( VectorIndices ( nVectors ) )
    end if

    call VectorIndices ( 1 ) % Initialize ( CSC % FAST_EIGENSPEED_PLUS_U )
    call VectorIndices ( 2 ) % Initialize ( CSC % FAST_EIGENSPEED_MINUS_U )

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
    end if !-- FieldOption

    Vector ( 1 : CSC % N_VECTORS_CS ) &
      = [ 'FastEigenspeedPlus ', &
          'FastEigenspeedMinus' ]

    !-- FieldSet

    call CSC % FieldSet_C_Form % Initialize &
           ( C, &
             FieldOption = Field, &
             VectorOption = Vector, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             UnitOption = Unit, &
             VectorIndicesOption = VectorIndices, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    !-- Primitive fields

    if ( present ( iaPrimitiveOption ) ) then
      CSC % nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( CSC % iaPrimitive, source = iaPrimitiveOption )
    else
      CSC % nPrimitive  =  CSC % N_PRIMITIVE_CS
      allocate ( CSC % iaPrimitive ( CSC % nPrimitive ) )
    end if !-- iaPrimitiveOption

    !-- Balanced fields

    if ( present ( iaBalancedOption ) ) then
      CSC % nBalanced  =  size ( iaBalancedOption )
      allocate ( CSC % iaBalanced, source = iaBalancedOption )
    else
      CSC % nBalanced  =  CSC % N_BALANCED_CS
      allocate ( CSC % iaBalanced ( CSC % nBalanced ) )
    end if !-- iaBalancedOption

  end subroutine InitializeAllocate_CS


  subroutine SetStream ( SC, CSC )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    class ( CurrentSet_C_Form ), intent ( in ) :: &
      CSC

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected

    allocate ( iaSelected ( 0 ) )

    call SC % AddFieldSet &
           ( CSC, &
             iaSelectedOption  =  iaSelected )

  end subroutine SetStream


  subroutine Show_FSC ( FSC )

    class ( CurrentSet_C_Form ), intent ( in ) :: &
      FSC

    integer ( KDI ) :: &
      iF, &  !-- iField
      iS     !-- iSelected

    call FSC % FieldSet_C_Form % Show ( )

    call Show ( FSC % nPrimitive,   'nPrimitive', FSC % IGNORABILITY )
    call Show ( FSC % iaPrimitive, 'iaPrimitive', FSC % IGNORABILITY )

    call Show ( FSC % nBalanced,   'nBalanced', FSC % IGNORABILITY )
    call Show ( FSC % iaBalanced, 'iaBalanced', FSC % IGNORABILITY )

  end subroutine Show_FSC


  impure elemental subroutine Finalize ( CSC )

    type ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC

    if ( allocated ( CSC % iaBalanced ) ) &
      deallocate ( CSC % iaBalanced )
    if ( allocated ( CSC % iaPrimitive ) ) &
      deallocate ( CSC % iaPrimitive )

  end subroutine Finalize

  
end module CurrentSet_C__Form
