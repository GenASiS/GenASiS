module Gravitation_NH_C__Form

  !-- Gravitation_NewtonHeader_Chart_Form

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_N  = 4, &
      N_VECTORS_N = 1

  type, public, extends ( Geometry_F_C_Form ) :: Gravitation_NH_C_Form
    integer ( KDI ) :: &
      N_FIELDS_N = N_FIELDS_N, &
      N_VECTORS_N = N_VECTORS_N
    integer ( KDI ) :: &
      POTENTIAL = 0, &
      FORCE_D_1 = 0, &
      FORCE_D_2 = 0, &
      FORCE_D_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      FORCE_D
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    final :: &
      Finalize
  end type Gravitation_NH_C_Form


contains


  subroutine InitializeAllocate_FS &
               ( FSC, C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Gravitation_NH_C_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
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
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iV, &  !-- iVector
      oF, &  !-- oField
      oV, &  !-- oVector
      nFields, &
      nVectors
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a Gravitation_NH_C' 
    
    Name  =  'Gravitation'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  FSC % N_FIELDS_F

    FSC % POTENTIAL  =  oF + 1
    FSC % FORCE_D_1  =  oF + 2
    FSC % FORCE_D_2  =  oF + 3
    FSC % FORCE_D_3  =  oF + 4

    nFields  =  oF  +  FSC % N_FIELDS_N
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    FSC % FORCE_D  =  [ FSC % FORCE_D_1, &
                        FSC % FORCE_D_2, &
                        FSC % FORCE_D_3 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + FSC % N_FIELDS_N ) &
      = [ 'Potential', &
          'Force_D_1', &
          'Force_D_2', &
          'Force_D_3' ]

    !-- Units

    !-- Vector indices

    oV  =  FSC % N_VECTORS_F

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  oV  +  FSC % N_VECTORS_N  +  1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  oV  +  FSC % N_VECTORS_N
      allocate ( VectorIndices ( nVectors ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( FSC % FORCE_D )

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
    end if !-- FieldOption

    Vector ( oV  +  1 : oV  +  FSC % N_VECTORS_N ) &
      = [ 'Force_D' ]

    !-- Geometry_F

    call FSC % Geometry_F_C_Form % Initialize &
           ( C, &
             FieldOption = Field, &
             VectorOption = Vector, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndices, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_FS


  impure elemental subroutine Finalize ( GC )

    type ( Gravitation_NH_C_Form ), intent ( inout ) :: &
      GC

  end subroutine Finalize


end module Gravitation_NH_C__Form
