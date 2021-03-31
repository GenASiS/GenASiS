module FieldSet_MH__Form

  !-- FieldSet_ManifoldHeader__Form

  use Basics
  use Manifold_H__Form

  implicit none
  private

  type, public :: FieldSet_MH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iFieldSet    = 0, &
      nFields      = 0, &
      nVectors     = 0, &
      nStreams     = 0
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    logical ( KDL ) :: &
      DeviceMemory, &
      PinnedMemory, &
      DevicesCommunicate
    character ( LDF ) :: &
      Name = '', &
      Type = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector
    class ( Manifold_H_Form ), pointer :: &
      Manifold => null ( )
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    procedure, public, pass :: &
      Show => Show_FSM
    final :: &
      Finalize
  end type FieldSet_MH_Form

  type, public :: FieldSet_MH_Pointer
    class ( FieldSet_MH_Form ), pointer :: &
      Pointer => null ( )
  end type FieldSet_MH_Pointer


contains


  subroutine Initialize_H &
               ( FSM, M, Name, nFields, FieldOption, VectorOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption )

    class ( FieldSet_MH_Form ), intent ( inout ) :: &
      FSM
    class ( Manifold_H_Form ), intent ( inout ), target :: &
      M
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ) :: &
      nFields
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    integer ( KDI ) :: &
      iF, &  !-- iField
      iV     !-- iVector
    character ( 2 ) :: &
      FieldNumber, &
      VectorNumber

    FSM % IGNORABILITY  =  M % IGNORABILITY

    if ( FSM % Type  ==  '' ) &
      FSM % Type  =  'a FieldSet_M' 
    
    FSM % DeviceMemory  =  .false.
    if ( present ( DeviceMemoryOption ) ) &
      FSM % DeviceMemory  =  DeviceMemoryOption
    
    FSM % PinnedMemory  =  .false.
    if ( present ( PinnedMemoryOption ) ) &
      FSM % PinnedMemory  =  PinnedMemoryOption
    
    FSM % DevicesCommunicate  =  .false.
    if ( present ( DevicesCommunicateOption ) )  &
      FSM % DevicesCommunicate  =  DevicesCommunicateOption  

    FSM % Name  =  Name

    call Show ( 'Initializing ' // trim ( FSM % Type ), FSM % IGNORABILITY )
    call Show ( FSM % Name, 'Name', FSM % IGNORABILITY )
   
      M % nFieldSets  =  M % nFieldSets  +  1
    FSM % iFieldSet   =  M % nFieldSets

    associate ( nF  =>  FSM % nFields )
    nF  =  nFields
    allocate ( FSM % Field ( nF ) )
    allocate ( FSM % Unit ( nF ) )
    if ( present ( FieldOption ) ) then
      FSM % Field  =  FieldOption
    else
      do iF  =  1, nF
        write ( FieldNumber, fmt = '(i2.2)' ) iF
        FSM % Field ( iF )  =  'Field_' // FieldNumber
      end do  !-- iF
    end if  !-- FieldOption
    if ( present ( UnitOption ) ) &
      FSM % Unit  =  UnitOption
    end associate  !-- nF

    if ( present ( VectorIndicesOption ) ) then

      associate ( nV  =>  FSM % nVectors )
      nV  =  size ( VectorIndicesOption )

      allocate ( FSM % VectorIndices ( nV ) )
      do iV  =  1, nV
        call FSM % VectorIndices ( iV ) % Initialize &
               ( VectorIndicesOption ( iV ) )
      end do  !-- iV

      allocate ( FSM % Vector ( nV ) )
      if ( present ( VectorOption ) ) then
        FSM % Vector  =  VectorOption
      else
        do iV  =  1, nV
          write ( VectorNumber, fmt = '(i2.2)' ) iV
          FSM % Vector ( iV )  =  'Vector_' // VectorNumber
        end do  !-- iV
      end if  !-- VectorOption 
      end associate  !-- nV

    end if  !-- VectorIndicesOption

    FSM % Manifold  =>  M

  end subroutine Initialize_H


  subroutine Show_FSM ( FSM )

    class ( FieldSet_MH_Form ), intent ( in ) :: &
      FSM

    integer ( KDI ) :: &
      iF, &  !-- iField
      iV     !-- iVector
    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( FSM % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FSM % IGNORABILITY )

    associate ( M  =>  FSM % Manifold )
    call Show ( FSM % Name,      'Name',      FSM % IGNORABILITY )
    call Show (   M % Name,      'Manifold',  FSM % IGNORABILITY )
    call Show ( FSM % iFieldSet, 'iFieldSet', FSM % IGNORABILITY )
    end associate  !-- M

    call Show ( FSM % nFields, 'nFields', FSM % IGNORABILITY )
    do iF  =  1, FSM % nFields
      call Show ( FSM % Field ( iF ), 'Field',  FSM % IGNORABILITY )
      call Show ( iF,                 'iField', FSM % IGNORABILITY ) 
      call Show ( FSM % Unit ( iF ),  'Unit',   FSM % IGNORABILITY )
    end do !-- iF
    
    call Show ( FSM % nVectors, 'nVectors', FSM % IGNORABILITY )
    do iV  =  1, FSM % nVectors
      call Show ( FSM % Vector ( iV ),                 'Vector', &
                  FSM % IGNORABILITY )
      call Show ( FSM % VectorIndices ( iV ) % Value, 'VectorIndices', &
                  FSM % IGNORABILITY )
    end do  !-- iV

    call Show ( FSM % DeviceMemory,       'DeviceMemory', &
                FSM % IGNORABILITY )
    call Show ( FSM % DeviceMemory,       'DeviceMemory', &
                FSM % IGNORABILITY )
    call Show ( FSM % DevicesCommunicate, 'DevicesCommunicate', &
                FSM % IGNORABILITY )

  end subroutine Show_FSM


  impure elemental subroutine Finalize ( FSM )

    type ( FieldSet_MH_Form ), intent ( inout ) :: &
      FSM

    nullify ( FSM % Manifold )

    if ( allocated ( FSM % Vector ) ) &
      deallocate ( FSM % Vector )
    if ( allocated ( FSM % Field ) ) &
      deallocate ( FSM % Field )
    if ( allocated ( FSM % Unit ) ) &
      deallocate ( FSM % Unit )
    if ( allocated ( FSM % VectorIndices ) ) &
      deallocate ( FSM % VectorIndices )

    call Show ( 'Finalizing ' // trim ( FSM % Type ), FSM % IGNORABILITY )
    call Show ( FSM % Name, 'Name', FSM % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_MH__Form
