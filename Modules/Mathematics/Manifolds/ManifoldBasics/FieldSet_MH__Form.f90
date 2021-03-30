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
                 DevicesCommunicateOption, nVectorsOption )

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
    integer ( KDI ), intent ( in ), optional :: &
      nVectorsOption

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
    if ( present ( FieldOption ) ) then
      FSM % Field  =  FieldOption
    else
      do iF  =  1, nF
        write ( FieldNumber, fmt = '(i2.2)' ) iF
        FSM % Field ( iF )  =  'Field_' // FieldNumber
      end do  !-- iF
    end if  !-- FieldOption
    end associate  !-- nF

    associate ( nV  =>  FSM % nVectors )
    if ( present ( nVectorsOption ) ) then
      nV  =  nVectorsOption
    else
      nV  =  0
    end if !-- nVectorsOption
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

    FSM % Manifold  =>  M

  end subroutine Initialize_H


  subroutine Show_FSM ( FSM )

    class ( FieldSet_MH_Form ), intent ( in ) :: &
      FSM

    integer ( KDI ) :: &
      iF  !-- iField
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
    call Show ( FSM % Field, 'Fields', FSM % IGNORABILITY )
    
    call Show ( FSM % nVectors, 'nVectors', FSM % IGNORABILITY )
!    do iV  =  1, FSM % nVectors
!      call Show ( FSM % Vector ( iF ), 
!    end do !-- iV

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

    call Show ( 'Finalizing ' // trim ( FSM % Type ), FSM % IGNORABILITY )
    call Show ( FSM % Name, 'Name', FSM % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_MH__Form
