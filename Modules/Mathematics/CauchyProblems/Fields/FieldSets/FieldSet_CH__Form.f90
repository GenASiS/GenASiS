module FieldSet_CH__Form

  !-- FieldSet_ChartHeader__Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: FieldSet_CH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFields      = 0, &
      nVectors     = 0
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Type = '', &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector
    class ( Chart_H_Form ), pointer :: &
      Chart => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate
    generic, public :: &
      Initialize => InitializeAllocate
    procedure, public, pass :: &
      Show => Show_FSC
    final :: &
      Finalize
  end type FieldSet_CH_Form

!   ! type, public :: FieldSet_CH_Pointer
!   !   class ( FieldSet_CH_Form ), pointer :: &
!   !     Pointer => null ( )
!   ! end type FieldSet_CH_Pointer


contains


  subroutine InitializeAllocate &
               ( FSC, C, nFields, FieldOption, VectorOption, NameOption, &
                 UnitOption, VectorIndicesOption )

    class ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_H_Form ), intent ( inout ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nFields
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
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

    FSC % IGNORABILITY  =  C % IGNORABILITY

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a FieldSet_C' 
    
    FSC % Name  =  'FieldSet'
    if ( present ( NameOption ) ) &
      FSC % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( FSC % Type ), FSC % IGNORABILITY )
    call Show ( FSC % Name, 'Name', FSC % IGNORABILITY )
   
    FSC % Chart  =>  C

    associate ( nF  =>  FSC % nFields )
    nF  =  nFields
    allocate ( FSC % Field ( nF ) )
    allocate ( FSC % Unit ( nF ) )
    if ( present ( FieldOption ) ) then
      FSC % Field  =  FieldOption
    else
      do iF  =  1, nF
        write ( FieldNumber, fmt = '(i2.2)' ) iF
        FSC % Field ( iF )  =  'Field_' // FieldNumber
      end do  !-- iF
    end if  !-- FieldOption
    if ( present ( UnitOption ) ) &
      FSC % Unit  =  UnitOption
    end associate  !-- nF

    if ( present ( VectorIndicesOption ) ) then

      associate ( nV  =>  FSC % nVectors )
      nV  =  size ( VectorIndicesOption )

      allocate ( FSC % VectorIndices ( nV ) )
      do iV  =  1, nV
        call FSC % VectorIndices ( iV ) % Initialize &
               ( VectorIndicesOption ( iV ) )
      end do  !-- iV

      allocate ( FSC % Vector ( nV ) )
      if ( present ( VectorOption ) ) then
        FSC % Vector  =  VectorOption
      else
        do iV  =  1, nV
          write ( VectorNumber, fmt = '(i2.2)' ) iV
          FSC % Vector ( iV )  =  'Vector_' // VectorNumber
        end do  !-- iV
      end if  !-- VectorOption 
      end associate  !-- nV

    end if  !-- VectorIndicesOption

  end subroutine InitializeAllocate


  subroutine Show_FSC ( FSC )

    class ( FieldSet_CH_Form ), intent ( in ) :: &
      FSC

    integer ( KDI ) :: &
      iF, &  !-- iField
      iV     !-- iVector
    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( FSC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FSC % IGNORABILITY )

    associate ( C  =>  FSC % Chart )
    call Show ( FSC % Name, 'Name',  FSC % IGNORABILITY )
    call Show (   C % Name, 'Chart', FSC % IGNORABILITY )
    end associate  !-- C

    call Show ( FSC % nFields, 'nFields', FSC % IGNORABILITY )
    do iF  =  1, FSC % nFields
      call Show ( iF,                 'iField', FSC % IGNORABILITY ) 
      call Show ( FSC % Field ( iF ), 'Field',  FSC % IGNORABILITY )
      call Show ( FSC % Unit ( iF ),  'Unit',   FSC % IGNORABILITY )
    end do !-- iF
    
    call Show ( FSC % nVectors, 'nVectors', FSC % IGNORABILITY )
    do iV  =  1, FSC % nVectors
      call Show ( FSC % Vector ( iV ),                 'Vector', &
                  FSC % IGNORABILITY )
      call Show ( FSC % VectorIndices ( iV ) % Value, 'VectorIndices', &
                  FSC % IGNORABILITY )
    end do  !-- iV

  end subroutine Show_FSC


  impure elemental subroutine Finalize ( FSC )

    type ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC

    nullify ( FSC % Chart )

    if ( allocated ( FSC % Vector ) ) &
      deallocate ( FSC % Vector )
    if ( allocated ( FSC % Field ) ) &
      deallocate ( FSC % Field )
    if ( allocated ( FSC % Unit ) ) &
      deallocate ( FSC % Unit )
    if ( allocated ( FSC % VectorIndices ) ) &
      deallocate ( FSC % VectorIndices )

    call Show ( 'Finalizing ' // trim ( FSC % Type ), FSC % IGNORABILITY )
    call Show ( FSC % Name, 'Name', FSC % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_CH__Form
