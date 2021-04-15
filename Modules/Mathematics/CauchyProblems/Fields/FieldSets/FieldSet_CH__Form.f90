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
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected
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
    class ( FieldSet_CH_Form ), pointer :: &
      Primary => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_H
    procedure, private, pass :: &
      InitializeClone_H
    generic, public :: &
      Initialize_H => InitializeAllocate_H, InitializeClone_H
    procedure, public, pass :: &
      Show => Show_FSC
    final :: &
      Finalize_FS
  end type FieldSet_CH_Form

  type, public :: FieldSet_C_Element
    !-- FieldSet_Chart_Element
    class ( FieldSet_CH_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type FieldSet_C_Element

  type, public :: FieldSet_C_Pointer
    !-- FieldSet_Chart_Pointer
    class ( FieldSet_CH_Form ), pointer :: &
      Pointer => null ( )
  end type FieldSet_C_Pointer


contains


  subroutine InitializeAllocate_H &
               ( FSC, C, FieldOption, VectorOption, NameOption, UnitOption, &
                 VectorIndicesOption, nFieldsOption )

    class ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      iF, &  !-- iField
      iV     !-- iVector
    character ( 2 ) :: &
      FieldNumber, &
      VectorNumber

    FSC % IGNORABILITY  =  C % IGNORABILITY

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a FieldSet_C' 
    
    FSC % Name  =  'Fields'
    if ( present ( NameOption ) ) &
      FSC % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( FSC % Type ), FSC % IGNORABILITY )
    call Show ( FSC % Name, 'Name', FSC % IGNORABILITY )
   
    FSC % Chart  =>  C

    associate ( nF  =>  FSC % nFields )
    nF  =  1
    if ( present ( nFieldsOption ) ) &
      nF  =  nFieldsOption
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
    if ( .not. allocated ( FSC % iaSelected ) ) then
      allocate ( FSC % iaSelected ( nF ) )
      FSC % iaSelected  =  [ ( iF, iF = 1, nF ) ]       
    end if
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

  end subroutine InitializeAllocate_H


  subroutine InitializeClone_H &
               ( FSC_T, FSC_S, NameOption, iaSelectedOption )

    class ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC_T  !-- FSC_Target
    class ( FieldSet_CH_Form ), intent ( in ), target :: &
      FSC_S  !-- FSC_Source
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iV_S, &  !-- iVector
      nV_T   !-- nVectors_T
    integer ( KDI ), dimension ( 3 ) :: &
      iaV_T
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices_T
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Vector_T

    FSC_T % Primary  =>  FSC_S

    associate ( nF_S  =>  FSC_S % nFields )

    if ( present ( iaSelectedOption ) ) then
      allocate ( FSC_T % iaSelected, source = iaSelectedOption )
    else
      allocate ( FSC_T % iaSelected ( nF_S ) )
      FSC_T % iaSelected  =  [ ( iS, iS = 1, nF_S ) ]
    end if

    Name  =  FSC_S % Name
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Count vectors are among the selected
    nV_T  =  0
    do iV_S  =  1, FSC_S % nVectors
      associate ( iaV_S  =>  FSC_S % VectorIndices ( iV_S ) % Value )
      iaV_T ( 1 )  =  findloc ( FSC_T % iaSelected, iaV_S ( 1 ), dim = 1 )
      iaV_T ( 2 )  =  findloc ( FSC_T % iaSelected, iaV_S ( 2 ), dim = 1 )
      iaV_T ( 3 )  =  findloc ( FSC_T % iaSelected, iaV_S ( 3 ), dim = 1 )
      if ( all ( iaV_T > 0 ) ) &
        nV_T  =  nV_T + 1
      end associate !-- iaV
    end do !-- iV

    allocate ( Vector_T ( nV_T ) )
    allocate ( VectorIndices_T ( nV_T ) )

    !-- Populate vector names and indices
    nV_T  =  0
    do iV_S  =  1, FSC_S % nVectors
      associate ( iaV_S  =>  FSC_S % VectorIndices ( iV_S ) % Value )
      iaV_T ( 1 )  =  findloc ( FSC_T % iaSelected, iaV_S ( 1 ), dim = 1 )
      iaV_T ( 2 )  =  findloc ( FSC_T % iaSelected, iaV_S ( 2 ), dim = 1 )
      iaV_T ( 3 )  =  findloc ( FSC_T % iaSelected, iaV_S ( 3 ), dim = 1 )
      if ( all ( iaV_T > 0 ) ) then
        nV_T  =  nV_T + 1
        Vector_T ( nV_T )  =  FSC_S % Vector ( iV_S )
        call VectorIndices_T ( nV_T ) % Initialize ( iaV_T )
      end if
      end associate !-- iaV
    end do !-- iV

    call FSC_T % Initialize_H &
           ( C = FSC_S % Chart, &
             FieldOption = FSC_S % Field, &
             VectorOption = Vector_T, &
             NameOption = Name, &
             UnitOption = FSC_S % Unit, &
             VectorIndicesOption = VectorIndices_T, &
             nFieldsOption = size ( FSC_T % iaSelected ) )

    end associate !-- nF_S

  end subroutine InitializeClone_H


  subroutine Show_FSC ( FSC )

    class ( FieldSet_CH_Form ), intent ( in ) :: &
      FSC

    integer ( KDI ) :: &
      iF, &  !-- iField
      iS, &  !-- iSelected
      iV     !-- iVector
    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( FSC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FSC % IGNORABILITY )

    call Show ( FSC % Name, 'Name',  FSC % IGNORABILITY )
    if ( associated ( FSC % Primary ) ) &
      call Show ( FSC % Primary % Name, 'Primary', FSC % IGNORABILITY )

    call Show ( FSC % Chart % Name, 'Chart', FSC % IGNORABILITY )

    call Show ( FSC % nFields, 'nFields', FSC % IGNORABILITY )
    do iS  =  1, FSC % nFields
      iF  =  FSC % iaSelected ( iS )
      call Show ( iS,                 'iField', FSC % IGNORABILITY ) 
      call Show ( FSC % Field ( iF ), 'Field',     FSC % IGNORABILITY )
      call Show ( FSC % Unit ( iF ),  'Unit',      FSC % IGNORABILITY )
    end do !-- iF
    
    call Show ( FSC % nVectors, 'nVectors', FSC % IGNORABILITY )
    do iV  =  1, FSC % nVectors
      call Show ( FSC % Vector ( iV ),                 'Vector', &
                  FSC % IGNORABILITY )
      call Show ( FSC % VectorIndices ( iV ) % Value, 'VectorIndices', &
                  FSC % IGNORABILITY )
    end do  !-- iV

  end subroutine Show_FSC


  impure elemental subroutine Finalize_FS ( FSC )

    type ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC

    nullify ( FSC % Primary )
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
   
  end subroutine Finalize_FS


  impure elemental subroutine Finalize_E ( FSE )
    
    type ( FieldSet_C_Element ), intent ( inout ) :: &
      FSE

    if ( allocated ( FSE % Element ) ) &
      deallocate ( FSE % Element )

  end subroutine Finalize_E


end module FieldSet_CH__Form
