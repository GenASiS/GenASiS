module FieldSet_C__Form

  !-- FieldSet_Chart__Form

  use Basics
  use Manifolds
  use Storage_FSC__Form
  use GhostExchange_FSC__Form
  use Boundaries_FSC__Form

  implicit none
  private

  type, public :: FieldSet_C_Form
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
    class ( Storage_FSC_Form ), allocatable :: &
      Storage_FSC
    class ( GhostExchange_FSC_Form ), allocatable :: &
      GhostExchange_FSC
    class ( Boundaries_FSC_Form ), allocatable :: &
      Boundaries_FSC
    class ( Chart_H_Form ), pointer :: &
      Chart => null ( )
    class ( FieldSet_C_Form ), pointer :: &
      Primary => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, private, pass :: &
      InitializeClone
    generic, public :: &
      Initialize => InitializeAllocate_FS, InitializeClone
    procedure, public, pass :: &
      SetBoundaryConditionsFace
    procedure, public, pass :: &
      SetBoundaryConditionsEdge
    procedure, private, pass :: &
      Show_FS
    generic, public :: &
      Show => Show_FS
    procedure, public, pass :: &
      Clear => Clear_FS
    procedure, public, pass :: &
      Copy => Copy_FS
    procedure, public, pass :: &
      ExchangeGhostData
    procedure, public, pass :: &
      StartGhostExchange
    procedure, public, pass :: &
      FinishGhostExchange
    procedure, public, pass :: &
      ApplyBoundaryConditions
    procedure, public, pass :: &
      UpdateDevice => UpdateDevice_FS
    procedure, public, pass :: &
      UpdateHost => UpdateHost_FS
    final :: &
      Finalize_FS
  end type FieldSet_C_Form

  type, public :: FieldSet_C_Element
    !-- FieldSet_Chart_Element
    class ( FieldSet_C_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type FieldSet_C_Element


contains


  subroutine InitializeAllocate_FS &
               ( FSC, C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
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
      iF, &  !-- iField
      iV     !-- iVector
    character ( 2 ) :: &
      FieldNumber, &
      VectorNumber

    FSC % IGNORABILITY  =  C % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FSC % IGNORABILITY  =  IgnorabilityOption

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

    allocate ( FSC % Storage_FSC )
    associate ( SFSC => FSC % Storage_FSC )    
    if ( associated ( FSC % Primary ) ) then
      associate ( SFSC_S  =>  FSC % Primary % Storage_FSC ) 
      call SFSC % Initialize &
            ( SFSC_S, FSC % Vector, FSC % Name, FSC % VectorIndices, &
              FSC % iaSelected )
      end associate
    else
      call SFSC % Initialize &
             ( C, FSC % Field, FSC % Vector, FSC % Name, FSC % Unit, &
               FSC % VectorIndices, FSC % nFields, DeviceMemoryOption, &
               PinnedMemoryOption )
    end if
    end associate !-- SFSC

    allocate ( FSC % GhostExchange_FSC )
    associate ( GE => FSC % GhostExchange_FSC )    
    call GE % Initialize &
           ( DevicesCommunicateOption )
    end associate !-- GE

    allocate ( FSC % Boundaries_FSC )
    associate ( B => FSC % Boundaries_FSC )
    if ( associated ( FSC % Primary ) ) then
      call B % Initialize ( FSC % Primary % Boundaries_FSC )
    else
      call B % Initialize ( C )
    end if
    end associate !-- B

  end subroutine InitializeAllocate_FS


  subroutine InitializeClone &
               ( FSC_T, FSC_S, iaSelected, NameOption, IgnorabilityOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC_T  !-- FSC_Target
    class ( FieldSet_C_Form ), intent ( in ), target :: &
      FSC_S  !-- FSC_Source
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaSelected
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iV_S, &  !-- iVector
      nV_T   !-- nVectors_T
    integer ( KDI ), dimension ( 3 ) :: &
      iaV_T
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices_T
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Vector_T

    if ( associated ( FSC_S % Primary ) ) then
      FSC_T % Primary  =>  FSC_S % Primary
    else
      FSC_T % Primary  =>  FSC_S
    end if

    associate ( nF_S  =>  FSC_S % nFields )

    allocate ( FSC_T % iaSelected, source = iaSelected )

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
        call VectorIndices_T ( nV_T ) % Initialize ( iaV_S )
      end if
      end associate !-- iaV
    end do !-- iV

    call FSC_T % Initialize &
           ( C = FSC_S % Chart, &
             FieldOption = FSC_S % Field, &
             VectorOption = Vector_T, &
             NameOption = Name, &
             UnitOption = FSC_S % Unit, &
             VectorIndicesOption = VectorIndices_T, &
             nFieldsOption = size ( FSC_T % iaSelected ), &
             IgnorabilityOption = IgnorabilityOption )

    end associate !-- nF_S

  end subroutine InitializeClone


  subroutine SetBoundaryConditionsFace &
               ( FSC, BoundaryCondition, iDimension, BoundaryOption, &
                 iBoundaryOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    character ( * ), dimension ( 2 ), intent ( in ) :: &
      BoundaryCondition  !-- [ Inner, Outer ]
    integer ( KDI ), intent ( in ) :: &
      iDimension
    character ( * ), intent ( in ), optional :: &
      BoundaryOption
    integer ( KDI ), intent ( in ), optional :: &
      iBoundaryOption

    if ( iDimension  <=  FSC % Chart % nDimensions ) &
      call FSC % Boundaries_FSC % SetFace &
             ( FSC % Chart, BoundaryCondition, iDimension, BoundaryOption, &
               iBoundaryOption )

  end subroutine SetBoundaryConditionsFace


  subroutine SetBoundaryConditionsEdge &
               ( FSC, BoundaryCondition, iDimension, BoundaryOption, &
                 iBoundaryOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    character ( * ), dimension ( 4 ), intent ( in ) :: &
      BoundaryCondition  !-- [ InnerInner, OuterInner, InnerOuter, OuterOuter ]
    integer ( KDI ), intent ( in ) :: &
      iDimension
    character ( * ), intent ( in ), optional :: &
      BoundaryOption
    integer ( KDI ), intent ( in ), optional :: &
      iBoundaryOption

    call FSC % Boundaries_FSC % SetEdge &
           ( FSC % Chart, BoundaryCondition, iDimension, BoundaryOption, &
             iBoundaryOption )

  end subroutine SetBoundaryConditionsEdge


  subroutine Show_FS ( FSC )

    class ( FieldSet_C_Form ), intent ( in ) :: &
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
      call Show ( iS,                 'iSelected', FSC % IGNORABILITY ) 
      call Show ( iF,                 'iField',    FSC % IGNORABILITY ) 
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

    call Show ( FSC % Storage_FSC % DeviceMemory, &
                'DeviceMemory', FSC % IGNORABILITY ) 
    call Show ( FSC % Storage_FSC % PinnedMemory, &
                'PinnedMemory', FSC % IGNORABILITY ) 
    call Show ( FSC % GhostExchange_FSC % DevicesCommunicate, &
                'DevicesCommunicate', FSC % IGNORABILITY ) 

    call FSC % Boundaries_FSC % Show ( FSC % Chart, FSC % IGNORABILITY )

  end subroutine Show_FS


  subroutine Clear_FS ( FSC )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField

    associate ( FSV  =>  FSC % Storage_FSC % Storage % Value )
    do iS  =  1,  FSC % nFields
      iF  =  FSC % iaSelected ( iS )
      call Clear ( FSV ( :, iF ), &
                   UseDeviceOption = FSC % Storage_FSC % DeviceMemory )
    end do !-- iS
    end associate !-- FSV

  end subroutine Clear_FS


  subroutine Copy_FS ( FSC_S, FSC_T )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC_S, &
      FSC_T

    integer ( KDI ) :: &
      iS, &    !-- iSelected
      iF_S, &  !-- iField
      iF_T

    associate &
      ( FSV_S  =>  FSC_S % Storage_FSC % Storage % Value, &
        FSV_T  =>  FSC_T % Storage_FSC % Storage % Value )
    do iS  =  1,  FSC_S % nFields
      iF_S  =  FSC_S % iaSelected ( iS )
      iF_T  =  FSC_T % iaSelected ( iS )
      call Copy ( FSV_S ( :, iF_S ), FSV_T ( :, iF_T ), &
                  UseDeviceOption = FSC_S % Storage_FSC % DeviceMemory )
    end do
    end associate !-- FSV_S, etc.

  end subroutine Copy_FS


  subroutine ExchangeGhostData ( FSC, TimerLevelOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    associate &
      (   GE  =>  FSC % GhostExchange_FSC, &
        SFSC  =>  FSC % Storage_FSC, &
           C  =>  FSC % Chart )
    call GE % Exchange ( SFSC, C, TimerLevelOption )
    end associate !-- GE, etc.

  end subroutine ExchangeGhostData


  subroutine StartGhostExchange ( FSC, TimerLevelOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    associate &
      (   GE  =>  FSC % GhostExchange_FSC, &
        SFSC  =>  FSC % Storage_FSC, &
           C  =>  FSC % Chart )
    call GE % StartExchange ( C, SFSC, TimerLevelOption )
    end associate !-- GE, etc.

  end subroutine StartGhostExchange


  subroutine FinishGhostExchange ( FSC )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC

    associate &
      (   GE  =>  FSC % GhostExchange_FSC, &
        SFSC  =>  FSC % Storage_FSC, &
           C  =>  FSC % Chart )
    call GE % FinishExchange ( SFSC, C )
    end associate !-- GE, etc.

  end subroutine FinishGhostExchange


  subroutine ApplyBoundaryConditions ( FSC )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC

    call FSC % Boundaries_FSC % Apply ( FSC % Storage_FSC, FSC % Chart )

  end subroutine ApplyBoundaryConditions


  subroutine UpdateDevice_FS ( FSC, TimerLevelOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    call FSC % Storage_FSC % UpdateDevice ( TimerLevelOption )

  end subroutine UpdateDevice_FS


  subroutine UpdateHost_FS ( FSC, TimerLevelOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    call FSC % Storage_FSC % UpdateHost ( TimerLevelOption )

  end subroutine UpdateHost_FS


  impure elemental subroutine Finalize_FS ( FSC )

    type ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC

    nullify ( FSC % Primary )
    nullify ( FSC % Chart )

    if ( allocated ( FSC % Boundaries_FSC ) ) &
      deallocate ( FSC % Boundaries_FSC )
    if ( allocated ( FSC % GhostExchange_FSC ) ) &
      deallocate ( FSC % GhostExchange_FSC )
    if ( allocated ( FSC % Storage_FSC ) ) &
      deallocate ( FSC % Storage_FSC )
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


end module FieldSet_C__Form
