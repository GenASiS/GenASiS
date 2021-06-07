module FieldSet_C__Form

  !-- FieldSet_Chart__Form

  use Basics
  use Manifolds
  use Storage_FSC__Form
  use GhostExchange_FSC__Form

  implicit none
  private

  type, public :: FieldSet_C_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFields      = 0, &
      nVectors     = 0, &
      nBoundaries  = 0
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
      Vector, &
      Boundary
    character ( LDL ), dimension ( :, : ), allocatable :: &
      BoundaryCondition
    class ( Storage_FSC_Form ), allocatable :: &
      Storage_FSC
    class ( GhostExchange_FSC_Form ), allocatable :: &
      GhostExchange_FSC
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

    private :: &
      SetDefaultBoundaries, &
      ShowBoundaryConditions

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

    if ( associated ( FSC % Primary ) ) then
      allocate &
        ( FSC % Boundary, source = FSC % Primary % Boundary )
      allocate &
        ( FSC % BoundaryCondition, source = FSC % Primary % BoundaryCondition )
    else
      call SetDefaultBoundaries ( FSC )
    end if

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

    integer ( KDI ) :: &
      iB  !-- iBoundary

    associate ( C  =>  FSC % Chart )

    if ( C % Connectivity % nFaces == 0 ) then
      call Show ( 'Faces not included in Connectivity', CONSOLE % ERROR )
      call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsFace', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iDimension > C % nDimensions ) then
      call Show ( 'Selected iDimension > nDimensions', CONSOLE % ERROR )
      call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsFace', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    iB = 1
    if ( present ( iBoundaryOption ) ) then
      if ( iBoundaryOption > FSC % nBoundaries ) then
        call Show ( 'Selected iBoundary > nBoundaries', CONSOLE % ERROR )
        call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( iBoundaryOption == 1 ) then
        if ( present ( BoundaryOption ) ) then
          call Show ( 'Boundary name not allowed for iBoundary == 1', &
                      CONSOLE % ERROR )
          call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      else
        if ( .not. present ( BoundaryOption ) ) then
          call Show ( 'Boundary name required for iBoundary > 1', &
                      CONSOLE % ERROR )
          call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      end if
      iB = iBoundaryOption
    end if

    if ( present ( BoundaryOption ) ) then
      if ( .not.present ( iBoundaryOption ) ) then
        call Show ( 'Argument iBoundary required when Boundary name present', &
                    CONSOLE % ERROR )
        call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )         
      else
        FSC % Boundary ( iBoundaryOption ) = BoundaryOption
      end if
    end if

    associate &
      ( Cy  => C % Connectivity, &
        iD => iDimension )
    FSC % BoundaryCondition ( Cy % iaInner ( iD ), iB ) &
      = BoundaryCondition ( 1 )
    FSC % BoundaryCondition ( Cy % iaOuter ( iD ), iB ) &
      = BoundaryCondition ( 2 )
    end associate !-- Cy, etc.

    end associate !-- C

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

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- jDimension, etc.
      iB  !-- iBoundary

    associate ( C  =>  FSC % Chart )

    iD = iDimension
    jD = mod ( iD, 3 ) + 1
    kD = mod ( jD, 3 ) + 1
 
    if ( C % Connectivity % nEdges == 0 ) then
      call Show ( 'Edges not included in Connectivity', CONSOLE % ERROR )
      call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsEdge', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( jD > C % nDimensions .or. kD > C % nDimensions ) then
      call Show ( 'Selected jDimension or kDimension > nDimensions', &
                  CONSOLE % ERROR )
      call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsEdge', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    iB = 1
    if ( present ( iBoundaryOption ) ) then
      if ( iBoundaryOption > FSC % nBoundaries ) then
        call Show ( 'Selected iBoundary > nBoundaries', CONSOLE % ERROR )
        call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( iBoundaryOption == 1 ) then
        if ( present ( BoundaryOption ) ) then
          call Show ( 'Boundary name not allowed for iBoundary == 1', &
                      CONSOLE % ERROR )
          call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      else
        if ( .not. present ( BoundaryOption ) ) then
          call Show ( 'Boundary name required for iBoundary > 1', &
                      CONSOLE % ERROR )
          call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      end if
      iB = iBoundaryOption
    end if

    if ( present ( BoundaryOption ) ) then
      if ( .not.present ( iBoundaryOption ) ) then
        call Show ( 'Argument iBoundary required when BoundaryName present', &
                    CONSOLE % ERROR )
        call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )         
      else
        FSC % Boundary ( iBoundaryOption ) = BoundaryOption
      end if
    end if

    associate ( Cy => C % Connectivity )
    FSC % BoundaryCondition ( Cy % iaInnerInner ( iD ), iB ) &
      = BoundaryCondition ( 1 )
    FSC % BoundaryCondition ( Cy % iaOuterInner ( iD ), iB ) &
      = BoundaryCondition ( 2 )
    FSC % BoundaryCondition ( Cy % iaInnerOuter ( iD ), iB ) &
      = BoundaryCondition ( 3 )
    FSC % BoundaryCondition ( Cy % iaOuterOuter ( iD ), iB ) &
      = BoundaryCondition ( 4 )
    end associate !-- Cy

    end associate !-- C

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

    call ShowBoundaryConditions ( FSC )

  end subroutine Show_FS


  subroutine Clear_FS ( FSC )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC

    associate ( FSV  =>  FSC % Storage_FSC % Storage % Value )
    call Clear ( FSV, UseDeviceOption = FSC % Storage_FSC % DeviceMemory )
    end associate !-- FSV

  end subroutine Clear_FS


  subroutine Copy_FS ( FSC_S, FSC_T )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC_S, &
      FSC_T

    associate &
      ( FSV_S  =>  FSC_S % Storage_FSC % Storage % Value, &
        FSV_T  =>  FSC_T % Storage_FSC % Storage % Value )
    call Copy ( FSV_S, FSV_T, &
                UseDeviceOption = FSC_S % Storage_FSC % DeviceMemory )
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

    if ( allocated ( FSC % GhostExchange_FSC ) ) &
      deallocate ( FSC % GhostExchange_FSC )
    if ( allocated ( FSC % Storage_FSC ) ) &
      deallocate ( FSC % Storage_FSC )
    if ( allocated ( FSC % BoundaryCondition ) ) &
      deallocate ( FSC % BoundaryCondition )
    if ( allocated ( FSC % Boundary ) ) &
      deallocate ( FSC % Boundary )
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


  subroutine SetDefaultBoundaries ( FSC, nExcisionsOption )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      nExcisionsOption

    FSC % nBoundaries = 1
    if ( present ( nExcisionsOption ) ) &
      FSC % nBoundaries = 1 + nExcisionsOption 

    associate &
      ( Cy => FSC % Chart % Connectivity )
    allocate &
      ( FSC % BoundaryCondition ( Cy % nConnections, FSC % nBoundaries ) )
    allocate &
      ( FSC % Boundary ( FSC % nBoundaries ) )

    FSC % Boundary = ''
    FSC % Boundary ( 1 ) = 'Extent' 

    FSC % BoundaryCondition = ''
    FSC % BoundaryCondition ( :, 1 ) = 'PERIODIC'

    end associate !-- C

  end subroutine SetDefaultBoundaries


  subroutine ShowBoundaryConditions ( FSC )

    class ( FieldSet_C_Form ), intent ( in ) :: &
      FSC

    integer ( KDI ) :: &
      iB, &  !-- iBoundary
      iD, jD, kD  !-- iDimension, etc.

    associate &
      ( Cy  => FSC % Chart % Connectivity, &
        BC => FSC % BoundaryCondition ( :, : ), &
        BN => FSC % Boundary ( : ), &
        nD => FSC % Chart % nDimensions )

    call Show ( 'Boundary conditions', FSC % IGNORABILITY )
    call Show ( FSC % nBoundaries, 'nBoundaries', FSC % IGNORABILITY )

    do iB = 1, FSC % nBoundaries
      call Show ( BN ( iB ), 'Boundary', FSC % IGNORABILITY )
      call Show ( iB, 'iBoundary', FSC % IGNORABILITY )
  
      if ( Cy % nFaces > 0 ) then
          do iD = 1, nD
            call Show ( iD, 'Faces, iDimension', &
                        FSC % IGNORABILITY )
            associate &
              ( iaI => Cy % iaInner ( iD ), &
                iaO => Cy % iaOuter ( iD ) )
            call Show ( [ BC ( iaI, iB ), BC ( iaO, iB ) ], &
                        '[ Inner, Outer ]', FSC % IGNORABILITY )
            end associate !-- iaI, etc.
          end do !-- iD
      end if

      if ( Cy % nEdges > 0 ) then
          do iD = 1, nD
            jD = mod ( iD, 3 ) + 1
            kD = mod ( jD, 3 ) + 1
            if ( jD > nD .or. kD > nD ) &
              cycle
            call Show ( iD, 'Edges parallel to iDimension', &
                        FSC % IGNORABILITY )
            associate &
              ( iaII => Cy % iaInnerInner ( iD ), &
                iaOI => Cy % iaOuterInner ( iD ), &
                iaIO => Cy % iaInnerOuter ( iD ), &
                iaOO => Cy % iaOuterOuter ( iD ) )
            call Show ( [ BC ( iaII, iB ), BC ( iaOI, iB ), &
                          BC ( iaIO, iB ), BC ( iaOO, iB ) ], &
                        '[ InnerInner, OuterInner, InnerOuter, OuterOuter ]', &
                        FSC % IGNORABILITY )
            end associate !-- iaII, etc.
          end do !-- iD
      end if

    end do !-- iB

    end associate !-- Cy, etc.

  end subroutine ShowBoundaryConditions


end module FieldSet_C__Form
