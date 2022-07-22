module FieldSet_Form

  use Basics
  use Algebra
  use Manifolds
  use GhostExchange_Form
  use Boundaries_Form

  implicit none
  private

  type, public :: FieldSetForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFields      = 0, &
      nVectors     = 0, &
      nBoundaries  = 0
    integer ( KDI ) :: &
      iTimerGhost    = 0, &
      iTimerGhost_UH = 0, &  !-- UpdateHost
      iTimerGhost_EG = 0, &  !-- Exchange
      iTimerGhost_UD = 0     !-- UpdateDevice
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected
    logical ( KDL ) :: &
      DeviceMemory, &
      PinnedMemory, &
      DevicesCommunicate
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Type = '', &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector
    type ( StorageForm ), pointer :: &
      Storage_GS => null ( )
    type ( StorageForm ), dimension ( : ), allocatable :: &
      Storage
    class ( GhostExchangeForm ), dimension ( : ), allocatable :: &
      GhostExchange
    class ( BoundariesForm ), dimension ( : ), allocatable :: &
      Boundaries
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    class ( FieldSetForm ), pointer :: &
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
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      TimerGhost
    procedure, public, pass :: &
      Clear => Clear_FS
    procedure, public, pass ( FS_S ) :: &
      Copy => Copy_FS
    procedure, private, pass :: &
      MultiplyAdd_FS
    procedure, private, pass :: &
      MultiplyAddInPlace_FS
    generic, public :: &
      MultiplyAdd => MultiplyAdd_FS, MultiplyAddInPlace_FS
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
      Finalize
  end type FieldSetForm

  type, public :: FieldSetElement
    class ( FieldSetForm ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type FieldSetElement


contains


  subroutine InitializeAllocate_FS &
               ( FS, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, AssociateFieldsOption, &
                 UnitOption, VectorIndicesOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( FieldSetForm ), intent ( inout ), target :: &
      FS
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption, &
      AssociateFieldsOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iF, &  !-- iField
      iV     !-- iVector
    character ( 2 ) :: &
      FieldNumber, &
      VectorNumber

    FS % IGNORABILITY  =  A % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FS % IGNORABILITY  =  IgnorabilityOption

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a FieldSet' 
    
    FS % Name  =  'Fields'
    if ( present ( NameOption ) ) &
      FS % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( FS % Type ), FS % IGNORABILITY )
    call Show ( FS % Name, 'Name', FS % IGNORABILITY )

    !-- Metadata

    FS % DeviceMemory  =  .false.
    if ( present ( DeviceMemoryOption ) ) &
      FS % DeviceMemory  =  DeviceMemoryOption
    
    FS % PinnedMemory  =  .false.
    if ( present ( PinnedMemoryOption ) ) &
      FS % PinnedMemory  =  PinnedMemoryOption

    FS % DevicesCommunicate  =  .false.
    if ( present ( DevicesCommunicateOption ) )  &
      FS % DevicesCommunicate  =  DevicesCommunicateOption

    FS % Atlas  =>  A
    associate ( nC  =>  FS % Atlas % nCharts )

    associate ( nF  =>  FS % nFields )
    nF  =  1
    if ( present ( nFieldsOption ) ) &
      nF  =  nFieldsOption
    allocate ( FS % Field ( nF ) )
    allocate ( FS % Unit ( nF, nC ) )
    if ( present ( FieldOption ) ) then
      FS % Field  =  FieldOption
    else
      do iF  =  1, nF
        write ( FieldNumber, fmt = '(i2.2)' ) iF
        FS % Field ( iF )  =  'Field_' // FieldNumber
      end do  !-- iF
    end if  !-- FieldOption
    if ( present ( UnitOption ) ) &
      FS % Unit  =  UnitOption
    if ( .not. allocated ( FS % iaSelected ) ) then
      allocate ( FS % iaSelected ( nF ) )
      FS % iaSelected  =  [ ( iF, iF = 1, nF ) ]       
    end if
    end associate  !-- nF

    if ( present ( VectorIndicesOption ) ) then

      associate ( nV  =>  FS % nVectors )
      nV  =  size ( VectorIndicesOption )

      allocate ( FS % VectorIndices ( nV ) )
      do iV  =  1, nV
        call FS % VectorIndices ( iV ) % Initialize &
               ( VectorIndicesOption ( iV ) )
      end do  !-- iV

      allocate ( FS % Vector ( nV ) )
      if ( present ( VectorOption ) ) then
        FS % Vector  =  VectorOption
      else
        do iV  =  1, nV
          write ( VectorNumber, fmt = '(i2.2)' ) iV
          FS % Vector ( iV )  =  'Vector_' // VectorNumber
        end do  !-- iV
      end if  !-- VectorOption 
      end associate  !-- nV

    end if  !-- VectorIndicesOption

    !-- Storage, GhostExchange, Boundaries

    allocate ( FS % Storage ( nC ) )
    allocate ( FS % GhostExchange ( nC ) )
    allocate ( FS % Boundaries ( nC ) )

    do iC  =  1,  nC
      select type ( C  =>  FS % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      associate ( S  =>  FS % Storage ( iC ) )    
      if ( associated ( FS % Primary ) ) then
        associate ( S_S  =>  FS % Primary % Storage ( iC ) ) 
        call S % Initialize &
               ( S_S, &
                VectorOption = FS % Vector, &
                NameOption = FS % Name, &
                VectorIndicesOption = FS % VectorIndices, &
                iaSelectedOption = FS % iaSelected )
        end associate !-- S_S
      else
        call S % Initialize &
               ( [ C % nCellsLocal, FS % nFields ], &
                 VariableOption = FS % Field, &
                 VectorOption = FS % Vector, &
                 NameOption = FS % Name, &
                 ClearOption = .true., &
                 PinnedOption = FS % PinnedMemory, &
                 UnitOption = FS % Unit ( :, iC ), &
                 VectorIndicesOption = FS % VectorIndices )
        if ( FS % DeviceMemory ) &
          call S % AllocateDevice ( AssociateFieldsOption )
      end if  !-- associated Primary 
      end associate !-- S

      associate ( GE  =>  FS % GhostExchange ( iC ) )    
      call GE % Initialize ( )
      end associate !-- GE

      associate ( B => FS % Boundaries ( iC ) )
      if ( associated ( FS % Primary ) ) then
        call B % Initialize ( FS % Primary % Boundaries ( iC ) )
      else
        call B % Initialize ( C )
      end if
      end associate !-- B

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'FieldSet_Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeAllocate_FS', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- nC

    select type ( A  =>  FS % Atlas )
    class is ( Atlas_SCG_Form )
      FS % nBoundaries  =  FS % Boundaries ( 1 ) % nBoundaries
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'FieldSet_Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeAllocate_FS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    !-- For convenience with a single Chart_GS
    select type ( A  =>  FS % Atlas )
    class is ( Atlas_SCG_Form )
      FS % Storage_GS  =>  FS % Storage ( 1 )
    end select !-- A

  end subroutine InitializeAllocate_FS


  subroutine InitializeClone &
               ( FS_T, FS_S, iaSelected, NameOption, IgnorabilityOption )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_T  !-- FS_Target
    class ( FieldSetForm ), intent ( in ), target :: &
      FS_S  !-- FS_Source
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaSelected
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iS, &    !-- iSelected
      iV_S, &  !-- iVector
      nV_T     !-- nVectors_T
    integer ( KDI ), dimension ( 3 ) :: &
      iaV_T
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices_T
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Vector_T

    if ( associated ( FS_S % Primary ) ) then
      FS_T % Primary  =>  FS_S % Primary
    else
      FS_T % Primary  =>  FS_S
    end if

    associate ( nF_S  =>  FS_S % nFields )

    allocate ( FS_T % iaSelected, source = iaSelected )

    Name  =  FS_S % Name
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Count vectors among the selected
    nV_T  =  0
    do iV_S  =  1, FS_S % nVectors
      associate ( iaV_S  =>  FS_S % VectorIndices ( iV_S ) % Value )
      iaV_T ( 1 )  =  findloc ( FS_T % iaSelected, iaV_S ( 1 ), dim = 1 )
      iaV_T ( 2 )  =  findloc ( FS_T % iaSelected, iaV_S ( 2 ), dim = 1 )
      iaV_T ( 3 )  =  findloc ( FS_T % iaSelected, iaV_S ( 3 ), dim = 1 )
      if ( all ( iaV_T > 0 ) ) &
        nV_T  =  nV_T + 1
      end associate !-- iaV
    end do !-- iV

    allocate ( Vector_T ( nV_T ) )
    allocate ( VectorIndices_T ( nV_T ) )

    !-- Populate vector names and indices
    nV_T  =  0
    do iV_S  =  1, FS_S % nVectors
      associate ( iaV_S  =>  FS_S % VectorIndices ( iV_S ) % Value )
      iaV_T ( 1 )  =  findloc ( FS_T % iaSelected, iaV_S ( 1 ), dim = 1 )
      iaV_T ( 2 )  =  findloc ( FS_T % iaSelected, iaV_S ( 2 ), dim = 1 )
      iaV_T ( 3 )  =  findloc ( FS_T % iaSelected, iaV_S ( 3 ), dim = 1 )
      if ( all ( iaV_T > 0 ) ) then
        nV_T  =  nV_T + 1
        Vector_T ( nV_T )  =  FS_S % Vector ( iV_S )
        call VectorIndices_T ( nV_T ) % Initialize ( iaV_S )
      end if
      end associate !-- iaV
    end do !-- iV

    call FS_T % Initialize &
           ( A = FS_S % Atlas, &
             FieldOption = FS_S % Field, &
             VectorOption = Vector_T, &
             NameOption = Name, &
             DeviceMemoryOption = FS_S % DeviceMemory, &
             PinnedMemoryOption = FS_S % PinnedMemory, &
             DevicesCommunicateOption = FS_S % DevicesCommunicate, &
             UnitOption = FS_S % Unit, &
             VectorIndicesOption = VectorIndices_T, &
             nFieldsOption = size ( FS_T % iaSelected ), &
             IgnorabilityOption = IgnorabilityOption )

    end associate !-- nF_S

  end subroutine InitializeClone


  subroutine SetBoundaryConditionsFace &
               ( FS, BoundaryCondition, iC, iD, BoundaryOption, &
                 iBoundaryOption )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    character ( * ), dimension ( 2 ), intent ( in ) :: &
      BoundaryCondition  !-- [ Inner, Outer ]
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
    character ( * ), intent ( in ), optional :: &
      BoundaryOption
    integer ( KDI ), intent ( in ), optional :: &
      iBoundaryOption

    associate &
      ( C  =>  FS % Atlas % Chart ( iC ) % Element, &
        B  =>  FS % Boundaries ( iC ) )
    if ( iD  <=  C % nDimensions ) then
      call B % SetFace &
             ( C, BoundaryCondition, iD, BoundaryOption, iBoundaryOption )
      if (       trim ( BoundaryCondition ( 1 ) )  ==  'PERIODIC'   &
           .and. trim ( BoundaryCondition ( 2 ) )  ==  'PERIODIC' ) &
      then
        C % Periodic ( iD )  =  .true.
      end if
    end if
    end associate !-- C, etc.

  end subroutine SetBoundaryConditionsFace


  subroutine SetBoundaryConditionsEdge &
               ( FS, BoundaryCondition, iC, iD, BoundaryOption, &
                 iBoundaryOption )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    character ( * ), dimension ( 4 ), intent ( in ) :: &
      BoundaryCondition  !-- [ InnerInner, OuterInner, InnerOuter, OuterOuter ]
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
    character ( * ), intent ( in ), optional :: &
      BoundaryOption
    integer ( KDI ), intent ( in ), optional :: &
      iBoundaryOption

    associate &
      ( C  =>  FS % Atlas % Chart ( iC ) % Element, &
        B  =>  FS % Boundaries ( iC ) )
    call B % SetEdge &
           ( C, BoundaryCondition, iD, BoundaryOption, iBoundaryOption )
    end associate !-- C, etc.

  end subroutine SetBoundaryConditionsEdge


  subroutine Show_FS ( FS )

    class ( FieldSetForm ), intent ( in ) :: &
      FS

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iF, &  !-- iField
      iS, &  !-- iSelected
      iV     !-- iVector
    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( FS % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FS % IGNORABILITY )

    call Show ( FS % Name, 'Name',  FS % IGNORABILITY )
    if ( associated ( FS % Primary ) ) &
      call Show ( FS % Primary % Name, 'Primary', FS % IGNORABILITY )

    call Show ( FS % Atlas % Name, 'Atlas', FS % IGNORABILITY )

    call Show ( FS % DeviceMemory, &
                'DeviceMemory', FS % IGNORABILITY ) 
    call Show ( FS % PinnedMemory, &
                'PinnedMemory', FS % IGNORABILITY ) 
    call Show ( FS % DevicesCommunicate, &
                'DevicesCommunicate', FS % IGNORABILITY ) 

    call Show ( FS % nFields, 'nFields', FS % IGNORABILITY )
    do iS  =  1, FS % nFields
      iF  =  FS % iaSelected ( iS )
      call Show ( iS,                'iSelected', FS % IGNORABILITY ) 
      call Show ( iF,                'iField',    FS % IGNORABILITY ) 
      call Show ( FS % Field ( iF ), 'Field',     FS % IGNORABILITY )
    end do !-- iF
    
    call Show ( FS % nVectors, 'nVectors', FS % IGNORABILITY )
    do iV  =  1, FS % nVectors
      call Show ( FS % Vector ( iV ),                 'Vector', &
                  FS % IGNORABILITY )
      call Show ( FS % VectorIndices ( iV ) % Value, 'VectorIndices', &
                  FS % IGNORABILITY )
    end do  !-- iV

    do iC  =  1, FS % Atlas % nCharts
      associate &
        ( C  =>  FS % Atlas % Chart ( iC ) % Element, &
          B  =>  FS % Boundaries ( iC ) )
      call Show ( C % Name, 'Chart', FS % IGNORABILITY )
      do iS  =  1, FS % nFields
        iF  =  FS % iaSelected ( iS )
        call Show ( FS % Unit ( iF, iC ),  'Unit',  FS % IGNORABILITY )
      end do !-- iS
      call B % Show ( C, FS % IGNORABILITY )
      end associate !-- C
    end do !-- iC

  end subroutine Show_FS


  function TimerGhost ( FS, Level ) result ( T )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = FS % iTimerGhost, &
               Name = trim ( FS % Name ) // '_Ghst', &
               Level = Level )

  end function TimerGhost


  subroutine Clear_FS ( FS, UseDeviceOption )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iS, &  !-- iSelected
      iF     !-- iField
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = FS % DeviceMemory
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    do iC  =  1,  FS % Atlas % nCharts
      associate ( FSV  =>  FS % Storage ( iC ) % Value )
      do iS  =  1,  FS % nFields
        iF  =  FS % iaSelected ( iS )
        call Clear ( FSV ( :, iF ), &
                     UseDeviceOption = UseDevice )
      end do !-- iS
      end associate !-- FSV
    end do !-- iC

  end subroutine Clear_FS


  subroutine Copy_FS ( FS_T, FS_S, UseDeviceOption )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_T
    class ( FieldSetForm ), intent ( in ) :: &
      FS_S
    logical ( KDL ), intent ( in ), optional :: &
       UseDeviceOption
    
    integer ( KDI ) :: &
      iC, &    !-- iChart
      iS, &    !-- iSelected
      iF_S, &  !-- iField
      iF_T
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = FS_S % DeviceMemory
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    do iC  =  1,  FS_S % Atlas % nCharts
      associate &
        ( FSV_S  =>  FS_S % Storage ( iC ) % Value, &
          FSV_T  =>  FS_T % Storage ( iC ) % Value )
      do iS  =  1,  FS_S % nFields
        iF_S  =  FS_S % iaSelected ( iS )
        iF_T  =  FS_T % iaSelected ( iS )
        call Copy ( FSV_S ( :, iF_S ), FSV_T ( :, iF_T ), &
                    UseDeviceOption = UseDevice )
      end do !-- iS
      end associate !-- FSV_S, etc.
    end do !-- iC

  end subroutine Copy_FS


  subroutine MultiplyAdd_FS ( FS_D, FS_A, FS_B, C, UseDeviceOption )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_D
    class ( FieldSetForm ), intent ( in ) :: &
      FS_A, &
      FS_B
    real ( KDR ), intent ( in ) :: &
      C
    logical ( KDL ), intent ( in ), optional :: &
       UseDeviceOption
    
    integer ( KDI ) :: &
      iC, &    !-- iChart
      iS, &    !-- iSelected
      iF_A, &  !-- iField
      iF_B, &
      iF_D
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = FS_D % DeviceMemory
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    do iC  =  1,  FS_A % Atlas % nCharts
      associate &
        ( A  =>  FS_A % Storage ( iC ) % Value, &
          B  =>  FS_B % Storage ( iC ) % Value, &
          D  =>  FS_D % Storage ( iC ) % Value )
      do iS  =  1,  FS_D % nFields
        iF_A  =  FS_A % iaSelected ( iS )
        iF_B  =  FS_B % iaSelected ( iS )
        iF_D  =  FS_D % iaSelected ( iS )
        call MultiplyAdd &
               ( A ( :, iF_A ), B ( :, iF_B ), C, D ( :, iF_D ), &
                 UseDeviceOption = UseDevice )
      end do !-- iS
      end associate !-- A, etc.
    end do !-- iC

  end subroutine MultiplyAdd_FS


  subroutine MultiplyAddInPlace_FS ( FS_A, FS_B, C, UseDeviceOption )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_A
    class ( FieldSetForm ), intent ( in ) :: &
      FS_B
    real ( KDR ), intent ( in ) :: &
      C
    logical ( KDL ), intent ( in ), optional :: &
       UseDeviceOption
    
    integer ( KDI ) :: &
      iC, &    !-- iChart
      iS, &    !-- iSelected
      iF_A, &  !-- iField
      iF_B
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = FS_A % DeviceMemory
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    do iC  =  1,  FS_A % Atlas % nCharts
      associate &
        ( A  =>  FS_A % Storage ( iC ) % Value, &
          B  =>  FS_B % Storage ( iC ) % Value )
      do iS  =  1,  FS_A % nFields
        iF_A  =  FS_A % iaSelected ( iS )
        iF_B  =  FS_B % iaSelected ( iS )
        call MultiplyAdd &
               ( A ( :, iF_A ), B ( :, iF_B ), C, &
                 UseDeviceOption = UseDevice )
      end do !-- iS
      end associate !-- A, etc.
    end do !-- iC

  end subroutine MultiplyAddInPlace_FS


  subroutine ExchangeGhostData ( FS, T_Option )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC  !-- iChart

    call FS % StartGhostExchange ( T_Option )
    call FS % FinishGhostExchange ( T_Option )

  end subroutine ExchangeGhostData


  subroutine StartGhostExchange ( FS, T_Option )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC  !-- iChart
    type ( TimerForm ), pointer :: &
      T_UH, &
      T_EG

    if ( FS % DeviceMemory .and. .not. FS % DevicesCommunicate &
         .and. present ( T_Option ) ) &
    then
      T_UH  =>  PROGRAM_HEADER % Timer &
                  ( Handle = FS % iTimerGhost_UH, &
                    Name = trim ( T_Option % Name ) // '_UpdtHst', &
                    Level = T_Option % Level + 1 )
      T_EG  =>  PROGRAM_HEADER % Timer &
                  ( Handle = FS % iTimerGhost_EG, &
                    Name = trim ( T_Option % Name ) // '_Exchng', &
                    Level = T_Option % Level + 1 )
    else
      T_UH  =>  null ( )
      T_EG  =>  null ( )
    end if

    if ( .not. FS % DevicesCommunicate ) then
      if ( associated ( T_UH ) ) call T_UH % Start ( ) 
      call FS % UpdateHost ( )
      if ( associated ( T_UH ) ) call T_UH % Stop ( ) 
    end if

    if ( associated ( T_EG ) ) call T_EG % Start ( ) 
    do iC  =  1,  FS % Atlas % nCharts
      associate &
        ( GE  =>  FS % GhostExchange ( iC ), &
           S  =>  FS % Storage ( iC ), &
           C  =>  FS % Atlas % Chart ( iC ) % Element )
      call GE % StartExchange ( C, S, FS % DevicesCommunicate )
      end associate !-- GE, etc.
    end do !-- iC
    if ( associated ( T_EG ) ) call T_EG % Stop ( ) 

  end subroutine StartGhostExchange


  subroutine FinishGhostExchange ( FS, T_Option )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
    type ( TimerForm ), pointer :: &
      T_EG, &
      T_UD

    if ( FS % DeviceMemory .and. .not. FS % DevicesCommunicate &
         .and. present ( T_Option ) ) &
    then
      T_EG  =>  PROGRAM_HEADER % Timer &
                  ( Handle = FS % iTimerGhost_EG, &
                    Name = trim ( T_Option % Name ) // '_Exchng', &
                    Level = T_Option % Level + 1 )
      T_UD  =>  PROGRAM_HEADER % Timer &
                  ( Handle = FS % iTimerGhost_UD, &
                    Name = trim ( T_Option % Name ) // '_UpdtDvc', &
                    Level = T_Option % Level + 1 )
    else
      T_EG  =>  null ( )
      T_UD  =>  null ( )
    end if

    if ( associated ( T_EG ) ) call T_EG % Start ( ) 
    do iC  =  1,  FS % Atlas % nCharts
      associate &
        (  B  =>  FS % Boundaries ( iC ), &
          GE  =>  FS % GhostExchange ( iC ), &
           S  =>  FS % Storage ( iC ), &
           C  =>  FS % Atlas % Chart ( iC ) % Element )
      call GE % FinishExchange ( S, C, C % Periodic, FS % DevicesCommunicate )
      end associate !-- GE, etc.
    end do !-- iC
    if ( associated ( T_EG ) ) call T_EG % Stop ( ) 

    if ( .not. FS % DevicesCommunicate ) then
      if ( associated ( T_UD ) ) call T_UD % Start ( ) 
      call FS % UpdateDevice ( )
      if ( associated ( T_UD ) ) call T_UD % Stop ( ) 
    end if

  end subroutine FinishGhostExchange


  subroutine ApplyBoundaryConditions ( FS )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1,  FS % Atlas % nCharts
      associate &
        ( B  =>  FS % Boundaries ( iC ), &
          S  =>  FS % Storage ( iC ), &
          C  =>  FS % Atlas % Chart ( iC ) % Element )
      call B % Apply ( S, C )
      end associate !-- B, etc.
    end do !-- iC

  end subroutine ApplyBoundaryConditions


  subroutine UpdateDevice_FS ( FS )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1,  FS % Atlas % nCharts
      call FS % Storage ( iC ) % UpdateDevice ( )
    end do !-- iC

  end subroutine UpdateDevice_FS


  subroutine UpdateHost_FS ( FS )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1,  FS % Atlas % nCharts
      call FS % Storage ( iC ) % UpdateHost ( )
    end do !-- iC

  end subroutine UpdateHost_FS


  impure elemental subroutine Finalize ( FS )

    type ( FieldSetForm ), intent ( inout ) :: &
      FS

    nullify ( FS % Primary )
    nullify ( FS % Atlas )
    nullify ( FS % Storage_GS )

    if ( allocated ( FS % Boundaries ) ) &
      deallocate ( FS % Boundaries )
    if ( allocated ( FS % GhostExchange ) ) &
      deallocate ( FS % GhostExchange )
    if ( allocated ( FS % Storage ) ) &
      deallocate ( FS % Storage )
    if ( allocated ( FS % Vector ) ) &
      deallocate ( FS % Vector )
    if ( allocated ( FS % Field ) ) &
      deallocate ( FS % Field )
    if ( allocated ( FS % Unit ) ) &
      deallocate ( FS % Unit )
    if ( allocated ( FS % VectorIndices ) ) &
      deallocate ( FS % VectorIndices )
    if ( allocated ( FS % iaSelected ) ) &
      deallocate ( FS % iaSelected )

    call Show ( 'Finalizing ' // trim ( FS % Type ), FS % IGNORABILITY )
    call Show ( FS % Name, 'Name', FS % IGNORABILITY )
   
  end subroutine Finalize


  impure elemental subroutine Finalize_E ( FSE )
    
    type ( FieldSetElement ), intent ( inout ) :: &
      FSE

    if ( allocated ( FSE % Element ) ) &
      deallocate ( FSE % Element )

  end subroutine Finalize_E


end module FieldSet_Form
