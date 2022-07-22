module Geometry_F__Form

  !-- Geometry_Flat__Form

  use Basics
  use Manifolds
  use FieldSets

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_F  = 25, &
      N_VECTORS_F =  0

  type, public, extends ( FieldSetForm ) :: Geometry_F_Form
    integer ( KDI ) :: &
      N_FIELDS_F = N_FIELDS_F, &
      N_VECTORS_F = N_VECTORS_F
    integer ( KDI ) :: &
      !-- Coordinate fields
      EDGE_I_U_1 = 0, &
      EDGE_I_U_2 = 0, &
      EDGE_I_U_3 = 0, &
      WIDTH_U_1  = 0, &
      WIDTH_U_2  = 0, &
      WIDTH_U_3  = 0, &
      CENTER_U_1 = 0, &
      CENTER_U_2 = 0, &
      CENTER_U_3 = 0
    integer ( KDI ) :: &
      !-- Finite volume fields
      AREA_I_D_1   = 0, &
      AREA_I_D_2   = 0, &
      AREA_I_D_3   = 0, &
      VOLUME       = 0, &
      AVERAGE_1_U_1 = 0, &
      AVERAGE_1_U_2 = 0, &
      AVERAGE_1_U_3 = 0, &
      AVERAGE_2_U_1 = 0, &
      AVERAGE_2_U_2 = 0, &
      AVERAGE_2_U_3 = 0
    integer ( KDI ) :: &
      !-- Flat metric fields
      METRIC_F_DD_11 = 0, &
      METRIC_F_DD_22 = 0, &
      METRIC_F_DD_33 = 0, &
      METRIC_F_UU_11 = 0, &
      METRIC_F_UU_22 = 0, &
      METRIC_F_UU_33 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      EDGE_I_U, &
      WIDTH_U, &
      CENTER_U, &
      AREA_I_D, &
      AVERAGE_1_U, &
      AVERAGE_2_U, &
      METRIC_F_DD, &
      METRIC_F_UU
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass ( G ) :: &
      SetStream
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      ComputeReconstruction
    final :: &
      Finalize
  end type Geometry_F_Form

    private :: &
      SetUnits, &
      SetCoordinates, &
      ComputeFromCoordinates

      private :: &
        Compute_FV_R_Kernel, &
        Compute_FV_C_Kernel, &
        Compute_FV_S_Kernel, &
        Compute_M_R_Kernel, &
        Compute_M_C_Kernel, &
        Compute_M_S_Kernel

    interface
      
      module subroutine Compute_FV_R_Kernel &
               ( W_1, W_2, W_3, E_I_1, E_I_2, E_I_3, nD, &
                 A_I_1, A_I_2, A_I_3, V, &
                 AV_1_1, AV_1_2, AV_1_3, AV_2_1, AV_2_2, AV_2_3 )
        !-- Compute_FiniteVolume_Rectangular_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3, &
          E_I_1, E_I_2, E_I_3
        integer ( KDI ), intent ( in ) :: &
          nD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          A_I_1, A_I_2, A_I_3, &
          V, &
          AV_1_1, AV_1_2, AV_1_3, &
          AV_2_1, AV_2_2, AV_2_3
      end subroutine Compute_FV_R_Kernel

      module subroutine Compute_FV_C_Kernel &
               ( W_1, W_2, W_3, E_I_1, E_I_2, E_I_3, nD, &
                 A_I_1, A_I_2, A_I_3, V, &
                 AV_1_1, AV_1_2, AV_1_3, AV_2_1, AV_2_2, AV_2_3 )
        !-- Compute_FiniteVolume_Cylindrical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3, &
          E_I_1, E_I_2, E_I_3
        integer ( KDI ), intent ( in ) :: &
          nD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          A_I_1, A_I_2, A_I_3, &
          V, &
          AV_1_1, AV_1_2, AV_1_3, &
          AV_2_1, AV_2_2, AV_2_3
      end subroutine Compute_FV_C_Kernel

      module subroutine Compute_FV_S_Kernel &
               ( W_1, W_2, W_3, E_I_1, E_I_2, E_I_3, nD, &
                 A_I_1, A_I_2, A_I_3, V, &
                 AV_1_1, AV_1_2, AV_1_3, AV_2_1, AV_2_2, AV_2_3 )
        !-- Compute_FiniteVolume_Spherical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3, &
          E_I_1, E_I_2, E_I_3
        integer ( KDI ), intent ( in ) :: &
          nD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          A_I_1, A_I_2, A_I_3, &
          V, &
          AV_1_1, AV_1_2, AV_1_3, &
          AV_2_1, AV_2_2, AV_2_3
      end subroutine Compute_FV_S_Kernel

      module subroutine Compute_M_R_Kernel &
               ( M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 UseDeviceOption )
        !-- Compute_Metric_Rectangular_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_11, M_UU_22, M_UU_33
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_M_R_Kernel

      module subroutine Compute_M_C_Kernel &
               ( RP, nD, M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 UseDeviceOption )
        !-- Compute_Metric_Cylindrical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          RP
        integer ( KDI ), intent ( in ) :: &
          nD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_11, M_UU_22, M_UU_33
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_M_C_Kernel

      module subroutine Compute_M_S_Kernel &
               ( R, Th, nD, M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, &
                 M_UU_33, UseDeviceOption )
        !-- Compute_Metric_Spherical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          R, Th
        integer ( KDI ), intent ( in ) :: &
          nD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_11, M_UU_22, M_UU_33
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_M_S_Kernel

    end interface


contains


  subroutine InitializeAllocate_FS &
               ( FS, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, AssociateFieldsOption, &
                 UnitOption, VectorIndicesOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( Geometry_F_Form ), intent ( inout ), target :: &
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
      nFields
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a Geometry_F' 
    
    Name  =  'Geometry'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    FS % EDGE_I_U_1      =   1
    FS % EDGE_I_U_2      =   2
    FS % EDGE_I_U_3      =   3
    FS % WIDTH_U_1       =   4
    FS % WIDTH_U_2       =   5
    FS % WIDTH_U_3       =   6
    FS % CENTER_U_1      =   7
    FS % CENTER_U_2      =   8
    FS % CENTER_U_3      =   9
    FS % AREA_I_D_1      =  10
    FS % AREA_I_D_2      =  11
    FS % AREA_I_D_3      =  12
    FS % VOLUME          =  13
    FS % AVERAGE_1_U_1   =  14
    FS % AVERAGE_1_U_2   =  15
    FS % AVERAGE_1_U_3   =  16
    FS % AVERAGE_2_U_1   =  17
    FS % AVERAGE_2_U_2   =  18
    FS % AVERAGE_2_U_3   =  19
    FS % METRIC_F_DD_11  =  20
    FS % METRIC_F_DD_22  =  21
    FS % METRIC_F_DD_33  =  22
    FS % METRIC_F_UU_11  =  23
    FS % METRIC_F_UU_22  =  24
    FS % METRIC_F_UU_33  =  25

    nFields  =  FS % N_FIELDS_F
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    FS % EDGE_I_U  &
      =  [ FS % EDGE_I_U_1, FS % EDGE_I_U_2, FS % EDGE_I_U_3 ]
    FS % WIDTH_U  &
      =  [ FS % WIDTH_U_1, FS % WIDTH_U_2, FS % WIDTH_U_3  ]
    FS % CENTER_U  &
      =  [ FS % CENTER_U_1, FS % CENTER_U_2, FS % CENTER_U_3 ]
    FS % AREA_I_D  &
      =  [ FS % AREA_I_D_1, FS % AREA_I_D_2, FS % AREA_I_D_3 ]
    FS % AVERAGE_1_U  &
      =  [ FS % AVERAGE_1_U_1, FS % AVERAGE_1_U_2, FS % AVERAGE_1_U_3 ]
    FS % AVERAGE_2_U  &
      =  [ FS % AVERAGE_2_U_1, FS % AVERAGE_2_U_2, FS % AVERAGE_2_U_3 ]
    FS % METRIC_F_DD  &
      =  [ FS % METRIC_F_DD_11, FS % METRIC_F_DD_22, FS % METRIC_F_DD_33 ]
    FS % METRIC_F_UU  &
      =  [ FS % METRIC_F_UU_11, FS % METRIC_F_UU_22, FS % METRIC_F_UU_33 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( 1 : FS % N_FIELDS_F ) &
      = [ 'Edge_I_U_1    ', &
          'Edge_I_U_2    ', &
          'Edge_I_U_3    ', &
          'Width_U_1     ', &
          'Width_U_2     ', &
          'Width_U_3     ', &
          'Center_U_1    ', &
          'Center_U_2    ', &
          'Center_U_3    ', &
          'Area_I_D_1    ', &
          'Area_I_D_2    ', &
          'Area_I_D_3    ', &
          'Volume        ', &
          'Average_1_U_1 ', &
          'Average_1_U_2 ', &
          'Average_1_U_3 ', &
          'Average_2_U_1 ', &
          'Average_2_U_2 ', &
          'Average_2_U_3 ', &
          'Metric_F_DD_11', &
          'Metric_F_DD_22', &
          'Metric_F_DD_33', &
          'Metric_F_UU_11', &
          'Metric_F_UU_22', &
          'Metric_F_UU_33' ]

    !-- Units

    if ( present ( UnitOption ) ) then
      allocate ( Unit, source = UnitOption )
    else
      allocate ( Unit ( nFields, A % nCharts ) )
    end if !-- FieldOption

    call SetUnits ( Unit, FS, A )

    !-- FieldSet

    call FS % FieldSetForm % Initialize &
           ( A, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             AssociateFieldsOption = AssociateFieldsOption, &
             UnitOption = Unit, &
             VectorIndicesOption = VectorIndicesOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    call FS % Compute ( )

  end subroutine InitializeAllocate_FS


  subroutine SetStream ( S, G, iaAdditionalOption )

    class ( StreamForm ), intent ( inout ) :: &
      S
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaAdditionalOption

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected

    if ( present ( iaAdditionalOption ) ) then
      allocate ( iaSelected ( 6 + size ( iaAdditionalOption ) ) )
      iaSelected ( 7 : )  =  iaAdditionalOption
    else
      allocate ( iaSelected ( 6 ) )
    end if

    iaSelected ( 1 : 6 )  &
      =  [ G % CENTER_U_1,     G % CENTER_U_2,     G % CENTER_U_3, &
           G % METRIC_F_DD_11, G % METRIC_F_DD_22, G % METRIC_F_DD_33 ]

    call S % AddFieldSet ( G, iaSelectedOption = iaSelected )

  end subroutine SetStream


  subroutine Compute ( G )

    class ( Geometry_F_Form ), intent ( inout ) :: &
      G

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    associate ( A  =>  G % Atlas )
    do iC  =  1,  A % nCharts
      associate ( C  =>  A % Chart ( iC ) % Element )

      do iD  =  1,  C % nDimensions
        call SetCoordinates ( G, iC, iD )
      end do !-- iD
      call ComputeFromCoordinates ( G, iC )

      end associate !-- C
    end do !-- iC
    end associate !-- A

  end subroutine Compute


  subroutine ComputeReconstruction ( G, M_I, iC, iD )

    class ( Geometry_F_Form ), intent ( inout ) :: &
      G
    type ( FieldSetForm ), intent ( inout ) :: &
      M_I
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions

    associate &
      ( A  =>  G % Atlas )
    associate &
      ( C        =>  A % Chart ( iC ) % Element, &
        GV       =>  G % Storage ( iC ) % Value, &
        M_DD_11  =>  M_I % Storage ( iC ) % Value ( :, 1 ), &
        M_DD_22  =>  M_I % Storage ( iC ) % Value ( :, 2 ), &
        M_DD_33  =>  M_I % Storage ( iC ) % Value ( :, 3 ), &
        M_UU_11  =>  M_I % Storage ( iC ) % Value ( :, 4 ), &
        M_UU_22  =>  M_I % Storage ( iC ) % Value ( :, 5 ), &
        M_UU_33  =>  M_I % Storage ( iC ) % Value ( :, 6 ) )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call Compute_M_R_Kernel &
             ( M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
               UseDeviceOption = M_I % DeviceMemory )
    case ( 'CYLINDRICAL' )
      select case ( iD )
      case ( 1 )
        call Compute_M_C_Kernel &
               ( GV ( :, G % EDGE_I_U_1 ), &
                 C % nDimensions, &
                 M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 UseDeviceOption = M_I % DeviceMemory )
      case ( 2, 3 )
        call Compute_M_C_Kernel &
               ( GV ( :, G % CENTER_U_1 ), &
                 C % nDimensions, &
                 M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 UseDeviceOption = M_I % DeviceMemory )
      end select !-- iD
    case ( 'SPHERICAL' )
      select case ( iD )
      case ( 1 )
        call Compute_M_S_Kernel &
               ( GV ( :, G % EDGE_I_U_1 ), &
                 GV ( :, G % CENTER_U_2 ), &
                 C % nDimensions, &
                 M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 UseDeviceOption = M_I % DeviceMemory )
      case ( 2 )
        call Compute_M_S_Kernel &
               ( GV ( :, G % CENTER_U_1 ), &
                 GV ( :, G % EDGE_I_U_2 ), &
                 C % nDimensions, &
                 M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 UseDeviceOption = M_I % DeviceMemory )
      case ( 3 )
        call Compute_M_S_Kernel &
               ( GV ( :, G % CENTER_U_1 ), &
                 GV ( :, G % CENTER_U_2 ), &
                 C % nDimensions, &
                 M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 UseDeviceOption = M_I % DeviceMemory )
      end select !-- iD
    case default
      call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', &
                  CONSOLE % ERROR )
      call Show ( 'Geometry_F__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeReconstruction', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    end associate !-- C, etc.
    end associate !-- A

  end subroutine ComputeReconstruction


  impure elemental subroutine Finalize ( G )

    type ( Geometry_F_Form ), intent ( inout ) :: &
      G

  end subroutine Finalize

  
  subroutine SetUnits ( FieldUnit, G, A )

    type ( QuantityForm ), dimension ( :, : ), intent ( inout ) :: &
      FieldUnit
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( Atlas_H_Form ), intent ( in ) :: &
      A

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1,  A % nCharts

      associate &
        ( C  =>  A % Chart ( iC ) % Element )
      associate &
        ( CoordinateUnit    =>  C % CoordinateUnit, &
          CoordinateSystem  =>  C % CoordinateSystem )

      FieldUnit ( G % EDGE_I_U_1 : G % EDGE_I_U_3, iC ) &
        = CoordinateUnit
      FieldUnit ( G % WIDTH_U_1 : G % WIDTH_U_3, iC ) &
        = CoordinateUnit
      FieldUnit ( G % CENTER_U_1 : G % CENTER_U_3, iC ) &
        = CoordinateUnit

      select case ( trim ( CoordinateSystem ) )
      case ( 'RECTANGULAR' )
        FieldUnit ( G % VOLUME, iC )  &
          =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )  &
             *  CoordinateUnit ( 3 )
        FieldUnit ( G % AREA_I_D_1, iC )  &
          =  CoordinateUnit ( 2 )  *  CoordinateUnit ( 3 )
        FieldUnit ( G % AREA_I_D_2, iC )  &
          =  CoordinateUnit ( 3 )  *  CoordinateUnit ( 1 )
        FieldUnit ( G % AREA_I_D_3, iC )  &
          =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
        FieldUnit ( G % METRIC_F_DD_11, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_DD_22, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_DD_33, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_UU_11, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_UU_22, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_UU_33, iC ) = UNIT % IDENTITY
      case ( 'CYLINDRICAL' )
        FieldUnit ( G % VOLUME, iC )  &
          =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
        FieldUnit ( G % AREA_I_D_1, iC )  &
          =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
        FieldUnit ( G % AREA_I_D_2, iC )  &
          =  CoordinateUnit ( 1 ) ** 2
        FieldUnit ( G % AREA_I_D_3, iC )  &
          =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
        FieldUnit ( G % METRIC_F_DD_11, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_DD_22, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_DD_33, iC ) = CoordinateUnit ( 1 ) ** (  2 )
        FieldUnit ( G % METRIC_F_UU_11, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_UU_22, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_UU_33, iC ) = CoordinateUnit ( 1 ) ** ( -2 )
      case ( 'SPHERICAL' )
        FieldUnit ( G % VOLUME, iC )  &
          = CoordinateUnit ( 1 ) ** 3
        FieldUnit ( G % AREA_I_D_1, iC )  &
          =  CoordinateUnit ( 1 ) ** 2
        FieldUnit ( G % AREA_I_D_2, iC )  &
          =  CoordinateUnit ( 1 ) ** 3
        FieldUnit ( G % AREA_I_D_3, iC )  &
          =  CoordinateUnit ( 1 ) ** 3
        FieldUnit ( G % METRIC_F_DD_11, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_DD_22, iC ) = CoordinateUnit ( 1 ) ** (  2 )
        FieldUnit ( G % METRIC_F_DD_33, iC ) = CoordinateUnit ( 1 ) ** (  2 )
        FieldUnit ( G % METRIC_F_UU_11, iC ) = UNIT % IDENTITY
        FieldUnit ( G % METRIC_F_UU_22, iC ) = CoordinateUnit ( 1 ) ** ( -2 )
        FieldUnit ( G % METRIC_F_UU_33, iC ) = CoordinateUnit ( 1 ) ** ( -2 )
      end select !-- CoordinateSystem

      end associate !-- CoordinateUnit
      end associate !-- C

    end do !-- iC

  end subroutine SetUnits

  
  subroutine SetCoordinates ( G, iCt, iD )

    class ( Geometry_F_Form ), intent ( inout ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      iCt, &  !-- iChart
      iD      !-- iDimension

    integer ( KDI ) :: &
      iaF, iaL, &  !-- iaFirst, iaLast
      iC, &        !-- iCell
      oC           !-- oCell
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Edge_I_3D, &
      Width_3D, &
      Center_3D

    associate &
      (  A  =>  G % Atlas, &
        GV  =>  G % Storage ( iCt ) % Value )
    select type ( C  =>  A % Chart ( iCt ) % Element )
      class is ( Chart_GS_Form )

    iaF  =  1  -  C % nGhostLayers ( iD ) 
    if ( C % Distributed ) then
      iaL  =  C % nCellsBrick ( iD )  +  C % nGhostLayers ( iD )
       oC  =  ( C % iaBrick ( iD )  -  1 )  *  C % nCellsBrick ( iD )
    else
      iaL  =  C % nCells ( iD )  +  C % nGhostLayers ( iD )
       oC  =  0
    end if

    call C % SetFieldPointer &
           ( GV ( :, G % EDGE_I_U ( iD ) ), Edge_I_3D )
    call C % SetFieldPointer &
           ( GV ( :, G % WIDTH_U ( iD ) ),  Width_3D )
    call C % SetFieldPointer &
           ( GV ( :, G % CENTER_U ( iD ) ), Center_3D )

    associate &
      (   Edge_1D  =>  C %   Edge ( iD ) % Value, &
         Width_1D  =>  C %  Width ( iD ) % Value, &
        Center_1D  =>  C % Center ( iD ) % Value )
    do iC  =  iaF, iaL
      select case ( iD )
      case ( 1 )
        Edge_I_3D ( iC, :, : )  =    Edge_1D ( oC + iC )
         Width_3D ( iC, :, : )  =   Width_1D ( oC + iC )
        Center_3D ( iC, :, : )  =  Center_1D ( oC + iC )
      case ( 2 )
        Edge_I_3D ( :, iC, : )  =    Edge_1D ( oC + iC )
         Width_3D ( :, iC, : )  =   Width_1D ( oC + iC )
        Center_3D ( :, iC, : )  =  Center_1D ( oC + iC )
      case ( 3 )
        Edge_I_3D ( :, :, iC )  =    Edge_1D ( oC + iC )
         Width_3D ( :, :, iC )  =   Width_1D ( oC + iC )
        Center_3D ( :, :, iC )  =  Center_1D ( oC + iC )
      end select !-- iD
    end do !-- iC
    end associate !-- Edge_1D, etc.

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Geometry_F__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCoordinates', 'subroutine', CONSOLE % ERROR )
    end select !-- C  
    end associate !-- A, etc.

  end subroutine SetCoordinates


  subroutine ComputeFromCoordinates ( G, iC )

    !-- Assumes coordinate fields are set

    class ( Geometry_F_Form ), intent ( inout ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      iC  !-- iChart

    associate &
      ( A  =>  G % Atlas )
    associate &
      (  C  =>  A % Chart ( iC ) % Element, &
        GV  =>  G % Storage ( iC ) % Value )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call Compute_FV_R_Kernel &
             ( GV ( :, G % WIDTH_U_1 ), &
               GV ( :, G % WIDTH_U_2 ), &
               GV ( :, G % WIDTH_U_3 ), &
               GV ( :, G % EDGE_I_U_1 ), &
               GV ( :, G % EDGE_I_U_2 ), &
               GV ( :, G % EDGE_I_U_3 ), &
               C % nDimensions, &
               GV ( :, G % AREA_I_D_1 ), &
               GV ( :, G % AREA_I_D_2 ), &
               GV ( :, G % AREA_I_D_3 ), &
               GV ( :, G % VOLUME ), &
               GV ( :, G % AVERAGE_1_U_1 ), &
               GV ( :, G % AVERAGE_1_U_2 ), &
               GV ( :, G % AVERAGE_1_U_3 ), &
               GV ( :, G % AVERAGE_2_U_1 ), &
               GV ( :, G % AVERAGE_2_U_2 ), &
               GV ( :, G % AVERAGE_2_U_3 ) )
      call Compute_M_R_Kernel &
             ( GV ( :, G % METRIC_F_DD_11 ), &
               GV ( :, G % METRIC_F_DD_22 ), &
               GV ( :, G % METRIC_F_DD_33 ), &
               GV ( :, G % METRIC_F_UU_11 ), &
               GV ( :, G % METRIC_F_UU_22 ), &
               GV ( :, G % METRIC_F_UU_33 ) )
    case ( 'CYLINDRICAL' )
      call Compute_FV_C_Kernel &
             ( GV ( :, G % WIDTH_U_1 ), &
               GV ( :, G % WIDTH_U_2 ), &
               GV ( :, G % WIDTH_U_3 ), &
               GV ( :, G % EDGE_I_U_1 ), &
               GV ( :, G % EDGE_I_U_2 ), &
               GV ( :, G % EDGE_I_U_3 ), &
               C % nDimensions, &
               GV ( :, G % AREA_I_D_1 ), &
               GV ( :, G % AREA_I_D_2 ), &
               GV ( :, G % AREA_I_D_3 ), &
               GV ( :, G % VOLUME ), &
               GV ( :, G % AVERAGE_1_U_1 ), &
               GV ( :, G % AVERAGE_1_U_2 ), &
               GV ( :, G % AVERAGE_1_U_3 ), &
               GV ( :, G % AVERAGE_2_U_1 ), &
               GV ( :, G % AVERAGE_2_U_2 ), &
               GV ( :, G % AVERAGE_2_U_3 ) )
      call Compute_M_C_Kernel &
             ( GV ( :, G % CENTER_U_1 ), &
               C % nDimensions, &
               GV ( :, G % METRIC_F_DD_11 ), &
               GV ( :, G % METRIC_F_DD_22 ), &
               GV ( :, G % METRIC_F_DD_33 ), &
               GV ( :, G % METRIC_F_UU_11 ), &
               GV ( :, G % METRIC_F_UU_22 ), &
               GV ( :, G % METRIC_F_UU_33 ) )
    case ( 'SPHERICAL' )
      call Compute_FV_S_Kernel &
             ( GV ( :, G % WIDTH_U_1 ), &
               GV ( :, G % WIDTH_U_2 ), &
               GV ( :, G % WIDTH_U_3 ), &
               GV ( :, G % EDGE_I_U_1 ), &
               GV ( :, G % EDGE_I_U_2 ), &
               GV ( :, G % EDGE_I_U_3 ), &
               C % nDimensions, &
               GV ( :, G % AREA_I_D_1 ), &
               GV ( :, G % AREA_I_D_2 ), &
               GV ( :, G % AREA_I_D_3 ), &
               GV ( :, G % VOLUME ), &
               GV ( :, G % AVERAGE_1_U_1 ), &
               GV ( :, G % AVERAGE_1_U_2 ), &
               GV ( :, G % AVERAGE_1_U_3 ), &
               GV ( :, G % AVERAGE_2_U_1 ), &
               GV ( :, G % AVERAGE_2_U_2 ), &
               GV ( :, G % AVERAGE_2_U_3 ) )
      call Compute_M_S_Kernel &
             ( GV ( :, G % CENTER_U_1 ), &
               GV ( :, G % CENTER_U_2 ), &
               C % nDimensions, &
               GV ( :, G % METRIC_F_DD_11 ), &
               GV ( :, G % METRIC_F_DD_22 ), &
               GV ( :, G % METRIC_F_DD_33 ), &
               GV ( :, G % METRIC_F_UU_11 ), &
               GV ( :, G % METRIC_F_UU_22 ), &
               GV ( :, G % METRIC_F_UU_33 ) )
    case default
      call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', &
                  CONSOLE % ERROR )
      call Show ( 'Geometry_F__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeFromCoordinates', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select
  
    end associate !-- C, etc.
    end associate !-- A

  end subroutine ComputeFromCoordinates


end module Geometry_F__Form
