module Geometry_F_C__Form

  !-- Geometry_Flat_Chart_Form

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_F  = 19, &
      N_VECTORS_F =  0

  type, public, extends ( FieldSet_C_Form ) :: Geometry_F_C_Form
    integer ( KDI ) :: &
      N_FIELDS_F = N_FIELDS_F, &
      N_VECTORS_F = N_VECTORS_F
    integer ( KDI ) :: &
      !-- Coordinate fields
      EDGE_I_U_1   = 0, &
      EDGE_I_U_2   = 0, &
      EDGE_I_U_3   = 0, &
      WIDTH_U_1    = 0, &
      WIDTH_U_2    = 0, &
      WIDTH_U_3    = 0, &
      CENTER_U_1   = 0, &
      CENTER_U_2   = 0, &
      CENTER_U_3   = 0
    integer ( KDI ) :: &
      !-- Finite volume fields
      AREA_I_D_1   = 0, &
      AREA_I_D_2   = 0, &
      AREA_I_D_3   = 0, &
      VOLUME       = 0
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
      AREA_I_D
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass ( GC ) :: &
      SetStream
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Geometry_F_C_Form

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
               ( W_1, W_2, W_3, nD, A_I_1, A_I_2, A_I_3, V )
        !-- Compute_FiniteVolume_Rectangular_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3
        integer ( KDI ), intent ( in ) :: &
          nD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          A_I_1, A_I_2, A_I_3, &
          V
      end subroutine Compute_FV_R_Kernel

      module subroutine Compute_FV_C_Kernel &
               ( W_1, W_2, W_3, E_I_1, nD, A_I_1, A_I_2, A_I_3, V )
        !-- Compute_FiniteVolume_Cylindrical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3, &
          E_I_1
        integer ( KDI ), intent ( in ) :: &
          nD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          A_I_1, A_I_2, A_I_3, &
          V
      end subroutine Compute_FV_C_Kernel

      module subroutine Compute_FV_S_Kernel &
               ( W_1, W_2, W_3, E_I_1, E_I_2, nD, A_I_1, A_I_2, A_I_3, V )
        !-- Compute_FiniteVolume_Spherical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3, &
          E_I_1, E_I_2
        integer ( KDI ), intent ( in ) :: &
          nD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          A_I_1, A_I_2, A_I_3, &
          V
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
               ( FSC, C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Geometry_F_C_Form ), intent ( inout ) :: &
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
      nFields
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a Geometry_F_C' 
    
    Name  =  'Geometry'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    FSC % EDGE_I_U_1      =   1
    FSC % EDGE_I_U_2      =   2
    FSC % EDGE_I_U_3      =   3
    FSC % WIDTH_U_1       =   4
    FSC % WIDTH_U_2       =   5
    FSC % WIDTH_U_3       =   6
    FSC % CENTER_U_1      =   7
    FSC % CENTER_U_2      =   8
    FSC % CENTER_U_3      =   9
    FSC % AREA_I_D_1      =  10
    FSC % AREA_I_D_2      =  11
    FSC % AREA_I_D_3      =  12
    FSC % VOLUME          =  13
    FSC % METRIC_F_DD_11  =  14
    FSC % METRIC_F_DD_22  =  15
    FSC % METRIC_F_DD_33  =  16
    FSC % METRIC_F_UU_11  =  17
    FSC % METRIC_F_UU_22  =  18
    FSC % METRIC_F_UU_33  =  19

    nFields  =  FSC % N_FIELDS_F
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    FSC % EDGE_I_U  =  [ FSC % EDGE_I_U_1, FSC % EDGE_I_U_2, FSC % EDGE_I_U_3 ]
    FSC % WIDTH_U   =  [ FSC % WIDTH_U_1,  FSC % WIDTH_U_2,  FSC % WIDTH_U_3  ]
    FSC % CENTER_U  =  [ FSC % CENTER_U_1, FSC % CENTER_U_2, FSC % CENTER_U_3 ]
    FSC % AREA_I_D  =  [ FSC % AREA_I_D_1, FSC % AREA_I_D_2, FSC % AREA_I_D_3 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( 1 : FSC % N_FIELDS_F ) &
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
      allocate ( Unit ( nFields ) )
    end if !-- FieldOption

    call SetUnits ( Unit, FSC, C )

    !-- FieldSet

    call FSC % FieldSet_C_Form % Initialize &
           ( C, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             UnitOption = Unit, &
             VectorIndicesOption = VectorIndicesOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    call FSC % Compute ( )

  end subroutine InitializeAllocate_FS


  subroutine SetStream ( SC, GC, iaAdditionalOption )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    class ( Geometry_F_C_Form ), intent ( in ) :: &
      GC
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
      =  [ GC % CENTER_U_1,     GC % CENTER_U_2,     GC % CENTER_U_3, &
           GC % METRIC_F_DD_11, GC % METRIC_F_DD_22, GC % METRIC_F_DD_33 ]

    call SC % AddFieldSet ( GC, iaSelectedOption = iaSelected )

  end subroutine SetStream


  subroutine Compute ( GC )

    class ( Geometry_F_C_Form ), intent ( inout ) :: &
      GC

    integer ( KDI ) :: &
      iD  !-- iDimension

    associate ( nD  =>  GC % Chart % nDimensions )
    do iD = 1, nD
      call SetCoordinates ( GC, iD )
    end do !-- iD
    call ComputeFromCoordinates ( GC % Storage_FSC % Storage, GC )
    end associate !-- nD

  end subroutine Compute


  impure elemental subroutine Finalize ( GC )

    type ( Geometry_F_C_Form ), intent ( inout ) :: &
      GC

  end subroutine Finalize

  
  subroutine SetUnits ( FieldUnit, GC, C )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      FieldUnit
    class ( Geometry_F_C_Form ), intent ( in ) :: &
      GC
    class ( Chart_H_Form ), intent ( in ) :: &
      C

    associate &
      ( CoordinateUnit    =>  C % CoordinateUnit, &
        CoordinateSystem  =>  C % CoordinateSystem )

    FieldUnit ( GC % EDGE_I_U_1 : GC % EDGE_I_U_3 ) &
      = CoordinateUnit
    FieldUnit ( GC % WIDTH_U_1 : GC % WIDTH_U_3 ) &
      = CoordinateUnit
    FieldUnit ( GC % CENTER_U_1 : GC % CENTER_U_3 ) &
      = CoordinateUnit

    select case ( trim ( CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      FieldUnit ( GC % VOLUME )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )  &
           *  CoordinateUnit ( 3 )
      FieldUnit ( GC % AREA_I_D_1 )  &
        =  CoordinateUnit ( 2 )  *  CoordinateUnit ( 3 )
      FieldUnit ( GC % AREA_I_D_2 )  &
        =  CoordinateUnit ( 3 )  *  CoordinateUnit ( 1 )
      FieldUnit ( GC % AREA_I_D_3 )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
      FieldUnit ( GC % METRIC_F_DD_11 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_DD_22 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_DD_33 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_UU_11 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_UU_22 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_UU_33 ) = UNIT % IDENTITY
    case ( 'CYLINDRICAL' )
      FieldUnit ( GC % VOLUME )  &
        =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
      FieldUnit ( GC % AREA_I_D_1 )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
      FieldUnit ( GC % AREA_I_D_2 )  &
        =  CoordinateUnit ( 1 ) ** 2
      FieldUnit ( GC % AREA_I_D_3 )  &
        =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
      FieldUnit ( GC % METRIC_F_DD_11 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_DD_22 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_DD_33 ) = CoordinateUnit ( 1 ) ** (  2 )
      FieldUnit ( GC % METRIC_F_UU_11 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_UU_22 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_UU_33 ) = CoordinateUnit ( 1 ) ** ( -2 )
    case ( 'SPHERICAL' )
      FieldUnit ( GC % VOLUME )  &
        = CoordinateUnit ( 1 ) ** 3
      FieldUnit ( GC % AREA_I_D_1 )  &
        =  CoordinateUnit ( 1 ) ** 2
      FieldUnit ( GC % AREA_I_D_2 )  &
        =  CoordinateUnit ( 1 ) ** 3
      FieldUnit ( GC % AREA_I_D_3 )  &
        =  CoordinateUnit ( 1 ) ** 3
      FieldUnit ( GC % METRIC_F_DD_11 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_DD_22 ) = CoordinateUnit ( 1 ) ** (  2 )
      FieldUnit ( GC % METRIC_F_DD_33 ) = CoordinateUnit ( 1 ) ** (  2 )
      FieldUnit ( GC % METRIC_F_UU_11 ) = UNIT % IDENTITY
      FieldUnit ( GC % METRIC_F_UU_22 ) = CoordinateUnit ( 1 ) ** ( -2 )
      FieldUnit ( GC % METRIC_F_UU_33 ) = CoordinateUnit ( 1 ) ** ( -2 )
    end select !-- CoordinateSystem

    end associate !-- CoordinateUnit

  end subroutine SetUnits

  
  subroutine Compute_CGS ( GFSC, GC )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      GFSC  
    class ( Geometry_F_C_Form ), intent ( in ) :: &
      GC

  end subroutine Compute_CGS


  subroutine SetCoordinates ( GC, iD )

    class ( Geometry_F_C_Form ), intent ( inout ) :: &
      GC
    integer ( KDI ), intent ( in ) :: &
      iD      !-- iDimension

    integer ( KDI ) :: &
      iaF, iaL, &  !-- iaFirst, iaLast
      iC, &        !-- iCell
      oC           !-- oCell
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Edge_I_3D, &
      Width_3D, &
      Center_3D

    select type ( C  =>  GC % Chart )
    class is ( Chart_GS_Form )

    associate ( GV  =>  GC % Storage_FSC % Storage % Value )

    iaF  =  1  -  C % nGhostLayers ( iD ) 
    if ( C % Distributed ) then
      iaL  =  C % nCellsBrick ( iD )  +  C % nGhostLayers ( iD )
       oC  =  ( C % iaBrick ( iD )  -  1 )  *  C % nCellsBrick ( iD )
    else
      iaL  =  C % nCells ( iD )  +  C % nGhostLayers ( iD )
       oC  =  0
    end if

    call C % SetFieldPointer &
           ( GV ( :, GC % EDGE_I_U ( iD ) ), Edge_I_3D )
    call C % SetFieldPointer &
           ( GV ( :, GC % WIDTH_U ( iD ) ),  Width_3D )
    call C % SetFieldPointer &
           ( GV ( :, GC % CENTER_U ( iD ) ), Center_3D )

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

    end associate !-- GV

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Geometry_F_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCoordinates', 'subroutine', CONSOLE % ERROR )
    end select !-- C

  end subroutine SetCoordinates


  subroutine ComputeFromCoordinates ( GS, GC, nValuesOption, oValueOption )

    !-- Assumes coordinate fields are set

    class ( StorageForm ), intent ( inout ) :: &
      GS
    class ( Geometry_F_C_Form ), intent ( in ) :: &
      GC
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oValue, &
      nValues, &
      nDimensions

    oValue  =  0
    if ( present ( oValueOption ) ) &
      oValue  =  oValueOption

    nValues  =  GS % nValues
    if ( present ( nValuesOption ) ) &
      nValues  =  nValuesOption

    associate ( C  =>  GC % Chart )

    nDimensions  =  C % nDimensions

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call Compute_FV_R_Kernel &
             ( GS % Value ( :, GC % WIDTH_U_1 ), &
               GS % Value ( :, GC % WIDTH_U_2 ), &
               GS % Value ( :, GC % WIDTH_U_3 ), &
               nDimensions, &
               GS % Value ( :, GC % AREA_I_D_1 ), &
               GS % Value ( :, GC % AREA_I_D_2 ), &
               GS % Value ( :, GC % AREA_I_D_3 ), &
               GS % Value ( :, GC % VOLUME ) )
      call Compute_M_R_Kernel &
             ( GS % Value ( :, GC % METRIC_F_DD_11 ), &
               GS % Value ( :, GC % METRIC_F_DD_22 ), &
               GS % Value ( :, GC % METRIC_F_DD_33 ), &
               GS % Value ( :, GC % METRIC_F_UU_11 ), &
               GS % Value ( :, GC % METRIC_F_UU_22 ), &
               GS % Value ( :, GC % METRIC_F_UU_33 ) )
    case ( 'CYLINDRICAL' )
      call Compute_FV_C_Kernel &
             ( GS % Value ( :, GC % WIDTH_U_1 ), &
               GS % Value ( :, GC % WIDTH_U_2 ), &
               GS % Value ( :, GC % WIDTH_U_3 ), &
               GS % Value ( :, GC % EDGE_I_U_1 ), &
               nDimensions, &
               GS % Value ( :, GC % AREA_I_D_1 ), &
               GS % Value ( :, GC % AREA_I_D_2 ), &
               GS % Value ( :, GC % AREA_I_D_3 ), &
               GS % Value ( :, GC % VOLUME ) )
      call Compute_M_C_Kernel &
             ( GS % Value ( :, GC % CENTER_U_1 ), &
               nDimensions, &
               GS % Value ( :, GC % METRIC_F_DD_11 ), &
               GS % Value ( :, GC % METRIC_F_DD_22 ), &
               GS % Value ( :, GC % METRIC_F_DD_33 ), &
               GS % Value ( :, GC % METRIC_F_UU_11 ), &
               GS % Value ( :, GC % METRIC_F_UU_22 ), &
               GS % Value ( :, GC % METRIC_F_UU_33 ) )
    case ( 'SPHERICAL' )
      call Compute_FV_S_Kernel &
             ( GS % Value ( :, GC % WIDTH_U_1 ), &
               GS % Value ( :, GC % WIDTH_U_2 ), &
               GS % Value ( :, GC % WIDTH_U_3 ), &
               GS % Value ( :, GC % EDGE_I_U_1 ), &
               GS % Value ( :, GC % EDGE_I_U_2 ), &
               nDimensions, &
               GS % Value ( :, GC % AREA_I_D_1 ), &
               GS % Value ( :, GC % AREA_I_D_2 ), &
               GS % Value ( :, GC % AREA_I_D_3 ), &
               GS % Value ( :, GC % VOLUME ) )
      call Compute_M_S_Kernel &
             ( GS % Value ( :, GC % CENTER_U_1 ), &
               GS % Value ( :, GC % CENTER_U_2 ), &
               nDimensions, &
               GS % Value ( :, GC % METRIC_F_DD_11 ), &
               GS % Value ( :, GC % METRIC_F_DD_22 ), &
               GS % Value ( :, GC % METRIC_F_DD_33 ), &
               GS % Value ( :, GC % METRIC_F_UU_11 ), &
               GS % Value ( :, GC % METRIC_F_UU_22 ), &
               GS % Value ( :, GC % METRIC_F_UU_33 ) )
    case default
      call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', &
                  CONSOLE % ERROR )
      call Show ( 'Geometry_F_CH_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeFromCoordinates', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    end associate !-- C

  end subroutine ComputeFromCoordinates


end module Geometry_F_C__Form
