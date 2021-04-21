module Geometry_F_CH__Form

  !-- Geometry_Flat_ChartHeader_Form

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_F  = 19, &
      N_VECTORS_F =  0

  type, public :: Geometry_F_CH_Form
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
    class ( FieldSet_CH_Form ), allocatable :: &
      FieldSet
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, public, pass ( GC ) :: &
      SetStream
    procedure, public, pass ( GC ):: &
      ComputeFromCoordinates
    procedure, public, pass :: &
      Show => Show_GC
    final :: &
      Finalize
  end type Geometry_F_CH_Form

  type, public :: Geometry_C_Element
    !-- Geometry_Chart_Element
    class ( Geometry_F_CH_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type Geometry_C_Element

    private :: &
      SetUnits, &
      Compute_FV_R_Kernel, &
      Compute_FV_C_Kernel, &
      Compute_FV_S_Kernel, &
      Compute_M_R_Kernel, &
      Compute_M_C_Kernel, &
      Compute_M_S_Kernel

    interface
      
      module subroutine Compute_FV_R_Kernel &
               ( A_I_1, A_I_2, A_I_3, V, W_1, W_2, W_3, nD, nV, oV )
        !-- Compute_FiniteVolume_Rectangular_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          A_I_1, A_I_2, A_I_3, &
          V
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3
        integer ( KDI ), intent ( in ) :: &
          nD, &
          nV, &
          oV
      end subroutine Compute_FV_R_Kernel

      module subroutine Compute_FV_C_Kernel &
               ( A_I_1, A_I_2, A_I_3, V, W_1, W_2, W_3, E_I_1, nD, nV, oV )
        !-- Compute_FiniteVolume_Cylindrical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          A_I_1, A_I_2, A_I_3, &
          V
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3, &
          E_I_1
        integer ( KDI ), intent ( in ) :: &
          nD, &
          nV, &
          oV
      end subroutine Compute_FV_C_Kernel

      module subroutine Compute_FV_S_Kernel &
               ( A_I_1, A_I_2, A_I_3, V, W_1, W_2, W_3, E_I_1, E_I_2, &
                 nD, nV, oV )
        !-- Compute_FiniteVolume_Spherical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          A_I_1, A_I_2, A_I_3, &
          V
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          W_1, W_2, W_3, &
          E_I_1, E_I_2
        integer ( KDI ), intent ( in ) :: &
          nD, &
          nV, &
          oV
      end subroutine Compute_FV_S_Kernel

      module subroutine Compute_M_R_Kernel &
               ( M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 nV, oV, UseDeviceOption )
        !-- Compute_Metric_Rectangular_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_11, M_UU_22, M_UU_33
        integer ( KDI ), intent ( in ) :: &
          nV, &
          oV
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_M_R_Kernel

      module subroutine Compute_M_C_Kernel &
               ( M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 RP, nD, nV, oV, UseDeviceOption )
        !-- Compute_Metric_Cylindrical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_11, M_UU_22, M_UU_33
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          RP
        integer ( KDI ), intent ( in ) :: &
          nD, &
          nV, &
          oV
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_M_C_Kernel

      module subroutine Compute_M_S_Kernel &
               ( M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 R, Th, nD, nV, oV, UseDeviceOption )
        !-- Compute_Metric_Spherical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_11, M_UU_22, M_UU_33
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          R, Th
        integer ( KDI ), intent ( in ) :: &
          nD, &
          nV, &
          oV
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_M_S_Kernel

    end interface


contains


  subroutine Initialize_H &
               ( GC, C, NameOption, nFieldsOption, FieldOption, UnitOption )

    class ( Geometry_F_CH_Form ), intent ( inout ) :: &
      GC
    class ( Chart_H_Form ), intent ( inout ) :: &
      C
    character ( * ), intent ( inout ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption
    character ( * ), dimension ( : ), intent ( out ), allocatable, optional :: &
      FieldOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( out ), allocatable, &
      optional :: &
        UnitOption

    integer ( KDI ) :: &
      nFields
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    Name  =  'Geometry'
    if ( present ( NameOption ) ) then
      if ( NameOption  ==  '' ) then
        NameOption  =  Name
      else
        Name  =  NameOption
      end if
    end if

    !-- Field indices

    GC % EDGE_I_U_1      =   1
    GC % EDGE_I_U_2      =   2
    GC % EDGE_I_U_3      =   3
    GC % WIDTH_U_1       =   4
    GC % WIDTH_U_2       =   5
    GC % WIDTH_U_3       =   6
    GC % CENTER_U_1      =   7
    GC % CENTER_U_2      =   8
    GC % CENTER_U_3      =   9
    GC % AREA_I_D_1      =  10
    GC % AREA_I_D_2      =  11
    GC % AREA_I_D_3      =  12
    GC % VOLUME          =  13
    GC % METRIC_F_DD_11  =  14
    GC % METRIC_F_DD_22  =  15
    GC % METRIC_F_DD_33  =  16
    GC % METRIC_F_UU_11  =  17
    GC % METRIC_F_UU_22  =  18
    GC % METRIC_F_UU_33  =  19

    nFields  =  GC % N_FIELDS_F
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    GC % EDGE_I_U  =  [ GC % EDGE_I_U_1, GC % EDGE_I_U_2, GC % EDGE_I_U_3 ]
    GC % WIDTH_U   =  [ GC % WIDTH_U_1,  GC % WIDTH_U_2,  GC % WIDTH_U_3  ]
    GC % CENTER_U  =  [ GC % CENTER_U_1, GC % CENTER_U_2, GC % CENTER_U_3 ]
    GC % AREA_I_D  =  [ GC % AREA_I_D_1, GC % AREA_I_D_2, GC % AREA_I_D_3 ]

    !-- Field names

    allocate ( Field ( nFields ) )

    Field ( 1 : GC % N_FIELDS_F ) &
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

    if ( present ( FieldOption ) ) &
      allocate ( FieldOption, source = Field )

    !-- Units

    allocate ( Unit ( nFields ) )

    call SetUnits ( Unit, GC, C )

    if ( present ( UnitOption ) ) &
      allocate ( UnitOption, source = Unit )

    !-- FieldSet

    if ( .not. allocated ( GC % FieldSet ) ) then

      allocate ( GC % FieldSet )
      associate ( FS  =>  GC % FieldSet )

      FS % Type  =  'a Geometry_F_C'    

      call FS % Initialize_H &
             ( C, &
               FieldOption = Field, &
               NameOption = Name, &
               UnitOption = Unit, &
               nFieldsOption = nFields )

      end associate !-- FS
    end if

  end subroutine Initialize_H


  subroutine SetStream ( SC, GC )

    class ( Stream_CH_Form ), intent ( inout ) :: &
      SC
    class ( Geometry_F_CH_Form ), intent ( in ) :: &
      GC

    call SC % AddFieldSet &
           ( GC % FieldSet, &
             iaSelectedOption &
               =  [ GC % CENTER_U_1, GC % CENTER_U_2, GC % CENTER_U_3, &
                    GC % METRIC_F_DD_11, GC % METRIC_F_DD_22, &
                    GC % METRIC_F_DD_33 ] )

  end subroutine SetStream


  subroutine ComputeFromCoordinates ( GS, GC, nValuesOption, oValueOption )

    !-- Assumes coordinate fields are set

    class ( StorageForm ), intent ( inout ) :: &
      GS
    class ( Geometry_F_CH_Form ), intent ( in ) :: &
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

    associate ( C  =>  GC % FieldSet % Chart )

    nDimensions  =  C % nDimensions

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call Compute_FV_R_Kernel &
             ( GS % Value ( :, GC % AREA_I_D_1 ), &
               GS % Value ( :, GC % AREA_I_D_2 ), &
               GS % Value ( :, GC % AREA_I_D_3 ), &
               GS % Value ( :, GC % VOLUME ), &
               GS % Value ( :, GC % WIDTH_U_1 ), &
               GS % Value ( :, GC % WIDTH_U_2 ), &
               GS % Value ( :, GC % WIDTH_U_3 ), &
               nDimensions, nValues, oValue )
      call Compute_M_R_Kernel &
             ( GS % Value ( :, GC % METRIC_F_DD_11 ), &
               GS % Value ( :, GC % METRIC_F_DD_22 ), &
               GS % Value ( :, GC % METRIC_F_DD_33 ), &
               GS % Value ( :, GC % METRIC_F_UU_11 ), &
               GS % Value ( :, GC % METRIC_F_UU_22 ), &
               GS % Value ( :, GC % METRIC_F_UU_33 ), &
               nValues, oValue )
    case ( 'CYLINDRICAL' )
      call Compute_FV_C_Kernel &
             ( GS % Value ( :, GC % AREA_I_D_1 ), &
               GS % Value ( :, GC % AREA_I_D_2 ), &
               GS % Value ( :, GC % AREA_I_D_3 ), &
               GS % Value ( :, GC % VOLUME ), &
               GS % Value ( :, GC % WIDTH_U_1 ), &
               GS % Value ( :, GC % WIDTH_U_2 ), &
               GS % Value ( :, GC % WIDTH_U_3 ), &
               GS % Value ( :, GC % EDGE_I_U_1 ), &
               nDimensions, nValues, oValue )
      call Compute_M_C_Kernel &
             ( GS % Value ( :, GC % METRIC_F_DD_11 ), &
               GS % Value ( :, GC % METRIC_F_DD_22 ), &
               GS % Value ( :, GC % METRIC_F_DD_33 ), &
               GS % Value ( :, GC % METRIC_F_UU_11 ), &
               GS % Value ( :, GC % METRIC_F_UU_22 ), &
               GS % Value ( :, GC % METRIC_F_UU_33 ), &
               GS % Value ( :, GC % CENTER_U_1 ), &
               nDimensions, nValues, oValue )
    case ( 'SPHERICAL' )
      call Compute_FV_S_Kernel &
             ( GS % Value ( :, GC % AREA_I_D_1 ), &
               GS % Value ( :, GC % AREA_I_D_2 ), &
               GS % Value ( :, GC % AREA_I_D_3 ), &
               GS % Value ( :, GC % VOLUME ), &
               GS % Value ( :, GC % WIDTH_U_1 ), &
               GS % Value ( :, GC % WIDTH_U_2 ), &
               GS % Value ( :, GC % WIDTH_U_3 ), &
               GS % Value ( :, GC % EDGE_I_U_1 ), &
               GS % Value ( :, GC % EDGE_I_U_2 ), &
               nDimensions, nValues, oValue )
      call Compute_M_S_Kernel &
             ( GS % Value ( :, GC % METRIC_F_DD_11 ), &
               GS % Value ( :, GC % METRIC_F_DD_22 ), &
               GS % Value ( :, GC % METRIC_F_DD_33 ), &
               GS % Value ( :, GC % METRIC_F_UU_11 ), &
               GS % Value ( :, GC % METRIC_F_UU_22 ), &
               GS % Value ( :, GC % METRIC_F_UU_33 ), &
               GS % Value ( :, GC % CENTER_U_1 ), &
               GS % Value ( :, GC % CENTER_U_2 ), &
               nDimensions, nValues, oValue )
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


  subroutine Show_GC ( GC )

    class ( Geometry_F_CH_Form ), intent ( inout ) :: &
      GC

    call GC % FieldSet % Show ( )

  end subroutine Show_GC


  impure elemental subroutine Finalize ( GC )

    type ( Geometry_F_CH_Form ), intent ( inout ) :: &
      GC

    if ( allocated ( GC % FieldSet ) ) &
      deallocate ( GC % FieldSet )

  end subroutine Finalize

  
  impure elemental subroutine Finalize_E ( GE )
    
    type ( Geometry_C_Element ), intent ( inout ) :: &
      GE

    if ( allocated ( GE % Element ) ) &
      deallocate ( GE % Element )

  end subroutine Finalize_E


  subroutine SetUnits ( FieldUnit, GC, C )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      FieldUnit
    class ( Geometry_F_CH_Form ), intent ( in ) :: &
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

  
end module Geometry_F_CH__Form
