module Geometry_F_CH__Form

  !-- Geometry_Flat_ChartHeader_Form

  use Basics
  use ManifoldBasics
  use Chart_H__Form
  use FieldSet_CH__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_FLAT  = 19, &
      N_VECTORS_FLAT =  0

  type, public, extends ( FieldSet_CH_Form ) :: Geometry_F_CH_Form
    integer ( KDI ) :: &
      N_FIELDS_FLAT = N_FIELDS_FLAT, &
      N_VECTORS_FLAT = N_VECTORS_FLAT
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
      Initialize_F
    generic, public :: &
      Initialize => Initialize_F
    procedure, public, pass ( GC ) :: &
      SetStream
    final :: &
      Finalize
  end type Geometry_F_CH_Form

    private :: &
      SetUnits

contains


  subroutine Initialize_F &
               ( GC, C, GM, FieldOption, VectorOption, UnitOption, &
                 VectorIndicesOption, nFieldsOption )

    class ( Geometry_F_CH_Form ), intent ( inout ) :: &
      GC
    class ( Chart_H_Form ), intent ( inout ), target :: &
      C
    class ( FieldSet_MH_Form ), intent ( in ), target :: &
      GM
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      nFields
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( GC % Type == '' ) &
      GC % Type = 'a Geometry_F_C'

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

    nFields  =  GC % N_FIELDS_FLAT
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    GC % EDGE_I_U  =  [ GC % EDGE_I_U_1, GC % EDGE_I_U_2, GC % EDGE_I_U_3 ]
    GC % WIDTH_U   =  [ GC % WIDTH_U_1,  GC % WIDTH_U_2,  GC % WIDTH_U_3  ]
    GC % CENTER_U  =  [ GC % CENTER_U_1, GC % CENTER_U_2, GC % CENTER_U_3 ]
    GC % AREA_I_D  =  [ GC % AREA_I_D_1, GC % AREA_I_D_2, GC % AREA_I_D_3 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      Field  =  FieldOption
    else
      allocate ( Field ( nFields ) )
    end if

    Field ( 1 : GC % N_FIELDS_FLAT ) &
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
      Unit  =  UnitOption
    else
      allocate ( Unit ( nFields ) )
    end if

    call SetUnits ( Unit, GC, C )

    !-- Parent initialization

    call GC % FieldSet_CH_Form % Initialize &
           ( C, GM, nFields, FieldOption = Field, VectorOption = VectorOption, &
             UnitOption = Unit, VectorIndicesOption = VectorIndicesOption )

  end subroutine Initialize_F


  subroutine SetStream ( G_Stream, GC, G_Source )

    class ( StorageForm ), intent ( inout ) :: &
      G_Stream
    class ( Geometry_F_CH_Form ), intent ( in ) :: &
      GC
    class ( StorageForm ), intent ( in ) :: &
      G_Source

    call G_Stream % Initialize &
           ( G_Source, &
             iaSelectedOption &
               =  [ GC % CENTER_U_1, GC % CENTER_U_2, GC % CENTER_U_3, &
                    GC % METRIC_F_DD_11, GC % METRIC_F_DD_22, &
                    GC % METRIC_F_DD_33 ] )

  end subroutine SetStream


  impure elemental subroutine Finalize ( GC )

    type ( Geometry_F_CH_Form ), intent ( inout ) :: &
      GC

  end subroutine Finalize

  
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
