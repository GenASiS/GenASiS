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
    final :: &
      Finalize
  end type Geometry_F_CH_Form


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

    nFields  =  GC % N_FIELDS_FLAT
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

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

    if ( present ( UnitOption ) ) then
      Unit  =  UnitOption
    else
      allocate ( Unit ( nFields ) )
    end if

    call GC % FieldSet_CH_Form % Initialize &
           ( C, GM, nFields, FieldOption = Field, VectorOption = VectorOption, &
             UnitOption = Unit, VectorIndicesOption = VectorIndicesOption )

  end subroutine Initialize_F


  impure elemental subroutine Finalize ( GC )

    type ( Geometry_F_CH_Form ), intent ( inout ) :: &
      GC

  end subroutine Finalize

  
end module Geometry_F_CH__Form
