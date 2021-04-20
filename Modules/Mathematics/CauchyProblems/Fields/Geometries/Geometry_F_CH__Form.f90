module Geometry_F_CH__Form

  !-- Geometry_Flat_ChartHeader_Form

  use Basics
  use Manifolds
  use FieldSets

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_FLAT  = 19, &
      N_VECTORS_FLAT =  0

  type, public :: Geometry_F_CH_Form
!     integer ( KDI ) :: &
!       N_FIELDS_FLAT = N_FIELDS_FLAT, &
!       N_VECTORS_FLAT = N_VECTORS_FLAT
!     integer ( KDI ) :: &
!       !-- Coordinate fields
!       EDGE_I_U_1   = 0, &
!       EDGE_I_U_2   = 0, &
!       EDGE_I_U_3   = 0, &
!       WIDTH_U_1    = 0, &
!       WIDTH_U_2    = 0, &
!       WIDTH_U_3    = 0, &
!       CENTER_U_1   = 0, &
!       CENTER_U_2   = 0, &
!       CENTER_U_3   = 0
!     integer ( KDI ) :: &
!       !-- Finite volume fields
!       AREA_I_D_1   = 0, &
!       AREA_I_D_2   = 0, &
!       AREA_I_D_3   = 0, &
!       VOLUME       = 0
!     integer ( KDI ) :: &
!       !-- Flat metric fields
!       METRIC_F_DD_11 = 0, &
!       METRIC_F_DD_22 = 0, &
!       METRIC_F_DD_33 = 0, &
!       METRIC_F_UU_11 = 0, &
!       METRIC_F_UU_22 = 0, &
!       METRIC_F_UU_33 = 0
!     integer ( KDI ), dimension ( 3 ) :: &
!       EDGE_I_U, &
!       WIDTH_U, &
!       CENTER_U, &
!       AREA_I_D
    class ( FieldSet_CH_Form ), allocatable :: &
      FieldSet
  contains
    procedure, public, pass :: &
      Initialize_H
!     procedure, public, pass ( GC ) :: &
!       SetStream
!     procedure, public, pass :: &
!       ComputeFromCoordinates
    final :: &
      Finalize
  end type Geometry_F_CH_Form

!     private :: &
!       SetUnits, &
!       Compute_FV_R_Kernel, &
!       Compute_FV_C_Kernel, &
!       Compute_FV_S_Kernel, &
!       Compute_M_R_Kernel, &
!       Compute_M_C_Kernel, &
!       Compute_M_S_Kernel

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
               ( GC, C, NameOption ) !FieldOption, VectorOption, UnitOption, &
!                 VectorIndicesOption, nFieldsOption )

    class ( Geometry_F_CH_Form ), intent ( inout ) :: &
      GC
    class ( Chart_H_Form ), intent ( inout ) :: &
      C
!     class ( FieldSet_MH_Form ), intent ( in ), target :: &
!       GM
!     character ( * ), dimension ( : ), intent ( in ), optional :: &
!       FieldOption, &
!       VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
!     type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
!       UnitOption
!     type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
!       VectorIndicesOption
!     integer ( KDI ), intent ( in ), optional :: &
!       nFieldsOption

!     integer ( KDI ) :: &
!       nFields
!     type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
!       Unit
    character ( LDL ) :: &
      Name
!     character ( LDL ), dimension ( : ), allocatable :: &
!       Field

    Name  =  'Geometry'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

!     !-- Field indices

!     GC % EDGE_I_U_1      =   1
!     GC % EDGE_I_U_2      =   2
!     GC % EDGE_I_U_3      =   3
!     GC % WIDTH_U_1       =   4
!     GC % WIDTH_U_2       =   5
!     GC % WIDTH_U_3       =   6
!     GC % CENTER_U_1      =   7
!     GC % CENTER_U_2      =   8
!     GC % CENTER_U_3      =   9
!     GC % AREA_I_D_1      =  10
!     GC % AREA_I_D_2      =  11
!     GC % AREA_I_D_3      =  12
!     GC % VOLUME          =  13
!     GC % METRIC_F_DD_11  =  14
!     GC % METRIC_F_DD_22  =  15
!     GC % METRIC_F_DD_33  =  16
!     GC % METRIC_F_UU_11  =  17
!     GC % METRIC_F_UU_22  =  18
!     GC % METRIC_F_UU_33  =  19

!     nFields  =  GC % N_FIELDS_FLAT
!     if ( present ( nFieldsOption ) ) &
!       nFields  =  nFieldsOption

!     GC % EDGE_I_U  =  [ GC % EDGE_I_U_1, GC % EDGE_I_U_2, GC % EDGE_I_U_3 ]
!     GC % WIDTH_U   =  [ GC % WIDTH_U_1,  GC % WIDTH_U_2,  GC % WIDTH_U_3  ]
!     GC % CENTER_U  =  [ GC % CENTER_U_1, GC % CENTER_U_2, GC % CENTER_U_3 ]
!     GC % AREA_I_D  =  [ GC % AREA_I_D_1, GC % AREA_I_D_2, GC % AREA_I_D_3 ]

!     !-- Field names

!     if ( present ( FieldOption ) ) then
!       Field  =  FieldOption
!     else
!       allocate ( Field ( nFields ) )
!     end if

!     Field ( 1 : GC % N_FIELDS_FLAT ) &
!       = [ 'Edge_I_U_1    ', &
!           'Edge_I_U_2    ', &
!           'Edge_I_U_3    ', &
!           'Width_U_1     ', &
!           'Width_U_2     ', &
!           'Width_U_3     ', &
!           'Center_U_1    ', &
!           'Center_U_2    ', &
!           'Center_U_3    ', &
!           'Area_I_D_1    ', &
!           'Area_I_D_2    ', &
!           'Area_I_D_3    ', &
!           'Volume        ', &
!           'Metric_F_DD_11', &
!           'Metric_F_DD_22', &
!           'Metric_F_DD_33', &
!           'Metric_F_UU_11', &
!           'Metric_F_UU_22', &
!           'Metric_F_UU_33' ]

!     !-- Units

!     if ( present ( UnitOption ) ) then
!       Unit  =  UnitOption
!     else
!       allocate ( Unit ( nFields ) )
!     end if

!     call SetUnits ( Unit, GC, C )

    !-- FieldSet

    if ( .not. allocated ( GC % FieldSet ) ) then
      allocate ( GC % FieldSet )
      associate ( FS  =>  GC % FieldSet )
      FS % Type  =  'a Geometry_F_C'
      call FS % Initialize_H &
             ( C, &
               NameOption = Name )!nFields, FieldOption = Field, VectorOption = VectorOption, &
!              UnitOption = Unit, VectorIndicesOption = VectorIndicesOption )
      end associate !-- FS
    end if

  end subroutine Initialize_H


!   subroutine SetStream ( G_Stream, GC, G_Source )

!     class ( StorageForm ), intent ( inout ) :: &
!       G_Stream
!     class ( Geometry_F_CH_Form ), intent ( in ) :: &
!       GC
!     class ( StorageForm ), intent ( in ) :: &
!       G_Source

!     call G_Stream % Initialize &
!            ( G_Source, &
!              iaSelectedOption &
!                =  [ GC % CENTER_U_1, GC % CENTER_U_2, GC % CENTER_U_3, &
!                     GC % METRIC_F_DD_11, GC % METRIC_F_DD_22, &
!                     GC % METRIC_F_DD_33 ] )

!   end subroutine SetStream


!   subroutine ComputeFromCoordinates ( GC, G, nValuesOption, oValueOption )

!     !-- Assumes coordinate fields are set

!     class ( Geometry_F_CH_Form ), intent ( inout ) :: &
!       GC
!     class ( StorageForm ), intent ( inout ) :: &
!       G
!     integer ( KDI ), intent ( in ), optional :: &
!       nValuesOption, &
!       oValueOption

!     integer ( KDI ) :: &
!       oValue, &
!       nValues, &
!       nDimensions

!     oValue  =  0
!     if ( present ( oValueOption ) ) &
!       oValue  =  oValueOption

!     nValues  =  G % nValues
!     if ( present ( nValuesOption ) ) &
!       nValues  =  nValuesOption

!     nDimensions  =  GC % Chart % nDimensions

!     select case ( trim ( GC % Chart % CoordinateSystem ) )
!     case ( 'RECTANGULAR' )
!       call Compute_FV_R_Kernel &
!              ( G % Value ( :, GC % AREA_I_D_1 ), &
!                G % Value ( :, GC % AREA_I_D_2 ), &
!                G % Value ( :, GC % AREA_I_D_3 ), &
!                G % Value ( :, GC % VOLUME ), &
!                G % Value ( :, GC % WIDTH_U_1 ), &
!                G % Value ( :, GC % WIDTH_U_2 ), &
!                G % Value ( :, GC % WIDTH_U_3 ), &
!                nDimensions, nValues, oValue )
!       call Compute_M_R_Kernel &
!              ( G % Value ( :, GC % METRIC_F_DD_11 ), &
!                G % Value ( :, GC % METRIC_F_DD_22 ), &
!                G % Value ( :, GC % METRIC_F_DD_33 ), &
!                G % Value ( :, GC % METRIC_F_UU_11 ), &
!                G % Value ( :, GC % METRIC_F_UU_22 ), &
!                G % Value ( :, GC % METRIC_F_UU_33 ), &
!                nValues, oValue )
!     case ( 'CYLINDRICAL' )
!       call Compute_FV_C_Kernel &
!              ( G % Value ( :, GC % AREA_I_D_1 ), &
!                G % Value ( :, GC % AREA_I_D_2 ), &
!                G % Value ( :, GC % AREA_I_D_3 ), &
!                G % Value ( :, GC % VOLUME ), &
!                G % Value ( :, GC % WIDTH_U_1 ), &
!                G % Value ( :, GC % WIDTH_U_2 ), &
!                G % Value ( :, GC % WIDTH_U_3 ), &
!                G % Value ( :, GC % EDGE_I_U_1 ), &
!                nDimensions, nValues, oValue )
!       call Compute_M_C_Kernel &
!              ( G % Value ( :, GC % METRIC_F_DD_11 ), &
!                G % Value ( :, GC % METRIC_F_DD_22 ), &
!                G % Value ( :, GC % METRIC_F_DD_33 ), &
!                G % Value ( :, GC % METRIC_F_UU_11 ), &
!                G % Value ( :, GC % METRIC_F_UU_22 ), &
!                G % Value ( :, GC % METRIC_F_UU_33 ), &
!                G % Value ( :, GC % CENTER_U_1 ), &
!                nDimensions, nValues, oValue )
!     case ( 'SPHERICAL' )
!       call Compute_FV_S_Kernel &
!              ( G % Value ( :, GC % AREA_I_D_1 ), &
!                G % Value ( :, GC % AREA_I_D_2 ), &
!                G % Value ( :, GC % AREA_I_D_3 ), &
!                G % Value ( :, GC % VOLUME ), &
!                G % Value ( :, GC % WIDTH_U_1 ), &
!                G % Value ( :, GC % WIDTH_U_2 ), &
!                G % Value ( :, GC % WIDTH_U_3 ), &
!                G % Value ( :, GC % EDGE_I_U_1 ), &
!                G % Value ( :, GC % EDGE_I_U_2 ), &
!                nDimensions, nValues, oValue )
!       call Compute_M_S_Kernel &
!              ( G % Value ( :, GC % METRIC_F_DD_11 ), &
!                G % Value ( :, GC % METRIC_F_DD_22 ), &
!                G % Value ( :, GC % METRIC_F_DD_33 ), &
!                G % Value ( :, GC % METRIC_F_UU_11 ), &
!                G % Value ( :, GC % METRIC_F_UU_22 ), &
!                G % Value ( :, GC % METRIC_F_UU_33 ), &
!                G % Value ( :, GC % CENTER_U_1 ), &
!                G % Value ( :, GC % CENTER_U_2 ), &
!                nDimensions, nValues, oValue )
!     case default
!       call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
!       call Show ( GC % Chart % CoordinateSystem, 'CoordinateSystem', &
!                   CONSOLE % ERROR )
!       call Show ( 'Geometry_F_Form', 'module', CONSOLE % ERROR )
!       call Show ( 'ComputeFromCoordinates', 'subroutine', CONSOLE % ERROR )
!       call PROGRAM_HEADER % Abort ( )
!     end select

!   end subroutine ComputeFromCoordinates


  impure elemental subroutine Finalize ( GC )

    type ( Geometry_F_CH_Form ), intent ( inout ) :: &
      GC

  end subroutine Finalize

  
!   subroutine SetUnits ( FieldUnit, GC, C )

!     type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
!       FieldUnit
!     class ( Geometry_F_CH_Form ), intent ( in ) :: &
!       GC
!     class ( Chart_H_Form ), intent ( in ) :: &
!       C

!     associate &
!       ( CoordinateUnit    =>  C % CoordinateUnit, &
!         CoordinateSystem  =>  C % CoordinateSystem )

!     FieldUnit ( GC % EDGE_I_U_1 : GC % EDGE_I_U_3 ) &
!       = CoordinateUnit
!     FieldUnit ( GC % WIDTH_U_1 : GC % WIDTH_U_3 ) &
!       = CoordinateUnit
!     FieldUnit ( GC % CENTER_U_1 : GC % CENTER_U_3 ) &
!       = CoordinateUnit

!     select case ( trim ( CoordinateSystem ) )
!     case ( 'RECTANGULAR' )
!       FieldUnit ( GC % VOLUME )  &
!         =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )  &
!            *  CoordinateUnit ( 3 )
!       FieldUnit ( GC % AREA_I_D_1 )  &
!         =  CoordinateUnit ( 2 )  *  CoordinateUnit ( 3 )
!       FieldUnit ( GC % AREA_I_D_2 )  &
!         =  CoordinateUnit ( 3 )  *  CoordinateUnit ( 1 )
!       FieldUnit ( GC % AREA_I_D_3 )  &
!         =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
!       FieldUnit ( GC % METRIC_F_DD_11 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_DD_22 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_DD_33 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_UU_11 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_UU_22 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_UU_33 ) = UNIT % IDENTITY
!     case ( 'CYLINDRICAL' )
!       FieldUnit ( GC % VOLUME )  &
!         =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
!       FieldUnit ( GC % AREA_I_D_1 )  &
!         =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
!       FieldUnit ( GC % AREA_I_D_2 )  &
!         =  CoordinateUnit ( 1 ) ** 2
!       FieldUnit ( GC % AREA_I_D_3 )  &
!         =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
!       FieldUnit ( GC % METRIC_F_DD_11 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_DD_22 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_DD_33 ) = CoordinateUnit ( 1 ) ** (  2 )
!       FieldUnit ( GC % METRIC_F_UU_11 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_UU_22 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_UU_33 ) = CoordinateUnit ( 1 ) ** ( -2 )
!     case ( 'SPHERICAL' )
!       FieldUnit ( GC % VOLUME )  &
!         = CoordinateUnit ( 1 ) ** 3
!       FieldUnit ( GC % AREA_I_D_1 )  &
!         =  CoordinateUnit ( 1 ) ** 2
!       FieldUnit ( GC % AREA_I_D_2 )  &
!         =  CoordinateUnit ( 1 ) ** 3
!       FieldUnit ( GC % AREA_I_D_3 )  &
!         =  CoordinateUnit ( 1 ) ** 3
!       FieldUnit ( GC % METRIC_F_DD_11 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_DD_22 ) = CoordinateUnit ( 1 ) ** (  2 )
!       FieldUnit ( GC % METRIC_F_DD_33 ) = CoordinateUnit ( 1 ) ** (  2 )
!       FieldUnit ( GC % METRIC_F_UU_11 ) = UNIT % IDENTITY
!       FieldUnit ( GC % METRIC_F_UU_22 ) = CoordinateUnit ( 1 ) ** ( -2 )
!       FieldUnit ( GC % METRIC_F_UU_33 ) = CoordinateUnit ( 1 ) ** ( -2 )
!     end select !-- CoordinateSystem

!     end associate !-- CoordinateUnit

!   end subroutine SetUnits

  
end module Geometry_F_CH__Form
