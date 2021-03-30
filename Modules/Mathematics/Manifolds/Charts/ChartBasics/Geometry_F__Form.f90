module Geometry_F__Form

  !-- Geometry_Flat_Form
  
  use Basics
  use FieldSet_CH__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_FLAT  = 19, &
      N_VECTORS_FLAT = 0

  type, public, extends ( StorageForm ) :: Geometry_F_Form
    integer ( KDI ) :: &
      IGNORABILITY   = 0, &
      N_FIELDS       = 0, &
      N_VECTORS      = 0, &
      N_FIELDS_FLAT  = N_FIELDS_FLAT, &
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
    character ( LDF ) :: &
      Type = '', &
      CoordinateSystem = ''
    class ( FieldSet_CH_Form ), pointer :: &
      Geometry_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_G_F
    generic, public :: &
      Initialize => InitializeAllocate_G_F
    procedure, public, pass :: &
      SetStream
    procedure, public, pass :: &
      ComputeFromCoordinates
    final :: &  !-- FIXME: Intel doesn't like final procedure name to be the
                !          same as the parent's final
      Finalize_G_F
  end type Geometry_F_Form

    private :: &
      InitializeBasics, &
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


  subroutine InitializeAllocate_G_F &
               ( G, GC, nValues, VariableOption, VectorOption, UnitOption, &
                 VectorIndicesOption )

    class ( Geometry_F_Form ), intent ( inout ) :: &
      G
    class ( FieldSet_CH_Form ), intent ( in ), target :: &
      GC
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector

    G % Geometry_C  =>  GC
    G % CoordinateSystem  =  GC % Chart % CoordinateSystem

    call InitializeBasics &
           ( G, GC % Name, Variable, Vector, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, UnitOption, VectorIndicesOption )

    call SetUnits ( VariableUnit, G, GC % Chart % CoordinateUnit )

    call G % StorageForm % Initialize &
          ( [ nValues, G % N_FIELDS ], &
            VariableOption = Variable, VectorOption = Vector, &
            NameOption = GC % Name, ClearOption = .true., &
            PinnedOption = GC % PinnedMemory, &
            UnitOption = VariableUnit, &
            VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_G_F


  subroutine SetStream ( G, G_Stream )

    class ( Geometry_F_Form ), intent ( inout ) :: &
      G
    class ( StorageForm ), intent ( inout ) :: &
      G_Stream

    call G_Stream % Initialize &
           ( G, iaSelectedOption &
                  = [ G % CENTER_U_1, G % CENTER_U_2, G % CENTER_U_3, &
                      G % METRIC_F_DD_11, G % METRIC_F_DD_22, &
                      G % METRIC_F_DD_33 ] )

  end subroutine SetStream


  subroutine ComputeFromCoordinates ( G, nValuesOption, oValueOption )

    !-- Assumes coordinate fields are set

    class ( Geometry_F_Form ), intent ( inout ) :: &
      G
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

    nValues  =  G % nValues
    if ( present ( nValuesOption ) ) &
      nValues  =  nValuesOption

    nDimensions  =  G % Geometry_C % Chart % nDimensions

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call Compute_FV_R_Kernel &
             ( G % Value ( :, G % AREA_I_D_1 ), &
               G % Value ( :, G % AREA_I_D_2 ), &
               G % Value ( :, G % AREA_I_D_3 ), &
               G % Value ( :, G % VOLUME ), &
               G % Value ( :, G % WIDTH_U_1 ), &
               G % Value ( :, G % WIDTH_U_2 ), &
               G % Value ( :, G % WIDTH_U_3 ), &
               nDimensions, nValues, oValue )
      call Compute_M_R_Kernel &
             ( G % Value ( :, G % METRIC_F_DD_11 ), &
               G % Value ( :, G % METRIC_F_DD_22 ), &
               G % Value ( :, G % METRIC_F_DD_33 ), &
               G % Value ( :, G % METRIC_F_UU_11 ), &
               G % Value ( :, G % METRIC_F_UU_22 ), &
               G % Value ( :, G % METRIC_F_UU_33 ), &
               nValues, oValue )
    case ( 'CYLINDRICAL' )
      call Compute_FV_C_Kernel &
             ( G % Value ( :, G % AREA_I_D_1 ), &
               G % Value ( :, G % AREA_I_D_2 ), &
               G % Value ( :, G % AREA_I_D_3 ), &
               G % Value ( :, G % VOLUME ), &
               G % Value ( :, G % WIDTH_U_1 ), &
               G % Value ( :, G % WIDTH_U_2 ), &
               G % Value ( :, G % WIDTH_U_3 ), &
               G % Value ( :, G % EDGE_I_U_1 ), &
               nDimensions, nValues, oValue )
      call Compute_M_C_Kernel &
             ( G % Value ( :, G % METRIC_F_DD_11 ), &
               G % Value ( :, G % METRIC_F_DD_22 ), &
               G % Value ( :, G % METRIC_F_DD_33 ), &
               G % Value ( :, G % METRIC_F_UU_11 ), &
               G % Value ( :, G % METRIC_F_UU_22 ), &
               G % Value ( :, G % METRIC_F_UU_33 ), &
               G % Value ( :, G % CENTER_U_1 ), &
               nDimensions, nValues, oValue )
    case ( 'SPHERICAL' )
      call Compute_FV_S_Kernel &
             ( G % Value ( :, G % AREA_I_D_1 ), &
               G % Value ( :, G % AREA_I_D_2 ), &
               G % Value ( :, G % AREA_I_D_3 ), &
               G % Value ( :, G % VOLUME ), &
               G % Value ( :, G % WIDTH_U_1 ), &
               G % Value ( :, G % WIDTH_U_2 ), &
               G % Value ( :, G % WIDTH_U_3 ), &
               G % Value ( :, G % EDGE_I_U_1 ), &
               G % Value ( :, G % EDGE_I_U_2 ), &
               nDimensions, nValues, oValue )
      call Compute_M_S_Kernel &
             ( G % Value ( :, G % METRIC_F_DD_11 ), &
               G % Value ( :, G % METRIC_F_DD_22 ), &
               G % Value ( :, G % METRIC_F_DD_33 ), &
               G % Value ( :, G % METRIC_F_UU_11 ), &
               G % Value ( :, G % METRIC_F_UU_22 ), &
               G % Value ( :, G % METRIC_F_UU_33 ), &
               G % Value ( :, G % CENTER_U_1 ), &
               G % Value ( :, G % CENTER_U_2 ), &
               nDimensions, nValues, oValue )
    case default
      call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
      call Show ( G % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'Geometry_F_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeFromCoordinates', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine ComputeFromCoordinates


  impure elemental subroutine Finalize_G_F ( G )

    type ( Geometry_F_Form ), intent ( inout ) :: &
      G

    if ( G % Name == '' ) return

    nullify ( G % Geometry_C )

    call Show ( 'Finalizing a ' // trim ( G % Type ), G % IGNORABILITY )
    call Show ( G % Name, 'Name', G % IGNORABILITY )

  end subroutine Finalize_G_F


  subroutine InitializeBasics &
               ( G, Name, Variable, Vector, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, VariableUnitOption, &
                 VectorIndicesOption )

    class ( Geometry_F_Form ), intent ( inout ) :: &
      G
    character ( LDF ), intent ( in ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable, &
      Vector
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    !-- FIXME: intent(out) here caused ICE with Intel Compiler 15
    !          Temporarily set to intent(inout)
    !type ( Integer_1D_Form ), dimension ( : ), allocatable, &
    !  intent ( out ) :: &
    type ( Integer_1D_Form ), dimension ( : ), allocatable, &
      intent ( inout ) :: &
        VectorIndices
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      iV  !-- iVector

    if ( G % Type == '' ) &
      G % Type = 'Geometry_F'

    G % IGNORABILITY = CONSOLE % INFO_4
    call Show ( 'Initializing a ' // trim ( G % Type ), G % IGNORABILITY )
    call Show ( Name, 'Name', G % IGNORABILITY )
    call Show ( G % CoordinateSystem, 'CoordinateSystem', G % IGNORABILITY )

    !-- variable indices

    G % EDGE_I_U_1      =   1
    G % EDGE_I_U_2      =   2
    G % EDGE_I_U_3      =   3
    G % WIDTH_U_1       =   4
    G % WIDTH_U_2       =   5
    G % WIDTH_U_3       =   6
    G % CENTER_U_1      =   7
    G % CENTER_U_2      =   8
    G % CENTER_U_3      =   9
    G % AREA_I_D_1      =  10
    G % AREA_I_D_2      =  11
    G % AREA_I_D_3      =  12
    G % VOLUME          =  13
    G % METRIC_F_DD_11  =  14
    G % METRIC_F_DD_22  =  15
    G % METRIC_F_DD_33  =  16
    G % METRIC_F_UU_11  =  17
    G % METRIC_F_UU_22  =  18
    G % METRIC_F_UU_33  =  19

    if ( G % N_FIELDS == 0 ) &
      G % N_FIELDS = G % N_FIELDS_FLAT

    G % EDGE_I_U  =  [ G % EDGE_I_U_1, G % EDGE_I_U_2, G % EDGE_I_U_3 ]
    G % WIDTH_U   =  [ G % WIDTH_U_1,  G % WIDTH_U_2,  G % WIDTH_U_3  ]
    G % CENTER_U  =  [ G % CENTER_U_1, G % CENTER_U_2, G % CENTER_U_3 ]
    G % AREA_I_D  =  [ G % AREA_I_D_1, G % AREA_I_D_2, G % AREA_I_D_3 ]

    !-- variable names

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( G % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : G % N_FIELDS_FLAT ) &
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

    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( G % N_FIELDS ) )
    end if
    
    !-- vectors

    if ( G % N_VECTORS == 0 ) &
      G % N_VECTORS = G % N_VECTORS_FLAT

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( G % N_VECTORS ) )
      Vector = ''
    end if

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = G % N_VECTORS_FLAT + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( G % N_VECTORS ) )
    end if

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, G, CoordinateUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      CoordinateUnit

    VariableUnit ( G % EDGE_I_U_1 : G % EDGE_I_U_3 ) &
      = CoordinateUnit
    VariableUnit ( G % WIDTH_U_1 : G % WIDTH_U_3 ) &
      = CoordinateUnit
    VariableUnit ( G % CENTER_U_1 : G % CENTER_U_3 ) &
      = CoordinateUnit

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      VariableUnit ( G % VOLUME )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )  &
           *  CoordinateUnit ( 3 )
      VariableUnit ( G % AREA_I_D_1 )  &
        =  CoordinateUnit ( 2 )  *  CoordinateUnit ( 3 )
      VariableUnit ( G % AREA_I_D_2 )  &
        =  CoordinateUnit ( 3 )  *  CoordinateUnit ( 1 )
      VariableUnit ( G % AREA_I_D_3 )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
      VariableUnit ( G % METRIC_F_DD_11 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_DD_22 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_DD_33 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_UU_11 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_UU_22 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_UU_33 ) = UNIT % IDENTITY
    case ( 'CYLINDRICAL' )
      VariableUnit ( G % VOLUME )  &
        =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
      VariableUnit ( G % AREA_I_D_1 )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
      VariableUnit ( G % AREA_I_D_2 )  &
        =  CoordinateUnit ( 1 ) ** 2
      VariableUnit ( G % AREA_I_D_3 )  &
        =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
      VariableUnit ( G % METRIC_F_DD_11 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_DD_22 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_DD_33 ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_F_UU_11 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_UU_22 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_UU_33 ) = CoordinateUnit ( 1 ) ** ( -2 )
    case ( 'SPHERICAL' )
      VariableUnit ( G % VOLUME )  &
        = CoordinateUnit ( 1 ) ** 3
      VariableUnit ( G % AREA_I_D_1 )  &
        =  CoordinateUnit ( 1 ) ** 2
      VariableUnit ( G % AREA_I_D_2 )  &
        =  CoordinateUnit ( 1 ) ** 3
      VariableUnit ( G % AREA_I_D_3 )  &
        =  CoordinateUnit ( 1 ) ** 3
      VariableUnit ( G % METRIC_F_DD_11 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_DD_22 ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_F_DD_33 ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_F_UU_11 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_F_UU_22 ) = CoordinateUnit ( 1 ) ** ( -2 )
      VariableUnit ( G % METRIC_F_UU_33 ) = CoordinateUnit ( 1 ) ** ( -2 )
    end select !-- CoordinateSystem

  end subroutine SetUnits

  
end module Geometry_F__Form
