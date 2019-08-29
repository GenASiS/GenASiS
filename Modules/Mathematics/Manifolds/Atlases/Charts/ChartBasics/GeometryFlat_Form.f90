!-- GeometryFlat extends Storage to include functionality related to
!   flat geometry, though with possibly curvilinear coordinates.

#include "Preprocessor"

module GeometryFlat_Form
  
  use iso_c_binding
  use Basics
  use ChartHeader_SL__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_FLAT  = 20, &
      N_VECTORS_FLAT = 1

  type, public, extends ( StorageForm ) :: GeometryFlatForm
    integer ( KDI ) :: &
      IGNORABILITY   = 0, &
      N_FIELDS       = 0, &
      N_VECTORS      = 0, &
      N_FIELDS_FLAT  = N_FIELDS_FLAT, &
      N_VECTORS_FLAT = N_VECTORS_FLAT
    integer ( KDI ) :: &
      VOLUME = 0
    integer ( KDI ), dimension ( 3 ) :: &
      CENTER_U      = 0, &
      WIDTH_LEFT_U  = 0, &
      WIDTH_RIGHT_U = 0, &
      AREA_INNER_D  = 0, &
      COARSENING    = 0
    integer ( KDI ) :: &
      METRIC_DD_22 = 0, &
      METRIC_DD_33 = 0, &
      METRIC_UU_22 = 0, &
      METRIC_UU_33 = 0
    character ( LDL ) :: &
      CoordinateSystem = ''
    character ( LDF ) :: &
      Type = ''
  contains
    procedure, public, pass :: &
      InitializeAllocate_G
    generic, public :: &
      Initialize => InitializeAllocate_G
    procedure, public, pass :: &
      SetMetricFixed
    procedure, public, pass :: &
      SetOutput
    procedure, private, pass ( G ) :: &
      ComputeReconstruction_CSL
    generic, public :: &
      ComputeReconstruction => ComputeReconstruction_CSL
    final :: &  !-- FIXME: Intel doesn't like final procedure name to be the
                !          same as the parent's final
      Finalize_G
  end type GeometryFlatForm

    private :: &
      InitializeBasics, &
      SetCoordinateSystem, &
      SetUnits, &
      SetFiniteVolumeRectangular, &
      SetFiniteVolumeCylindrical, &
      SetFiniteVolumeSpherical, &
      SetMetricRectangular, &
      SetMetricCylindrical, &
      SetMetricSpherical, &
      ComputeEdges

contains


  subroutine InitializeAllocate_G &
               ( G, CoordinateSystem, CoordinateUnit, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 PinnedOption, UnitOption, VectorIndicesOption )

    class ( GeometryFlatForm ), intent ( inout ) :: &
      G
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      CoordinateUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption, &
      PinnedOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDF ) :: &
      Name 
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector
    logical ( KDL ) :: &
      Clear

    call InitializeBasics &
           ( G, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call SetCoordinateSystem ( G, CoordinateSystem )

    call SetUnits ( VariableUnit, G, CoordinateUnit )

    Clear = .true.
    if ( present ( ClearOption ) ) Clear = ClearOption

    call G % StorageForm % Initialize &
           ( [ nValues, G % N_FIELDS ], &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = Clear, &
             PinnedOption = PinnedOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

    G % Value ( :, G % COARSENING ( 1 ) )  =  1.0_KDR
    G % Value ( :, G % COARSENING ( 2 ) )  =  1.0_KDR
    G % Value ( :, G % COARSENING ( 3 ) )  =  1.0_KDR

  end subroutine InitializeAllocate_G


  subroutine SetMetricFixed ( G, nDimensions, nValues, oValue )

    class ( GeometryFlatForm ), intent ( inout ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call SetFiniteVolumeRectangular &
             ( G % Value ( :, G % VOLUME ), &
               G % Value ( :, G % AREA_INNER_D ( 1 ) ), &
               G % Value ( :, G % AREA_INNER_D ( 2 ) ), &
               G % Value ( :, G % AREA_INNER_D ( 3 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
               nDimensions, nValues, oValue )
      call SetMetricRectangular &
             ( G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               nDimensions, nValues, oValue )
    case ( 'CYLINDRICAL' )
      call SetFiniteVolumeCylindrical &
             ( G % Value ( :, G % VOLUME ), &
               G % Value ( :, G % AREA_INNER_D ( 1 ) ), &
               G % Value ( :, G % AREA_INNER_D ( 2 ) ), &
               G % Value ( :, G % AREA_INNER_D ( 3 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               nDimensions, nValues, oValue )
      call SetMetricCylindrical &
             ( G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               nDimensions, nValues, oValue )
    case ( 'SPHERICAL' )
      call SetFiniteVolumeSpherical &
             ( G % Value ( :, G % VOLUME ), &
               G % Value ( :, G % AREA_INNER_D ( 1 ) ), &
               G % Value ( :, G % AREA_INNER_D ( 2 ) ), &
               G % Value ( :, G % AREA_INNER_D ( 3 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               G % Value ( :, G % CENTER_U ( 2 ) ), &
               nDimensions, nValues, oValue )
      call SetMetricSpherical &
             ( G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               G % Value ( :, G % CENTER_U ( 2 ) ), &
               nDimensions, nValues, oValue )
    case default
      call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
      call Show ( G % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'GeometryFlat_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetMetricFixed', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine SetMetricFixed


  subroutine SetOutput ( G, Output )

    class ( GeometryFlatForm ), intent ( inout ) :: &
      G
    class ( StorageForm ), intent ( inout ) :: &
      Output

    call Output % Initialize &
           ( G, iaSelectedOption &
                  = [ G % CENTER_U, G % METRIC_DD_22, G % METRIC_DD_33 ] )

  end subroutine SetOutput


  subroutine ComputeReconstruction_CSL ( G_I, CSL, G, nDimensions, iDimension )

    type ( StorageForm ), intent ( inout ) :: &
      G_I
    class ( ChartHeader_SL_Form ), intent ( in ) :: &
      CSL
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      iDimension

    integer ( KDI ) :: &
      iD
      
    if ( G % AllocatedDevice .neqv. G_I % AllocatedDevice ) then
      call Show ( 'Geometry not allocated consistently on Host or Device', &
                  CONSOLE % ERROR )
      call Show ( 'GeometryFlat_Form', 'module', &
                  CONSOLE % ERROR )
      call Show ( 'ComputeReconstruction_CSL', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    do iD = 1, 3
      if ( iD == iDimension ) &
        cycle
      call Copy ( G   % Value ( :, G % CENTER_U ( iD ) ), &
                  G_I % Value ( :, G % CENTER_U ( iD ) ), &
                  UseDeviceOption = G % AllocatedDevice )
    end do
    
    call ComputeEdges &
           ( G   % Value ( :, G % CENTER_U ( iDimension ) ), &
             G   % Value ( :, G % WIDTH_LEFT_U ( iDimension ) ), &
             G_I % Value ( :, G % CENTER_U ( iDimension ) ), &
             UseDeviceOption = G % AllocatedDevice )

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call SetMetricRectangular &
             ( G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               nDimensions, nValues = G % nValues, oValue = 0, &
               UseDeviceOption = G % AllocatedDevice )
    case ( 'CYLINDRICAL' )
      call SetMetricCylindrical &
             ( G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               G_I % Value ( :, G % CENTER_U ( 1 ) ), &
               nDimensions, nValues = G % nValues, oValue = 0, &
               UseDeviceOption = G % AllocatedDevice )
    case ( 'SPHERICAL' )
      call SetMetricSpherical &
             ( G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               G_I % Value ( :, G % CENTER_U ( 1 ) ), &
               G_I % Value ( :, G % CENTER_U ( 2 ) ), &
               nDimensions, nValues = G % nValues, oValue = 0, &
               UseDeviceOption = G % AllocatedDevice )
    end select

  end subroutine ComputeReconstruction_CSL


  impure elemental subroutine Finalize_G ( G )

    type ( GeometryFlatForm ), intent ( inout ) :: &
      G

    if ( G % Name == '' ) return

    call Show ( 'Finalizing a ' // trim ( G % Type ), G % IGNORABILITY )
    call Show ( G % Name, 'Name', G % IGNORABILITY )

  end subroutine Finalize_G


  subroutine InitializeBasics &
               ( G, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( GeometryFlatForm ), intent ( inout ) :: &
      G
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable, &
      Vector
    character ( LDF ), intent ( out ) :: &
      Name
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
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      iV  !-- iVector

    if ( G % Type == '' ) &
      G % Type = 'GeometryFlat'

    Name = 'Geometry'
    if ( present ( NameOption ) ) Name = NameOption

    G % IGNORABILITY = CONSOLE % INFO_4
    call Show ( 'Initializing a ' // trim ( G % Type ), G % IGNORABILITY )
    call Show ( Name, 'Name', G % IGNORABILITY )

    !-- variable indices

    G % VOLUME         =  1
    G % CENTER_U       =  [  2,  3,  4 ]
    G % WIDTH_LEFT_U   =  [  5,  6,  7 ]
    G % WIDTH_RIGHT_U  =  [  8,  9, 10 ]
    G % AREA_INNER_D   =  [ 11, 12, 13 ]
    G % COARSENING     =  [ 14, 15, 16 ]
    G % METRIC_DD_22   =  17
    G % METRIC_DD_33   =  18
    G % METRIC_UU_22   =  19
    G % METRIC_UU_33   =  20

    if ( G % N_FIELDS == 0 ) &
      G % N_FIELDS = G % N_FIELDS_FLAT

    !-- variable names

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( G % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : G % N_FIELDS_FLAT ) &
      = [ 'Volume        ', &
          'Center_U_1    ', &
          'Center_U_2    ', &
          'Center_U_3    ', &
          'WidthLeft_U_1 ', &
          'WidthLeft_U_2 ', &
          'WidthLeft_U_3 ', &
          'WidthRight_U_1', &
          'WidthRight_U_2', &
          'WidthRight_U_3', &
          'AreaInner_D_1 ', &
          'AreaInner_D_2 ', &
          'AreaInner_D_3 ', &
          'Coarsening_1  ', &
          'Coarsening_2  ', &
          'Coarsening_3  ', &
          'Metric_DD_22  ', &
          'Metric_DD_33  ', &
          'Metric_UU_22  ', &
          'Metric_UU_33  ' ]

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

    Vector ( 1 : G % N_VECTORS_FLAT ) &
      = [ 'Position' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = G % N_VECTORS_FLAT + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( G % N_VECTORS ) )
    end if

    call VectorIndices ( 1 ) % Initialize ( G % CENTER_U )

  end subroutine InitializeBasics


  subroutine SetCoordinateSystem ( G, CoordinateSystem )

    class ( GeometryFlatForm ), intent ( inout ) :: &
      G
    character ( * ), intent ( in ) :: &
      CoordinateSystem

    call Show ( 'Setting Geometry coordinate system ', G % IGNORABILITY )
    call Show ( trim ( G % Name ), 'Name', G % IGNORABILITY )

    G % CoordinateSystem = trim ( CoordinateSystem )
    call Show ( trim ( G % CoordinateSystem ), 'Coordinate system', &
                G % IGNORABILITY )

  end subroutine SetCoordinateSystem


  subroutine SetUnits ( VariableUnit, G, CoordinateUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      CoordinateUnit

    VariableUnit ( G % CENTER_U ( 1 ) : G % CENTER_U ( 3 ) ) &
      = CoordinateUnit
    VariableUnit ( G % WIDTH_LEFT_U ( 1 ) : G % WIDTH_LEFT_U ( 3 ) ) &
      = CoordinateUnit
    VariableUnit ( G % WIDTH_RIGHT_U ( 1 ) : G % WIDTH_RIGHT_U ( 3 ) ) &
      = CoordinateUnit

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      VariableUnit ( G % VOLUME )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )  &
           *  CoordinateUnit ( 3 )
      VariableUnit ( G % AREA_INNER_D ( 1 ) )  &
        =  CoordinateUnit ( 2 )  *  CoordinateUnit ( 3 )
      VariableUnit ( G % AREA_INNER_D ( 2 ) )  &
        =  CoordinateUnit ( 3 )  *  CoordinateUnit ( 1 )
      VariableUnit ( G % AREA_INNER_D ( 3 ) )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
      VariableUnit ( G % METRIC_DD_22 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_DD_33 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_UU_22 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_UU_33 ) = UNIT % IDENTITY
    case ( 'CYLINDRICAL' )
      VariableUnit ( G % VOLUME )  &
        =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
      VariableUnit ( G % AREA_INNER_D ( 1 ) )  &
        =  CoordinateUnit ( 1 )  *  CoordinateUnit ( 2 )
      VariableUnit ( G % AREA_INNER_D ( 2 ) )  &
        =  CoordinateUnit ( 1 ) ** 2
      VariableUnit ( G % AREA_INNER_D ( 3 ) )  &
        =  CoordinateUnit ( 1 ) ** 2  *  CoordinateUnit ( 2 )
      VariableUnit ( G % METRIC_DD_22 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_DD_33 ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_UU_22 ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_UU_33 ) = CoordinateUnit ( 1 ) ** ( -2 )
    case ( 'SPHERICAL' )
      VariableUnit ( G % VOLUME )  &
        = CoordinateUnit ( 1 ) ** 3
      VariableUnit ( G % AREA_INNER_D ( 1 ) )  &
        =  CoordinateUnit ( 1 ) ** 2
      VariableUnit ( G % AREA_INNER_D ( 2 ) )  &
        =  CoordinateUnit ( 1 ) ** 3
      VariableUnit ( G % AREA_INNER_D ( 3 ) )  &
        =  CoordinateUnit ( 1 ) ** 3
      VariableUnit ( G % METRIC_DD_22 ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_DD_33 ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_UU_22 ) = CoordinateUnit ( 1 ) ** ( -2 )
      VariableUnit ( G % METRIC_UU_33 ) = CoordinateUnit ( 1 ) ** ( -2 )
    end select !-- CoordinateSystem

  end subroutine SetUnits


  subroutine SetFiniteVolumeRectangular &
               ( V, A_I_1, A_I_2, A_I_3, &
                 W_L_1, W_L_2, W_L_3, W_R_1, W_R_2, W_R_3, &
                 nDimensions, nValues, oValue )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      V, &
      A_I_1, A_I_2, A_I_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      W_L_1, W_L_2, W_L_3, &
      W_R_1, W_R_2, W_R_3
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      dX, dY, dZ

    !$OMP parallel do private ( iV, dX, dY, dZ )
    do iV = oValue + 1, oValue + nValues

      dX  =  W_L_1 ( iV )  +  W_R_1 ( iV )
      dY  =  W_L_2 ( iV )  +  W_R_2 ( iV )
      dZ  =  W_L_3 ( iV )  +  W_R_3 ( iV )

      select case ( nDimensions )
      case ( 1 )
        V ( iV )      =  dX
        A_I_1 ( iV )  =  1.0_KDR
        A_I_2 ( iV )  =  dX
        A_I_3 ( iV )  =  dX
      case ( 2 )
        V ( iV )      =  dX * dY
        A_I_1 ( iV )  =  dY
        A_I_2 ( iV )  =  dX
        A_I_3 ( iV )  =  dX * dY
      case ( 3 )
        V ( iV )      =  dX * dY * dZ
        A_I_1 ( iV )  =  dY * dZ
        A_I_2 ( iV )  =  dZ * dX
        A_I_3 ( iV )  =  dX * dY
      end select

    end do
    !$OMP end parallel do

  end subroutine SetFiniteVolumeRectangular


  subroutine SetFiniteVolumeCylindrical &
               ( V, A_I_1, A_I_2, A_I_3, &
                 W_L_1, W_L_2, W_L_3, W_R_1, W_R_2, W_R_3, RP_C, &
                 nDimensions, nValues, oValue )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      V, &
      A_I_1, A_I_2, A_I_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      W_L_1, W_L_2, W_L_3, &
      W_R_1, W_R_2, W_R_3, &
      RP_C
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      Pi, &
      RP_I, RP_O, &
      dZ, dPh

    Pi  =  CONSTANT % PI

    !$OMP parallel do private ( iV, RP_I, RP_O, dZ, dPh )
    do iV = oValue + 1, oValue + nValues

      RP_I  =  RP_C ( iV )  -  W_L_1 ( iV )
      RP_O  =  RP_C ( iV )  +  W_R_1 ( iV )

        dZ  =  W_L_2 ( iV )  +  W_R_2 ( iV ) 
       dPh  =  W_L_3 ( iV )  +  W_R_3 ( iV ) 

      select case ( nDimensions )
      case ( 1 )
        V ( iV )      =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )  
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  RP_I  
        A_I_2 ( iV )  =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )
      case ( 2 )
        V ( iV )      =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dZ
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  RP_I  *  dZ  
        A_I_2 ( iV )  =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 ) * dZ
      case ( 3 )
        V ( iV )      =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dZ * dPh
        A_I_1 ( iV )  =  RP_I * dZ * dPh  
        A_I_2 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dPh
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 ) * dZ
      end select

    end do
    !$OMP end parallel do

  end subroutine SetFiniteVolumeCylindrical


  subroutine SetFiniteVolumeSpherical &
               ( V, A_I_1, A_I_2, A_I_3, &
                 W_L_1, W_L_2, W_L_3, W_R_1, W_R_2, W_R_3, R_C, Th_C, &
                 nDimensions, nValues, oValue )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      V, &
      A_I_1, A_I_2, A_I_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      W_L_1, W_L_2, W_L_3, &
      W_R_1, W_R_2, W_R_3, &
      R_C, Th_C
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      Pi, &
      R_I, R_O, &
      Th_I, Th_O, &
      dPh

    Pi  =  CONSTANT % PI

    !$OMP parallel do private ( iV, R_I, R_O, Th_I, Th_O, dPh )
    do iV = oValue + 1, oValue + nValues

      R_I  =  R_C ( iV )  -  W_L_1 ( iV )
      R_O  =  R_C ( iV )  +  W_R_1 ( iV )

      Th_I  =  Th_C ( iV )  -  W_L_2 ( iV )
      Th_O  =  Th_C ( iV )  +  W_R_2 ( iV )

      dPh  =  W_L_3 ( iV )  +  W_R_3 ( iV ) 

      select case ( nDimensions )
      case ( 1 )
        V ( iV )      =  4.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 )
        A_I_1 ( iV )  =  4.0_KDR  *  Pi  *  R_I ** 2
        A_I_2 ( iV )  =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 )
        A_I_3 ( iV )  =  2.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )
      case ( 2 )
        V ( iV )      =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  R_I ** 2  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
        A_I_2 ( iV )  =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  sin ( Th_I )
        A_I_3 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
      case ( 3 )
        V ( iV )      =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )  *  dPh
        A_I_1 ( iV )  =  R_I ** 2  *  ( cos ( Th_I )  -  cos ( Th_O ) )  * dPh 
        A_I_2 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  sin ( Th_I )  *  dPh
        A_I_3 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
      end select

    end do
    !$OMP end parallel do

  end subroutine SetFiniteVolumeSpherical


  subroutine SetMetricRectangular &
               ( M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
                 nDimensions, nValues, oValue, UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV  !-- iValue
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  1.0_KDR
        M_UU_22 ( iV )  =  1.0_KDR
        M_UU_33 ( iV )  =  1.0_KDR
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  1.0_KDR
        M_UU_22 ( iV )  =  1.0_KDR
        M_UU_33 ( iV )  =  1.0_KDR
      end do
      !$OMP end parallel do

    end if

  end subroutine SetMetricRectangular


  subroutine SetMetricCylindrical &
               ( M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
                 RP, nDimensions, nValues, oValue, UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      RP
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV  !-- iValue
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  RP ( iV ) ** 2 
        M_UU_22 ( iV )  =  1.0_KDR
        if ( RP ( iV )  >  0.0_KDR ) then
          M_UU_33 ( iV )  =  1.0_KDR  /  RP ( iV ) ** 2
        else
          M_UU_33 ( iV )  =  0.0_KDR
        end if
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  RP ( iV ) ** 2 
        M_UU_22 ( iV )  =  1.0_KDR
        if ( RP ( iV )  >  0.0_KDR ) then
          M_UU_33 ( iV )  =  1.0_KDR  /  RP ( iV ) ** 2
        else
          M_UU_33 ( iV )  =  0.0_KDR
        end if
      end do
      !$OMP end parallel do

    end if

  end subroutine SetMetricCylindrical


  subroutine SetMetricSpherical &
               ( M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
                 R, Th, nDimensions, nValues, oValue, UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, Th
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      Sin_Th
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then      
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, Sin_Th )
      do iV = oValue + 1, oValue + nValues

        select case ( nDimensions )
        case ( 1 )
          Sin_Th  =  1.0_KDR
        case ( 2 )
          Sin_Th  =  sin ( Th ( iV ) )
        case ( 3 )
          Sin_Th  =  sin ( Th ( iV ) )
        end select

        M_DD_22 ( iV )  =  R ( iV ) ** 2
        M_DD_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** 2
        if ( R ( iV )  *  Sin_Th  >  0.0_KDR ) then
          M_UU_22 ( iV )  =  R ( iV ) ** ( -2 )
          M_UU_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** ( -2 )
        else
          M_UU_22 ( iV )  =  0.0_KDR
          M_UU_33 ( iV )  =  0.0_KDR
        end if

      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV, Sin_Th )
      do iV = oValue + 1, oValue + nValues

        select case ( nDimensions )
        case ( 1 )
          Sin_Th  =  1.0_KDR
        case ( 2 )
          Sin_Th  =  sin ( Th ( iV ) )
        case ( 3 )
          Sin_Th  =  sin ( Th ( iV ) )
        end select

        M_DD_22 ( iV )  =  R ( iV ) ** 2
        M_DD_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** 2
        if ( R ( iV )  *  Sin_Th  >  0.0_KDR ) then
          M_UU_22 ( iV )  =  R ( iV ) ** ( -2 )
          M_UU_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** ( -2 )
        else
          M_UU_22 ( iV )  =  0.0_KDR
          M_UU_33 ( iV )  =  0.0_KDR
        end if

      end do
      !$OMP end parallel do
      
    end if

  end subroutine SetMetricSpherical
  

  subroutine ComputeEdges ( X, dX_L, X_I, UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, &
      dX_L
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      X_I
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nV = size ( X )
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do private ( iV ) &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, nV
        X_I ( iV )  =  X ( iV )  -  dX_L ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else      
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nV
        X_I ( iV )  =  X ( iV )  -  dX_L ( iV )
      end do
      !$OMP end parallel do
    end if
    
  end subroutine ComputeEdges

  
end module GeometryFlat_Form
