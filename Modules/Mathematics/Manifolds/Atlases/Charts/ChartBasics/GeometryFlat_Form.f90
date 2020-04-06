!-- GeometryFlat extends Storage to include functionality related to
!   flat geometry, though with possibly curvilinear coordinates.

module GeometryFlat_Form
  
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
    integer ( KDI ) :: &
      iTimerSetFiniteVolume, &
      iTimerSetMetric, &
      iTimerCopyGeometry, &
      iTimerComputeEdges
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
      SetUnits!, &
!    FIXME: required to generate .smod with GCC 6.1.0  
    public :: &
      SetFiniteVolumeRectangularKernel, &
      SetFiniteVolumeCylindricalKernel, &
      SetFiniteVolumeSphericalKernel, &
      SetMetricRectangularKernel, &
      SetMetricCylindricalKernel, &
      SetMetricSphericalKernel, &
      ComputeEdgesKernel
    
    interface
      
      module subroutine SetFiniteVolumeRectangularKernel &
                   ( V, A_I_1, A_I_2, A_I_3, &
                     W_L_1, W_L_2, W_L_3, W_R_1, W_R_2, W_R_3, &
                     nDimensions, nValues, oValue )
        use Basics
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
      end subroutine SetFiniteVolumeRectangularKernel

      module subroutine SetFiniteVolumeCylindricalKernel &
                   ( V, A_I_1, A_I_2, A_I_3, &
                     W_L_1, W_L_2, W_L_3, W_R_1, W_R_2, W_R_3, RP_C, &
                     nDimensions, nValues, oValue )
        use Basics
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
      end subroutine SetFiniteVolumeCylindricalKernel


      module subroutine SetFiniteVolumeSphericalKernel &
                   ( V, A_I_1, A_I_2, A_I_3, &
                     W_L_1, W_L_2, W_L_3, W_R_1, W_R_2, W_R_3, R_C, Th_C, &
                     nDimensions, nValues, oValue )
        use Basics
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
      end subroutine SetFiniteVolumeSphericalKernel


      module subroutine SetMetricRectangularKernel &
                   ( M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
                     nDimensions, nValues, oValue, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M_DD_22, M_DD_33, &
          M_UU_22, M_UU_33
        integer ( KDI ), intent ( in ) :: &
          nDimensions, &
          nValues, &
          oValue
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine SetMetricRectangularKernel


      module subroutine SetMetricCylindricalKernel &
                   ( M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
                     RP, nDimensions, nValues, oValue, UseDeviceOption )
        use Basics
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
      end subroutine SetMetricCylindricalKernel


      module subroutine SetMetricSphericalKernel &
                   ( M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
                     R, Th, nDimensions, nValues, oValue, UseDeviceOption )
        use Basics
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
      end subroutine SetMetricSphericalKernel
      

      module subroutine ComputeEdgesKernel ( X, dX_L, X_I, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          X, &
          dX_L
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          X_I
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeEdgesKernel
    
    end interface

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
    
    call PROGRAM_HEADER % AddTimer &
           ( 'SetMetricFixed', G % iTimerSetMetric, Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'SetFiniteVolume', G % iTimerSetFiniteVolume, Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'CopyGeometry', G % iTimerCopyGeometry, Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeEdges', G % iTimerComputeEdges, Level = 6 )

  end subroutine InitializeAllocate_G


  subroutine SetMetricFixed ( G, nDimensions, nValues, oValue )

    class ( GeometryFlatForm ), intent ( inout ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue
      
    associate & 
      ( T_FV => PROGRAM_HEADER % Timer ( G % iTimerSetFiniteVolume ), &
        T_SM => PROGRAM_HEADER % Timer ( G % iTimerSetMetric ) )

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call T_FV % Start ( )
      call SetFiniteVolumeRectangularKernel &
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
      call T_FV % Stop ( )
      call T_SM % Start ( )
      call SetMetricRectangularKernel &
             ( G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               nDimensions, nValues, oValue )
      call T_SM % Stop ( )
    case ( 'CYLINDRICAL' )
      call T_FV % Start ( )
      call SetFiniteVolumeCylindricalKernel &
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
      call T_FV % Stop ( )
      call T_SM % Start ( )
      call SetMetricCylindricalKernel &
             ( G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               nDimensions, nValues, oValue )
      call T_SM % Stop ( )
    case ( 'SPHERICAL' )
      call T_FV % Start ( )
      call SetFiniteVolumeSphericalKernel &
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
      call T_FV % Stop ( )
      call T_SM % Start ( )
      call SetMetricSphericalKernel &
             ( G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               G % Value ( :, G % CENTER_U ( 2 ) ), &
               nDimensions, nValues, oValue )
      call T_SM % Stop ( )
    case default
      call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
      call Show ( G % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'GeometryFlat_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetMetricFixed', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select
    
    end associate !-- T_FV, T_SM

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
    
    
    associate & 
      ( T_C   => PROGRAM_HEADER % Timer ( G % iTimerCopyGeometry ), &
        T_CE  => PROGRAM_HEADER % Timer ( G % iTimerComputeEdges ), & 
        T_SM  => PROGRAM_HEADER % Timer ( G % iTimerSetMetric ) )  
        
    call T_C % Start ( )

    do iD = 1, 3
      if ( iD == iDimension ) &
        cycle
      call Copy ( G   % Value ( :, G % CENTER_U ( iD ) ), &
                  G_I % Value ( :, G % CENTER_U ( iD ) ), &
                  UseDeviceOption = G % AllocatedDevice )
    end do
    
    call T_C % Stop ( )
    
    call T_CE % Start ( )
    call ComputeEdgesKernel &
           ( G   % Value ( :, G % CENTER_U ( iDimension ) ), &
             G   % Value ( :, G % WIDTH_LEFT_U ( iDimension ) ), &
             G_I % Value ( :, G % CENTER_U ( iDimension ) ), &
             UseDeviceOption = G % AllocatedDevice )
    call T_CE % Stop ( )
    
    call T_SM % Start ( )
    
    select case ( trim ( G % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call SetMetricRectangularKernel &
             ( G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               nDimensions, nValues = G % nValues, oValue = 0, &
               UseDeviceOption = G % AllocatedDevice )
    case ( 'CYLINDRICAL' )
      call SetMetricCylindricalKernel &
             ( G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               G_I % Value ( :, G % CENTER_U ( 1 ) ), &
               nDimensions, nValues = G % nValues, oValue = 0, &
               UseDeviceOption = G % AllocatedDevice )
    case ( 'SPHERICAL' )
      call SetMetricSphericalKernel &
             ( G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               G_I % Value ( :, G % CENTER_U ( 1 ) ), &
               G_I % Value ( :, G % CENTER_U ( 2 ) ), &
               nDimensions, nValues = G % nValues, oValue = 0, &
               UseDeviceOption = G % AllocatedDevice )
    end select
    
    call T_SM % Stop ( )
    
    end associate !-- T_C, T_CE, T_SM

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

  
end module GeometryFlat_Form
