!-- GeometryFlat extends VariableGroup to include functionality related to
!   flat geometry, though with possibly curvilinear coordinates.

module GeometryFlat_Form

  use Basics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_FLAT  = 11, &
      N_VECTORS_FLAT = 1

  type, public, extends ( VariableGroupForm ) :: GeometryFlatForm
    integer ( KDI ) :: &
      IGNORABILITY   = 0, &
      N_FIELDS       = 0, &
      N_VECTORS      = 0, &
      N_FIELDS_FLAT  = N_FIELDS_FLAT, &
      N_VECTORS_FLAT = N_VECTORS_FLAT
    integer ( KDI ) :: &
      VOLUME_JACOBIAN = 0
    integer ( KDI ), dimension ( 3 ) :: &
      CENTER = 0, &
      WIDTH  = 0
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
      SetMetric
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass ( G ) :: &
      ComputeReconstruction
    final :: &  !-- FIXME: Intel doesn't like final procedure name to be the
                !          same as the parent's final
      Finalize_G
  end type GeometryFlatForm

    private :: &
      InitializeBasics, &
      SetCoordinateSystem, &
      SetUnits, &
      SetMetricCartesian, &
      SetMetricCylindrical, &
      SetMetricSpherical, &
      ComputeEdges

contains


  subroutine InitializeAllocate_G &
               ( G, CoordinateSystem, CoordinateUnit, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

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
      ClearOption
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

    call G % VariableGroupForm % Initialize &
           ( [ nValues, G % N_FIELDS ], &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_G


  subroutine SetMetric ( G, nDimensions, nValues, oValue )

    class ( GeometryFlatForm ), intent ( inout ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'CARTESIAN' )
      call SetMetricCartesian &
             ( G % Value ( :, G % VOLUME_JACOBIAN ), &
               G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               nDimensions, nValues, oValue )
    case ( 'CYLINDRICAL' )
      call SetMetricCylindrical &
             ( G % Value ( :, G % VOLUME_JACOBIAN ), &
               G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               G % Value ( :, G % CENTER ( 1 ) ), &
               nDimensions, nValues, oValue )
    case ( 'SPHERICAL' )
      call SetMetricSpherical &
             ( G % Value ( :, G % VOLUME_JACOBIAN ), &
               G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               G % Value ( :, G % METRIC_UU_22 ), &
               G % Value ( :, G % METRIC_UU_33 ), &
               G % Value ( :, G % CENTER ( 1 ) ), &
               G % Value ( :, G % CENTER ( 2 ) ), &
               nDimensions, nValues, oValue )
    case default
      call Show ( 'GeometryFlatForm % SetCellGeometry not implemented for ' &
                  // trim ( G % CoordinateSystem ), CONSOLE % ERROR )
    end select

  end subroutine SetMetric


  subroutine SetOutput ( G, Output )

    class ( GeometryFlatForm ), intent ( inout ) :: &
      G
    class ( VariableGroupForm ), intent ( inout ) :: &
      Output

    call Output % Initialize &
           ( G, iaSelectedOption &
                  = [ G % CENTER, G % METRIC_DD_22, G % METRIC_DD_33 ] )

  end subroutine SetOutput


  subroutine ComputeReconstruction ( G_I, G, nDimensions, iDimension )

    type ( VariableGroupForm ), intent ( inout ) :: &
      G_I
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      iDimension

    integer ( KDI ) :: &
      iD

    do iD = 1, 3
      if ( iD == iDimension ) &
        cycle
      call Copy ( G   % Value ( :, G % CENTER ( iD ) ), &
                  G_I % Value ( :, G % CENTER ( iD ) ) )
    end do

    call ComputeEdges &
           ( G   % Value ( :, G % CENTER ( iDimension ) ), &
             G   % Value ( :, G % WIDTH ( iDimension ) ), &
             G_I % Value ( :, G % CENTER ( iDimension ) ) )

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'CARTESIAN' )
      call SetMetricCartesian &
             ( G_I % Value ( :, G % VOLUME_JACOBIAN ), &
               G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               nDimensions, nValues = G % nValues, oValue = 0 )
    case ( 'CYLINDRICAL' )
      call SetMetricCylindrical &
             ( G_I % Value ( :, G % VOLUME_JACOBIAN ), &
               G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               G_I % Value ( :, G % CENTER ( 1 ) ), &
               nDimensions, nValues = G % nValues, oValue = 0 )
    case ( 'SPHERICAL' )
      call SetMetricSpherical &
             ( G_I % Value ( :, G % VOLUME_JACOBIAN ), &
               G_I % Value ( :, G % METRIC_DD_22 ), &
               G_I % Value ( :, G % METRIC_DD_33 ), &
               G_I % Value ( :, G % METRIC_UU_22 ), &
               G_I % Value ( :, G % METRIC_UU_33 ), &
               G_I % Value ( :, G % CENTER ( 1 ) ), &
               G_I % Value ( :, G % CENTER ( 2 ) ), &
               nDimensions, nValues = G % nValues, oValue = 0 )
    end select

  end subroutine ComputeReconstruction


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

    G % VOLUME_JACOBIAN = 1
    G % CENTER          = [ 2, 3, 4 ]
    G % WIDTH           = [ 5, 6, 7 ]
    G % METRIC_DD_22    = 8
    G % METRIC_DD_33    = 9
    G % METRIC_UU_22    = 10
    G % METRIC_UU_33    = 11

    if ( G % N_FIELDS == 0 ) G % N_FIELDS = G % N_FIELDS_FLAT

    !-- variable names

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( G % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : G % N_FIELDS_FLAT ) &
      = [ 'VolumeJacobian', &
          'Center_1      ', &
          'Center_2      ', &
          'Center_3      ', &
          'Width_1       ', &
          'Width_2       ', &
          'Width_3       ', &
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

    if ( G % N_VECTORS == 0 ) G % N_VECTORS = G % N_VECTORS_FLAT

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

    call VectorIndices ( 1 ) % Initialize ( G % CENTER )

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

    VariableUnit ( G % CENTER ( 1 ) : G % CENTER ( 3 ) ) &
      = CoordinateUnit
    VariableUnit ( G % WIDTH ( 1 ) : G % WIDTH ( 3 ) ) &
      = CoordinateUnit

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'CARTESIAN' )
      VariableUnit ( G % VOLUME_JACOBIAN ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_DD_22    ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_DD_33    ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_UU_22    ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_UU_33    ) = UNIT % IDENTITY
    case ( 'CYLINDRICAL' )
      VariableUnit ( G % VOLUME_JACOBIAN ) = CoordinateUnit ( 1 )
      VariableUnit ( G % METRIC_DD_22    ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_DD_33    ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_UU_22    ) = UNIT % IDENTITY
      VariableUnit ( G % METRIC_UU_33    ) = CoordinateUnit ( 1 ) ** ( -2 )
    case ( 'SPHERICAL' )
      VariableUnit ( G % VOLUME_JACOBIAN ) = CoordinateUnit ( 1 ) ** 2
      VariableUnit ( G % METRIC_DD_22    ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_DD_33    ) = CoordinateUnit ( 1 ) ** (  2 )
      VariableUnit ( G % METRIC_UU_22    ) = CoordinateUnit ( 1 ) ** ( -2 )
      VariableUnit ( G % METRIC_UU_33    ) = CoordinateUnit ( 1 ) ** ( -2 )
    end select !-- CoordinateSystem

  end subroutine SetUnits


  subroutine SetMetricCartesian &
               ( VJ, M_DD_22, M_DD_33, M_UU_22, M_UU_33, &
                 nDimensions, nValues, oValue )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      VJ, &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue

    integer ( KDI ) :: &
      iV  !-- iValue

    !$OMP parallel do private ( iV )
    do iV = oValue + 1, oValue + nValues
      VJ ( iV ) = 1.0_KDR
      M_DD_22   ( iV ) = 1.0_KDR
      M_DD_33   ( iV ) = 1.0_KDR
      M_UU_22   ( iV ) = 1.0_KDR
      M_UU_33   ( iV ) = 1.0_KDR
    end do
    !$OMP end parallel do

  end subroutine SetMetricCartesian


  subroutine SetMetricCylindrical &
               ( VJ, M_DD_22, M_DD_33, M_UU_22, M_UU_33, X_1, &
                 nDimensions, nValues, oValue )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      VJ, &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_1
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue

    integer ( KDI ) :: &
      iV  !-- iValue

    if ( nDimensions == 1 .or. nDimensions == 2 ) then
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        VJ ( iV )  =  2.0_KDR  *  CONSTANT % PI  *  X_1 ( iV )
      end do
      !$OMP end parallel do
    else if ( nDimensions == 3 ) then
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        VJ ( iV )  =  X_1 ( iV )
      end do
      !$OMP end parallel do
    end if

    !$OMP parallel do private ( iV )
    do iV = oValue + 1, oValue + nValues
      M_DD_22 ( iV )  =  1.0_KDR
      M_DD_33 ( iV )  =  X_1 ( iV ) ** 2 
      M_UU_22 ( iV )  =  1.0_KDR
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = oValue + 1, oValue + nValues
      if ( X_1 ( iV )  >  0.0_KDR ) then
        M_UU_33 ( iV )  =  1.0_KDR  /  X_1 ( iV ) ** 2
      else
        M_UU_33 ( iV )  =  0.0_KDR
      end if
    end do
    !$OMP end parallel do

  end subroutine SetMetricCylindrical


  subroutine SetMetricSpherical &
               ( VJ, M_DD_22, M_DD_33, M_UU_22, M_UU_33, X_1, X_2, &
                 nDimensions, nValues, oValue )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      VJ, &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_1, X_2
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nValues, &
      oValue

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ), dimension ( size ( VJ ) ) :: &
      R_Squared, &
      SinTheta

    !$OMP parallel do private ( iV )
    do iV = oValue + 1, oValue + nValues
      R_Squared ( iV )  =  X_1 ( iV ) ** 2
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = oValue + 1, oValue + nValues
      !-- Vanish exactly on axis
      if ( abs ( sin ( X_2 ( iV ) ) ) > 1.0e-10 ) then
        SinTheta ( iV )  =  sin ( X_2 ( iV ) )
      else
        SinTheta ( iV )  =  0.0_KDR
      end if
    end do
    !$OMP end parallel do

    if ( nDimensions == 1 ) then
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        VJ ( iV )  =  4.0_KDR  *  CONSTANT % PI  *  R_Squared ( iV )
      end do
      !$OMP end parallel do
    else if ( nDimensions == 2 ) then
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        VJ ( iV ) &
          =  2.0_KDR  *  CONSTANT % PI  *  R_Squared ( iV )  *  SinTheta ( iV )
      end do
      !$OMP end parallel do
    else if ( nDimensions == 3 ) then
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        VJ ( iV )  =  R_Squared ( iV )  *  SinTheta ( iV )
      end do
      !$OMP end parallel do
    end if

    if ( nDimensions == 1 ) then
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  R_Squared ( iV )
        M_DD_33 ( iV )  =  R_Squared ( iV )
      end do
      !$OMP end parallel do
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        if ( R_Squared ( iV )  >  0.0_KDR ) then
          M_UU_22 ( iV )  =  1.0_KDR  /  R_Squared ( iV )
          M_UU_33 ( iV )  =  1.0_KDR  /  R_Squared ( iV )
        else
          M_UU_22 ( iV )  =  0.0_KDR
          M_UU_33 ( iV )  =  0.0_KDR
        end if
      end do
      !$OMP end parallel do
    else
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  R_Squared ( iV )
        M_DD_33 ( iV )  =  R_Squared ( iV )  *  SinTheta ( iV ) ** 2
      end do
      !$OMP end parallel do
      !$OMP parallel do private ( iV )
      do iV = oValue + 1, oValue + nValues
        if ( R_Squared ( iV )  *  SinTheta ( iV ) ** 2  >  0.0_KDR ) then
          M_UU_22 ( iV ) &
            =  1.0_KDR  /  R_Squared ( iV )
          M_UU_33 ( iV ) &
            =  1.0_KDR  /  ( R_Squared ( iV )  *  SinTheta ( iV ) ** 2 )
        else
          M_UU_22 ( iV )  =  0.0_KDR
          M_UU_33 ( iV )  =  0.0_KDR
        end if
      end do
      !$OMP end parallel do
    end if

  end subroutine SetMetricSpherical


  subroutine ComputeEdges ( X, dX, X_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
       X, &
      dX
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      X_I

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( X )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      X_I ( iV )  =  X ( iV ) - 0.5_KDR * dX ( iV )
    end do
    !$OMP end parallel do
    
  end subroutine ComputeEdges


end module GeometryFlat_Form
