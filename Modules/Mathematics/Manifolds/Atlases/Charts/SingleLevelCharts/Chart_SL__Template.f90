!-- Chart_SL contains functionality common to local and distributed
!   single-level charts.

module Chart_SL__Template

  !-- Chart_SingleLevel_Template

  use Basics
  use AtlasBasics
  use ChartBasics
  use ChartStream_SL__Form

  implicit none
  private

  type, public, extends ( ChartHeader_SL_Form ), abstract :: &
    Chart_SL_Template
      class ( GeometryFlat_CSL_Form ), pointer :: &
        Geometry_CSL => null ( )
      type ( ChartStream_SL_ElementForm ), dimension ( : ), allocatable :: &
        Stream
  contains
    procedure, public, pass :: &
      SetGeometry
    procedure, public, pass :: &
      Geometry
    procedure, public, pass :: &
      OpenStream
    procedure, public, pass :: &
      AddFieldImage
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      CloseStreams
    procedure, public, pass ( C ) :: &
      CopyBoundaryTemplate
    procedure ( CB ), public, pass ( C ), deferred :: &
      CopyBoundary
    procedure, public, pass ( C ) :: &
      ReverseBoundaryTemplate
    procedure ( CB ), public, pass ( C ), deferred :: &
      ReverseBoundary
    procedure, public, pass :: &
      FinalizeTemplate_CSL
    procedure ( SGWC ), public, pass, deferred :: &
      SetGeometryWidthCenter
  end type Chart_SL_Template

  abstract interface

    subroutine CB ( Value, C, iDimension, iConnection )
      use Basics
      import Chart_SL_Template
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Value
      class ( Chart_SL_Template ), intent ( in ) :: &
        C
      integer ( KDI ), intent ( in ) :: &
        iDimension, &
        iConnection
    end subroutine CB

    subroutine SGWC ( C, Center, Width_L, Width_R, iD )
      use Basics
      import Chart_SL_Template
      class ( Chart_SL_Template ), intent ( inout ) :: &
        C
      type ( Real_1D_Form ), intent ( inout ) :: &
        Center, &
        Width_L, &
        Width_R
      integer ( KDI ), intent ( in ) :: &
        iD  !-- iDimension
    end subroutine SGWC

  end interface

    private :: &
      SetBoundaryLimits, &
      CopyBoundaryKernel, &
      ReverseBoundaryKernel

contains

  
  subroutine SetGeometry ( C, Geometry )

    class ( Chart_SL_Template ), intent ( inout ) :: &
      C
    class ( GeometryFlat_CSL_Form ), intent ( in ), target :: &
      Geometry

    integer ( KDI ) :: &
      iD, & !-- iDimension
      iC  !-- iCell
    type ( Real_1D_Form ), dimension ( ATLAS % MAX_DIMENSIONS ) :: &
      Center, &
      Width_L, &
      Width_R
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Center_3D, &
      Width_L_3D, &
      Width_R_3D
    class ( GeometryFlatForm ), pointer :: &
      G

    C % Geometry_CSL => Geometry

    G => C % Geometry ( )

    do iD = 1, C % nDimensions
      call C % SetGeometryWidthCenter &
             ( Center ( iD ), Width_L ( iD ), Width_R ( iD ), iD )
    end do !-- iD

    call Show ( Center ( 1 ) % Value, C % CoordinateUnit ( 1 ), &
                'Center 1', CONSOLE % INFO_7 )
    call Show ( Width_L ( 1 ) % Value, C % CoordinateUnit ( 1 ), &
                'Width_L 1', CONSOLE % INFO_7 )
    call Show ( Width_R ( 1 ) % Value, C % CoordinateUnit ( 1 ), &
                'Width_R 1', CONSOLE % INFO_7 )
    if ( C % nDimensions > 1 ) then
      call Show ( Center ( 2 ) % Value, C % CoordinateUnit ( 2 ), &
                'Center 2', CONSOLE % INFO_7 )
      call Show ( Width_L ( 2 ) % Value, C % CoordinateUnit ( 2 ), &
                'Width_L 2', CONSOLE % INFO_7 )
      call Show ( Width_R ( 2 ) % Value, C % CoordinateUnit ( 2 ), &
                'Width_R 2', CONSOLE % INFO_7 )
    end if
    if ( C % nDimensions > 2 ) then
      call Show ( Center ( 3 ) % Value, C % CoordinateUnit ( 3 ), &
                'Center 3', CONSOLE % INFO_7 )
      call Show ( Width_L ( 3 ) % Value, C % CoordinateUnit ( 3 ), &
                'Width_L 3', CONSOLE % INFO_7 )
      call Show ( Width_R ( 3 ) % Value, C % CoordinateUnit ( 3 ), &
                'Width_L 3', CONSOLE % INFO_7 )
    end if

    do iD = 1, C % nDimensions
      call C % SetVariablePointer &
             ( G % Value ( :, G % CENTER ( iD ) ), Center_3D )
      call C % SetVariablePointer &
             ( G % Value ( :, G % WIDTH_LEFT ( iD ) ), Width_L_3D )
      call C % SetVariablePointer &
             ( G % Value ( :, G % WIDTH_RIGHT ( iD ) ), Width_R_3D )
      associate &
        ( Center_1D => Center ( iD ) % Value, &
          Width_L_1D  => Width_L ( iD ) % Value, &
          Width_R_1D  => Width_R ( iD ) % Value )
      do iC = C % iaFirst ( iD ), C % iaLast ( iD )
        select case ( iD )
        case ( 1 )
          Center_3D ( iC, :, : ) = Center_1D ( iC )
          Width_L_3D  ( iC, :, : ) = Width_L_1D  ( iC )
          Width_R_3D  ( iC, :, : ) = Width_R_1D  ( iC )
        case ( 2 )
          Center_3D ( :, iC, : ) = Center_1D ( iC )
          Width_L_3D  ( :, iC, : ) = Width_L_1D  ( iC )
          Width_R_3D  ( :, iC, : ) = Width_R_1D  ( iC )
        case ( 3 )
          Center_3D ( :, :, iC ) = Center_1D ( iC )
          Width_L_3D  ( :, :, iC ) = Width_L_1D  ( iC )
          Width_R_3D  ( :, :, iC ) = Width_R_1D  ( iC )
        end select !-- iD
      end do !-- iC
      end associate !-- Center_1D, etc.
    end do !-- iD

    call G % SetMetric ( C % nDimensions, G % nValues, oValue = 0 )

    nullify ( G )
    nullify ( Width_R_3D )
    nullify ( Width_L_3D )
    nullify ( Center_3D )

  end subroutine SetGeometry


  function Geometry ( C ) result ( G )

    class ( Chart_SL_Template ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), pointer :: &
      G
      
    class ( VariableGroupForm ), pointer :: &
      Field

    if ( .not. associated ( C % Geometry_CSL ) ) then
      call Show ( 'Chart Geometry not set', CONSOLE % ERROR )
      call Show ( 'To correct, call Chart_SL % SetGeometry', CONSOLE % ERROR )
      call Show ( 'Chart_SL__Template', 'module', CONSOLE % ERROR )
      call Show ( 'Geometry', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Communicator % Synchronize ( )
      call PROGRAM_HEADER % Abort ( )
    end if
    
    Field => C % Geometry_CSL % Field
    select type ( Field )
    class is ( GeometryFlatForm )
      G => Field
    end select !-- Field

  end function Geometry


  subroutine OpenStream ( C, GIS, Label, iStream, VerboseOption )

    class ( Chart_SL_Template ), intent ( inout ) :: &
      C
    type ( GridImageStreamForm ), intent ( in ) :: &
      GIS
    character ( * ), intent ( in ) :: &
      Label
    integer ( KDI ), intent ( in ) :: &
      iStream
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    integer ( KDI ) :: &
      iF  !-- iField
    logical ( KDL ) :: &
      Verbose

    Verbose = .false.
    if ( present ( VerboseOption ) ) &
      Verbose = VerboseOption
    
    if ( .not. allocated ( C % Stream ) ) &
      allocate ( C % Stream ( ATLAS % MAX_STREAMS ) )

    allocate ( C % Stream ( iStream ) % Element )
    associate ( S => C % Stream ( iStream ) % Element )

    call S % Initialize &
           ( C, GIS, trim ( C % Name ) // '_Stream_' // trim ( Label ) )

    do iF = 1, C % nFields
      select type ( FC => C % Field ( iF ) % Pointer )
      class is ( Field_CSL_Template )
        if ( Verbose ) then
          call S % AddField ( FC % Field )
        else
          call S % AddField ( FC % FieldOutput )
        end if
      end select !-- FC
    end do !-- iF

    end associate !-- S

  end subroutine OpenStream


  subroutine AddFieldImage ( C, Field, iStream )

    class ( Chart_SL_Template ), intent ( inout ), target :: &
      C
    class ( VariableGroupForm ), intent ( in ) :: &
      Field
    integer ( KDI ), intent ( in ) :: &
      iStream

    associate ( CS => C % Stream ( iStream ) % Element )
    call CS % AddField ( Field )
    end associate !-- CS

  end subroutine AddFieldImage


  subroutine Write ( C, iStream, DirectoryOption, TimeOption, &
                     CycleNumberOption )

    class ( Chart_SL_Template ), intent ( inout ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iStream
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption

    associate ( CS => C % Stream ( iStream ) % Element )
    call CS % Write &
           ( DirectoryOption, TimeOption, CycleNumberOption )
    end associate !-- CS
    
  end subroutine Write


  subroutine CloseStreams ( C )

    class ( Chart_SL_Template ), intent ( inout ) :: &
      C

    deallocate ( C % Stream )

  end subroutine CloseStreams


  subroutine CopyBoundaryTemplate &
               ( Value, C, nCells, iDimension, iConnection )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    class ( Chart_SL_Template ), intent ( in ) :: &
      C
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nCells
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    integer ( KDI ), dimension ( 3 ) :: &
      oBI, &  !-- oBoundaryInterior
      oBE, &  !-- oBoundaryExterior
      dBI, &  !-- dBoundaryInterior, i.e. direction
      dBE, &  !-- dBoundaryExterior, i.e. direction
      nB      !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V

    call SetBoundaryLimits &
           ( C, nCells, iDimension, iConnection, nB, dBE, dBI, oBE, oBI )

    call C % SetVariablePointer ( Value, V )

    call CopyBoundaryKernel ( V, nB, dBE, dBI, oBE, oBI )

    nullify ( V )

  end subroutine CopyBoundaryTemplate


  subroutine ReverseBoundaryTemplate &
               ( Value, C, nCells, iDimension, iConnection )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    class ( Chart_SL_Template ), intent ( in ) :: &
      C
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nCells
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    integer ( KDI ), dimension ( 3 ) :: &
      oBI, &  !-- oBoundaryInterior
      oBE, &  !-- oBoundaryExterior
      dBI, &  !-- dBoundaryInterior, i.e. direction
      dBE, &  !-- dBoundaryExterior, i.e. direction
      nB      !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V

    call SetBoundaryLimits &
           ( C, nCells, iDimension, iConnection, nB, dBE, dBI, oBE, oBI )

    call C % SetVariablePointer ( Value, V )

    call ReverseBoundaryKernel ( V, nB, dBE, oBE )

    nullify ( V )

  end subroutine ReverseBoundaryTemplate


  impure elemental subroutine FinalizeTemplate_CSL ( C )

    class ( Chart_SL_Template ), intent ( inout ) :: &
      C

    if ( allocated ( C % Stream ) ) deallocate ( C % Stream )

    nullify ( C % Geometry_CSL )

  end subroutine FinalizeTemplate_CSL


  subroutine SetBoundaryLimits &
               ( C, nCells, iDimension, iConnection, nB, dBE, dBI, oBE, oBI )

    class ( Chart_SL_Template ), intent ( in ) :: &
      C
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nCells
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection
    integer ( KDI ), dimension ( 3 ), intent ( out ) :: &
      nB, &   !-- nBoundary
      dBI, &  !-- dBoundaryInterior, i.e. direction
      dBE, &  !-- dBoundaryExterior, i.e. direction
      oBI, &  !-- oBoundaryInterior
      oBE     !-- oBoundaryExterior

    integer ( KDI ) :: &
      jD, kD   !-- jDimension, kDimension

    associate &
      ( Connectivity => C % Atlas % Connectivity, &
        iD => iDimension, &
        iC => iConnection )

    jD = mod ( iD, 3 ) + 1
    kD = mod ( jD, 3 ) + 1

    !-- In setting oBI and oBE, note kernel routine does not inherit lbound

    oBI = C % nGhostLayers
    dBI = +1
    if ( iC == Connectivity % iaOuter ( iD ) ) then
      oBI ( iD ) = oBI ( iD ) + nCells ( iD ) + 1
      dBI ( iD ) = -1
    end if !-- iC

    oBE = oBI
    dBE = dBI
    dBE ( iD ) = -dBI ( iD )
    if ( iC == Connectivity % iaInner ( iD ) ) then
      oBE ( iD ) = oBE ( iD ) + 1
    else if ( iC == Connectivity % iaOuter ( iD ) ) then
      oBE ( iD ) = oBE ( iD ) - 1
    end if !-- iC

    nB ( iD ) = C % nGhostLayers ( iD )
    nB ( jD ) = nCells ( jD )
    nB ( kD ) = nCells ( kD )

    end associate !-- Connectivity, etc.
    
  end subroutine SetBoundaryLimits


  subroutine CopyBoundaryKernel ( V, nB, dBE, dBI, oBE, oBI )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      V
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nB,  & 
      dBE, &
      dBI, &
      oBE, &
      oBI

    integer ( KDI ) :: &
      iV, jV, kV

    !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
    do kV = 1, nB ( 3 )
      do jV = 1, nB ( 2 )
        do iV = 1, nB ( 1 )
          V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
              oBE ( 2 )  +  dBE ( 2 ) * jV, &
              oBE ( 3 )  +  dBE ( 3 ) * kV ) &
            = V ( oBI ( 1 )  +  dBI ( 1 ) * iV, &
                  oBI ( 2 )  +  dBI ( 2 ) * jV, &
                  oBI ( 3 )  +  dBI ( 3 ) * kV )
        end do 
      end do
    end do
    !$OMP end parallel do

  end subroutine CopyBoundaryKernel


  subroutine ReverseBoundaryKernel ( V, nB, dBE, oBE )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      V
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nB,  & 
      dBE, &
      oBE

    integer ( KDI ) :: &
      iV, jV, kV

    !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
    do kV = 1, nB ( 3 )
      do jV = 1, nB ( 2 )
        do iV = 1, nB ( 1 )
          V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
              oBE ( 2 )  +  dBE ( 2 ) * jV, &
              oBE ( 3 )  +  dBE ( 3 ) * kV ) &
            = - V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
                    oBE ( 2 )  +  dBE ( 2 ) * jV, &
                    oBE ( 3 )  +  dBE ( 3 ) * kV )
        end do 
      end do
    end do
    !$OMP end parallel do

  end subroutine ReverseBoundaryKernel


end module Chart_SL__Template
