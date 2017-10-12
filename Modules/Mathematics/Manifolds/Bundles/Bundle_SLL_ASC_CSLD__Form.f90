!-- Bundle_SLL_ASC_CSLD represents a Bundle, each of whose fibers has a 
!   Chart_SLL, over a base Atlas_SC with a Chart_SLD.

module Bundle_SLL_ASC_CSLD__Form

  !-- Bundle_SingleLevelLocal_AtlasSingleChart_ChartSingleLevelDistributed_Form

  use Basics
  use Atlases
  use BundleHeader_Form
  use FibersWritten_CSL__Form

  implicit none
  private

  type, public, extends ( BundleHeaderForm ) :: Bundle_SLL_ASC_CSLD_Form
    integer ( KDI ) :: &
      nBaseValues    = 0, &
      nFibers        = 0, &
      nSections      = 0, &
      nFibersWrite   = 0, &
      sFibersWrite   = 0, &
      nSectionsWrite = 0, &
      sSectionsWrite = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaBaseCell, &
      iaBaseCellLabel
    class ( Chart_SLL_Form ), pointer :: &
      Fiber_CSLL => null ( )
    class ( Chart_SLD_Form ), pointer :: &
      Base_CSLD => null ( )
    class ( Atlas_SC_Form ), pointer :: &
      Base_ASC => null ( )
    type ( GeometryFlat_ASC_Form ), allocatable :: &
      Geometry
    type ( Atlas_SC_1D_Form ), allocatable :: &
      Fiber
    type ( FibersWritten_CSL_Form ), allocatable :: &
      FibersWritten
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      CreateChart
    procedure, public, pass :: &
      GeometryFiber
    procedure, public, pass :: &
      AtlasFiber
    procedure, public, pass :: &
      OpenStream
    procedure, public, pass :: &
      MarkFibersWritten
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      CloseStreams
    procedure, private, pass :: &
      ShowHeader
    final :: &
      Finalize
  end type Bundle_SLL_ASC_CSLD_Form

contains


  subroutine Initialize ( B, Base, Name )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      B
    class ( AtlasHeaderForm ), intent ( inout ), target :: &
      Base
    character ( * ), intent ( in )  :: &
      Name
      
    if ( B % Type == '' ) &
      B % Type = 'a Bundle_SLL_ASC_CSLD' 

    call B % BundleHeaderForm % Initialize ( Base, Name )

    select type ( AB => Base )
    class is ( Atlas_SC_Form )

      B % Base_ASC => AB
      
      select type ( CB => AB % Chart )
      class is ( Chart_SLD_Form )

        B % Base_CSLD   => CB
        B % nBaseValues =  CB % nProperCells + CB % nGhostCells
        B % nFibers     =  CB % nProperCells

        B % nFibersWrite = 5
        call PROGRAM_HEADER % GetParameter ( B % nFibersWrite, 'nFibersWrite' )

        B % sFibersWrite = B % nFibers / B % nFibersWrite
        if ( mod ( B % nFibers, B % nFibersWrite ) > 0 ) &
          B % sFibersWrite = B % sFibersWrite + 1

        allocate ( B % FibersWritten )
        associate ( FW => B % FibersWritten )
        call FW % Initialize ( CB, B % nBaseValues )
        call CB % AddField ( FW )
        end associate !-- FW
        
      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Bundle_SLL_ASC_CSLD__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- CB
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Bundle_SLL_ASC_CSLD__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- AB

  end subroutine Initialize


  subroutine CreateChart &
               ( B, SpacingOption, CoordinateSystemOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nCellsOption, nGhostLayersOption, nDimensionsOption )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ), target :: &
      B
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption

    integer ( KDI ) :: &
      iF, &    !-- iFiber
      iCN, &   !-- iCellNumber
      iCNL, &  !-- iCellNumberLabel
      iCell, jCell, kCell
    integer ( KDI ), dimension ( 3 ) :: &
      iaCell
    class ( GeometryFlatForm ), pointer :: &
      GF

    associate ( FM => B % FiberMaster )
    call FM % CreateChart &
           ( SpacingOption, CoordinateSystemOption, &
             CoordinateUnitOption, MinCoordinateOption, &
             MaxCoordinateOption, RatioOption, ScaleOption, &
             nCellsOption, nGhostLayersOption, nDimensionsOption )

    allocate ( B % Geometry )
    associate ( GAF => B % Geometry )
    call GAF % Initialize ( FM )

    select type ( CF => B % FiberMaster % Chart )
    class is ( Chart_SLL_Form )
      B % Fiber_CSLL => CF
      B % nSections  =  CF % nProperCells
    end select !-- CF

    B % nSectionsWrite = 5
    call PROGRAM_HEADER % GetParameter ( B % nSectionsWrite, 'nSectionsWrite' )

    B % sSectionsWrite = B % nSections / B % nSectionsWrite
    if ( mod ( B % nSections, B % nSectionsWrite ) > 0 ) &
      B % sSectionsWrite = B % sSectionsWrite + 1

    call B % Show ( )

    call FM % SetGeometry ( GAF )

    allocate ( B % Fiber )
    associate ( F => B % Fiber )
    call F % Initialize ( B % nFibers )

    do iF = 1, B % nFibers
      allocate ( F % Atlas ( iF ) % Element )
      associate ( FA => F % Atlas ( iF ) % Element )
      call FA % Initialize ( FM )
      call FA % SetGeometry ( GAF )
      end associate !-- FA
    end do !-- iF

    end associate !-- F
    end associate !-- GAF
    end associate !-- FM

    !-- Base cell index

    allocate ( B % iaBaseCell ( B % nFibers ) )
    allocate ( B % iaBaseCellLabel ( B % nFibers ) )
    associate ( CB => B % Base_CSLD )

    iCN  = 0
    iCNL = 0
    iF   = 0
    do kCell = CB % iaFirst ( 3 ), CB % iaLast ( 3 )
      do jCell = CB % iaFirst ( 2 ), CB % iaLast ( 2 )
        do iCell = CB % iaFirst ( 1 ), CB % iaLast ( 1 )

          iaCell = [ iCell, jCell, kCell ]

          iCN = iCN + 1

          !-- Only proper and ghost (not exterior) cells are written
          if ( .not. any ( CB % iaBrick == 1 .and. iaCell < 1 ) &
               .and. .not. any ( CB % iaBrick == CB % nBricks &
                                 .and. iaCell > CB % nCellsBrick ) ) &
          then
            iCNL = iCNL + 1
          end if

          !-- Only proper cells are fibers
          if ( any ( iaCell < 1 ) ) cycle
          if ( any ( iaCell > CB % nCellsBrick ) ) cycle 
              
          iF = iF + 1
          B % iaBaseCell      ( iF ) = iCN
          B % iaBaseCellLabel ( iF ) = iCNL
              
        end do !-- iCell
      end do !-- jCell
    end do !-- kCell
    end associate !-- CB

    nullify ( GF )

  end subroutine CreateChart


  function GeometryFiber ( B ) result ( G )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( AF => B % Geometry % Atlas )
    class is ( Atlas_SC_Form )
      G => AF % Geometry ( )
    end select !-- AF

  end function GeometryFiber


  function AtlasFiber ( B, iFiber ) result ( AF )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iFiber
    class ( Atlas_SC_Form ), pointer :: &
      AF

    AF => B % Fiber % Atlas ( iFiber ) % Element 

  end function AtlasFiber


  subroutine OpenStream ( B, GIS_Base, iStream )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      B
    type ( GridImageStreamForm ), intent ( in ) :: &
      GIS_Base
    integer ( KDI ), intent ( in ) :: &
      iStream

    integer ( KDI ) :: &
      iF  !-- iFiber
    class ( Atlas_SC_Form ), pointer :: &
      AF
    character ( LDN + 1 ) :: &
      BlockNumber
    character ( LDF ) :: &
      BundleDirectory, &
      BundleName

    if ( .not. allocated ( B % GridImageStream ) ) &
      allocate ( B % GridImageStream ( ATLAS % MAX_STREAMS ) )

    write ( BlockNumber, fmt = '(a1,i7.7)' ) '_', &
      GIS_Base % Communicator % Rank
!    BundleDirectory &
!      = trim ( GIS_Base % WorkingDirectory ) // trim ( GIS_Base % Name ) &
!        // '_Bundle_MeshBlock' // trim ( BlockNumber ) // '/' 
    BundleDirectory &
      = trim ( GIS_Base % WorkingDirectory ) // 'Bundle_MeshBlock' &
        // trim ( BlockNumber ) // '/' 
!    BundleName &
!      = trim ( GIS_Base % Name ) // '_Bundle_MeshBlock' // trim ( BlockNumber )
    BundleName &
      = trim ( GIS_Base % Name ) // '_Bundle'
!    call Show ( BundleDirectory, 'BundleDirectory' )
!    call Show ( BundleName, 'BundleName' )

    associate ( GIS_Bundle => B % GridImageStream ( iStream ) )
    call GIS_Bundle % Initialize &
           ( BundleName, WorkingDirectoryOption = BundleDirectory )
    end associate !-- GIS_Bundle

    nullify ( AF )

  end subroutine OpenStream


  subroutine MarkFibersWritten ( B, AllFibersOption )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      B
    logical ( KDL ), intent ( in ), optional :: &
      AllFibersOption

    integer ( KDI ) :: &
      iF, &  !-- iFiber 
      sF     !-- sFibers
    logical ( KDL ) :: &
      AllFibers

    AllFibers = .false.
    if ( present ( AllFibersOption ) ) &
      AllFibers = AllFibersOption

    associate &
      ( FW   => B % FibersWritten % Field % Value ( :, 1 ), &
        iaBC => B % iaBaseCell, &
        nF   => B % nFibers )
    call Clear ( FW )

    sF = B % sFibersWrite
    if ( AllFibers ) &
      sF = 1

    do iF = 1, nF, sF
      FW ( iaBC ( iF ) ) = 1.0_KDR
    end do !-- iF

    end associate !-- FW, etc.
        
  end subroutine MarkFibersWritten


  subroutine Write ( B, iStream, AllFibersOption, TimeOption, &
                     CycleNumberOption )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iStream
    logical ( KDL ), intent ( in ), optional :: &
      AllFibersOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption

    integer ( KDI ) :: &
      iF, &   !-- iFiber
      sF  !-- sFibers
    logical ( KDL ) :: &
      AllFibers
    character ( LDN + 1 ) :: &
      CellNumber
    character ( LDF ) :: &
      Directory
    class ( Atlas_SC_Form ), pointer :: &
      AF

    AllFibers = .false.
    if ( present ( AllFibersOption ) ) &
      AllFibers = AllFibersOption

    associate ( GIS_Bundle => B % GridImageStream ( iStream ) )
    call GIS_Bundle % Open ( GIS_Bundle % ACCESS_CREATE )

    associate &
      ( iaBCL => B % iaBaseCellLabel, &
        nF    => B % nFibers )

    sF = B % sFibersWrite
    if ( AllFibers ) &
      sF = 1

    do iF = 1, nF, sF              

      write ( CellNumber, fmt = '(a1,i7.7)' ) '_', iaBCL ( iF )
      Directory &
        = 'Fiber' // CellNumber // '/'
!      call Show ( Directory, 'Directory', CONSOLE % INFO_3 )
!      call Show ( iF, 'iFiber', CONSOLE % INFO_3 )

      AF => B % AtlasFiber ( iF ) 
      call AF % OpenStream ( GIS_Bundle, '', iStream = iStream ) 
      call AF % Write &
             ( iStream = iStream, DirectoryOption = Directory, &
               TimeOption = TimeOption, CycleNumberOption = CycleNumberOption )
      call AF % CloseStreams ( ) 

    end do !-- iF
    end associate !-- iaBC, etc.

    call GIS_Bundle % Close ( )
    end associate !-- GIS_Bundle

    nullify ( AF )

  end subroutine Write


  subroutine CloseStreams ( B )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      B

    integer ( KDI ) :: &
      iF  !-- iFiber
    class ( Atlas_SC_Form ), pointer :: &
      AF

    if ( allocated ( B % GridImageStream ) ) &
      deallocate ( B % GridImageStream )

  end subroutine CloseStreams


  subroutine ShowHeader ( B )

    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      B

    call B % BundleHeaderForm % Show ( )

    call Show ( B % nBaseValues, 'nBaseValues', B % IGNORABILITY )
    call Show ( B % nFibers, 'nFibers', B % IGNORABILITY )
    call Show ( B % nSections, 'nSections', B % IGNORABILITY )
    call Show ( B % nFibersWrite, 'nFibersWrite', B % IGNORABILITY )
    call Show ( B % sFibersWrite, 'sFibersWrite', B % IGNORABILITY )
    call Show ( B % nSectionsWrite, 'nSectionsWrite', B % IGNORABILITY )
    call Show ( B % sSectionsWrite, 'sSectionsWrite', B % IGNORABILITY )

  end subroutine ShowHeader


  subroutine Finalize ( B )

    type ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      B

    if ( allocated ( B % FibersWritten ) ) &
      deallocate ( B % FibersWritten )
    if ( allocated ( B % Fiber ) ) &
      deallocate ( B % Fiber )
    if ( allocated ( B % Geometry ) ) &
      deallocate ( B % Geometry )
    if ( allocated ( B % iaBaseCellLabel ) ) &
      deallocate ( B % iaBaseCellLabel )
    if ( allocated ( B % iaBaseCell ) ) &
      deallocate ( B % iaBaseCell )

    nullify ( B % Base_ASC )
    nullify ( B % Base_CSLD )
    nullify ( B % Fiber_CSLL )

  end subroutine Finalize


end module Bundle_SLL_ASC_CSLD__Form
