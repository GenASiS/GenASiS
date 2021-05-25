program GatherMeshBlocksMinMaxAve

  !-- Gather the Minimum, Maximum, and Average from all Meshblocks and write
  !   it to an ASCII file as table data

  use Basics
  
  implicit none
  
  integer ( KDI ) :: &
    iV, &
    iS, &
    iR, &
    Status, &
    Filenumber, &
    WriteUnit
  integer ( KDI ), dimension ( 3 ) :: &
    iaFirst, &
    iaLast, &
    nProperCells
  real ( KDR ), dimension ( :, : ), allocatable :: &
    G_MMA  !-- Gathered MinMaxAve per Storage
  real ( KDR ), dimension ( :, :, : ), pointer :: &
    V_1
  character ( LDF ) :: &
    Filename, &
    OutputDirectory
  character ( LDL ) :: &
    WriteFormat
  type ( CollectiveOperation_R_Form ), allocatable :: &
    CO
  type ( GridImageStreamForm ), allocatable, target :: &
    GIS
  class ( GridImageSiloTemplate ), allocatable :: &
    GI
    
  allocate ( PROGRAM_HEADER )
  
  associate ( PH => PROGRAM_HEADER )
  
  call PH % Initialize &
         ( 'GatherMeshBlocksMinMaxAve', AppendDimensionalityOption = .false. )
  
  associate ( C => PROGRAM_HEADER % Communicator )
         
  Filename = 'FilenameNotSet'
  call PH % GetParameter ( Filename, 'Filename' )
  
  Filenumber = 0
  call PH % GetParameter ( Filenumber, 'Filenumber' )
  
  OutputDirectory = '../Output/'
  call PH % GetParameter ( OutputDirectory, 'OutputDirectory' )
  
  allocate ( GIS )
  call GIS % Initialize &
        ( Filename, PH % Communicator, &
          WorkingDirectoryOption = OutputDirectory )
  
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = Filenumber )
  call GIS % ListContents ( ContentTypeOption = 'Directory' )
  
  if ( trim ( GIS % ContentList ( 1 ) ) == 'Curves' ) then
    allocate ( CurveImageForm :: GI )
  else if ( trim ( GIS % ContentList ( 1 ) ) == 'Grid' ) then
    allocate ( StructuredGridImageForm :: GI )
  end if
  
  call GI % Initialize ( GIS )
  call GI % SetDirectory ( trim ( GIS % ContentList ( 1 ) ) )
  call GI % Read ( )
  
  call GIS % Close ( )
 
  deallocate ( GIS )
  
  select type ( GI )
  type is ( CurveImageForm )
    iaFirst = 1 
    iaLast  = [ GI % nTotalCells, 1, 1 ]
    nProperCells = [ GI % nTotalCells, 1, 1 ]
  type is ( StructuredGridImageForm )
    iaFirst = 1 - GI % oProperCellNodeLow
    iaLast  = max ( GI % nCells - GI % oProperCellNodeLow, 1 )
    nProperCells &
      = max ( GI % nCells - GI % oProperCellNodeLow &
                - GI % oProperCellNodeHigh, 1 )
  end select
    
  do iS = 1, GI % nStorages
  
    associate ( S => GI % Storage ( iS ) )
  
    allocate ( CO )
    call CO % Initialize &
           ( C, &
             nOutgoing = [ 3 ], nIncoming = [ C % Size * 3 ], &
             RootOption = CONSOLE % DisplayRank )
    
    allocate ( G_MMA ( C % Size * 3, S % nVariables ) )
    
    do iV = 1, S % nVariables
      
      !-- Use pointer remapping to exclude computing error on ghost cells
      V_1 &
        ( iaFirst ( 1 ) : iaLast ( 1 ), &
          iaFirst ( 2 ) : iaLast ( 2 ), &  
          iaFirst ( 3 ) : iaLast ( 3 ) ) => S % Value ( :, iV )
          
      CO % Outgoing % Value ( 1 ) &
        = minval ( V_1 ( 1 : nProperCells ( 1 ), &
                         1 : nProperCells ( 2 ), &
                         1 : nProperCells ( 3 ) ) )
      
      CO % Outgoing % Value ( 2 ) &
        = maxval ( V_1 ( 1 : nProperCells ( 1 ), &
                         1 : nProperCells ( 2 ), &
                         1 : nProperCells ( 3 ) ) )
      
      CO % Outgoing % Value ( 3 ) &
        = sum ( V_1 ( 1 : nProperCells ( 1 ), &
                      1 : nProperCells ( 2 ), &
                      1 : nProperCells ( 3 ) ) ) &
            / ( sum ( nProperCells ) )
      
      call CO % Gather ( )
      
      if ( C % Rank == CONSOLE % DisplayRank ) then
        G_MMA ( :, iV ) = CO % Incoming % Value
        call Show ( CO % Incoming % Value, S % Variable ( iV ), &
                    IgnorabilityOption = CONSOLE % INFO_3 )
      end if
      
    end do
    
    if ( C % Rank == CONSOLE % DisplayRank ) then
      open ( newunit = WriteUnit, &
             file = trim ( OutputDirectory ) // trim ( S % Name ) &
                    // '_MinMaxAve.txt', &
             action = 'write' )
      
      WriteFormat = ''
      write ( WriteFormat, fmt = '(a1,i4,a9)' ) &
              '(', S % nVariables, 'E20.10E3)'
      
      do iR = 1, size ( G_MMA, dim = 1 ) 
        write ( WriteUnit, fmt = WriteFormat, ioStat = Status ) G_MMA ( iR, : )
      end do
    end if
    
    
    deallocate ( G_MMA )
    deallocate ( CO )
    
    end associate
  
  end do
  
  end associate !-- C
  
  deallocate ( PROGRAM_HEADER )
  
  end associate !-- PH

end program GatherMeshBlocksMinMaxAve
