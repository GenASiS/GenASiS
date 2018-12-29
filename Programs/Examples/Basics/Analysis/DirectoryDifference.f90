program DirectoryDifference

  !-- Calculate the L1 error between two Directory

  use Basics
  use L1_Norm_Function
  
  implicit none
  
  integer ( KDI ) :: &
    iV, &
    iS, &
    iFile, &
    FilenumberStart, &
    FilenumberEnd
  integer ( KDI ), dimension ( 3 ) :: &
    iaFirst, &
    iaLast, &
    nProperCells
  real ( KDR ), dimension ( : ), allocatable :: &
    L1_Error
  real ( KDR ), dimension ( :, :, : ), pointer :: &
    V_1, V_2, V_O
  character ( LDF ) :: &
    Filename, &
    Directory_1, &
    Directory_2, &
    OutputDirectory
  type ( StorageForm ), dimension ( : ), allocatable :: &
    Storage
  type ( GridImageStreamForm ), allocatable, target :: &
    GIS_1, &
    GIS_2, &
    GIS_O       !-- GIS for writing output
  class ( GridImageSiloTemplate ), allocatable :: &
    GI_1, &
    GI_2, &
    GI_O        !-- GI for writing output
    
  allocate ( PROGRAM_HEADER )
  
  associate ( PH => PROGRAM_HEADER )
  call PH % Initialize &
         ( 'DirectoryDifference', AppendDimensionalityOption = .false. )
  
  Filename = 'FilenameNotSet'
  call PH % GetParameter ( Filename, 'Filename' )
  
  FilenumberStart = 0
  call PH % GetParameter ( FilenumberStart, 'FilenumberStart' )
  
  FilenumberEnd = huge ( 1_KDI )
  call PH % GetParameter ( FilenumberEnd, 'FilenumberEnd' )
  
  Directory_1 = '../Output/'
  call PH % GetParameter ( Directory_1, 'Directory_1' )
  
  Directory_2 = '../Output/'
  call PH % GetParameter ( Directory_2, 'Directory_2' )
  
  OutputDirectory = '../Output/'
  call PH % GetParameter ( OutputDirectory, 'OutputDirectory' )
  
  allocate ( GIS_1 )
  call GIS_1 % Initialize &
        ( Filename, PH % Communicator, &
          WorkingDirectoryOption = Directory_1 )
          
  allocate ( GIS_2 )
  call GIS_2 % Initialize &
        ( Filename, PH % Communicator, &
          WorkingDirectoryOption = Directory_2 )
  
  allocate ( GIS_O )
  call GIS_O % Initialize &
        ( trim ( Filename ) // '_Difference', &
          PH % Communicator, &
          WorkingDirectoryOption = OutputDirectory )
  
  do iFile = FilenumberStart, FilenumberEnd
    
    call Show ( iFile, 'Opening Filenumber' )
  
    call GIS_1 % Open ( GIS_1 % ACCESS_READ, NumberOption = iFile )
    call GIS_1 % ListContents ( ContentTypeOption = 'Directory' )
    
    if ( trim ( GIS_1 % ContentList ( 1 ) ) == 'Curves' ) then
      allocate ( CurveImageForm :: GI_1 )
      allocate ( CurveImageForm :: GI_2 )
      allocate ( CurveImageForm :: GI_O )
    else if ( trim ( GIS_1 % ContentList ( 1 ) ) == 'Grid' ) then
      allocate ( StructuredGridImageForm :: GI_1 )
      allocate ( StructuredGridImageForm :: GI_2 )
      allocate ( GI_O, source = GI_1 )
    end if
    
    call Show ( GIS_1 % ContentList, 'ContentList 1' )
      
    call GI_1 % Initialize ( GIS_1 )
    call GI_1 % SetReadAttributes &
           ( Directory = trim ( GIS_1 % ContentList ( 1 ) ), oValue = 0 )
    call GI_1 % Read ( )
    
    call GIS_1 % Close ( )
   
    call Show ( GIS_1 % ContentList, 'ContentList 2' )
    
    call GIS_2 % Open ( GIS_2 % ACCESS_READ, NumberOption = iFile )
    call GIS_2 % ListContents ( ContentTypeOption = 'Directory' )
    
    call GI_2 % Initialize ( GIS_2 )
    call GI_2 % SetReadAttributes &
           ( Directory = trim ( GIS_2 % ContentList ( 1 ) ), oValue = 0 )
    call GI_2 % Read ( )
    
    call GIS_2 % Close ( )
    
    
    call Show ( '<<< GIS_2 close' )
    call GIS_O % Open &
           ( GIS_O % ACCESS_CREATE, SeriesOption = .true., &
             NumberOption = iFile )
    
    !call Show ( '<<< Stream opened' )
    !call GI_O % Initialize ( GIS_O )
    GI_O % Stream => GIS_O
             
    call Show ( '<<< GI Initialized' )
    allocate ( Storage ( GI_1 % nStorages ) )
    call Show ( '<<<< Storages allocated' )
    
    
    select type ( GI_1 )
    type is ( CurveImageForm )
      iaFirst = 1 
      iaLast  = [ GI_1 % nTotalCells, 1, 1 ]
      nProperCells = [ GI_1 % nTotalCells, 1, 1 ]
    type is ( StructuredGridImageForm )
      iaFirst = 1 - GI_1 % oProperCellNodeLow
      iaLast  = max ( GI_1 % nCells - GI_1 % oProperCellNodeLow, 1 )
      nProperCells &
        = max ( GI_1 % nCells - GI_1 % oProperCellNodeLow &
                  - GI_1 % oProperCellNodeHigh, 1 )
    end select
    
    !select type ( GI_O )
    !type is ( StructuredGridImageForm )
    !  select type ( GI_1 )
    !  type is ( StructuredGridImageForm )
    !    call Show ( GIS_1 % ContentList, 'ContentList 3' )
    !    call GI_O % SetGrid &
    !           ( Directory = 'Grid', &
    !             NodeCoordinate &
    !               = reshape ( [ GI_1 % NodeCoordinate_1, &
    !                             GI_1 % NodeCoordinate_2, &
    !                             GI_1 % NodeCoordinate_3 ], &
    !                           [ size ( GI_1 % NodeCoordinate_1 ), 3 ], &
    !                           pad = [ 0.0_KDR ] ), &
    !             nDimensions = GI_1 % nDimensions, &
    !             nProperCells = GI_1 % nTotalCells - GI_1 % nGhostCells, &
    !             nGhostCells = GI_1 % nGhostCells, &
    !             oValue = 0, nCells = GI_1 % nCells )
    !  end select 
    !end select
      
    !-- Calculate and write the differences  
    do iS = 1, GI_1 % nStorages
    
      associate &
        ( S_1 => GI_1 % Storage ( iS ), &
          S_2 => GI_2 % Storage ( iS ), &
          S_O => Storage ( iS ) )
      
      call Show ( '<<<< initializin storage' )
      call S_O % Initialize &
             ( shape ( S_1 % Value ), &
               VectorIndicesOption = S_1 % VectorIndices, &
               UnitOption = S_1 % Unit, VectorOption = S_1 % Vector, &
               VariableOption = S_1 % Variable, NameOption = S_1 % Name, &
               ClearOption = .true. )
      call Show ( '<<<< storage initialized' )
      !call GI_O % AddStorage ( Storage ( iS ) )
      !call Show ( '<<<< storage added' )
      
      allocate ( L1_Error ( S_1 % nVariables ) )
      
      do iV = 1, S_1 % nVariables
        !-- Use pointer remapping to exclude computing error on ghost cells
        V_1 &
          ( iaFirst ( 1 ) : iaLast ( 1 ), &
            iaFirst ( 2 ) : iaLast ( 2 ), &  
            iaFirst ( 3 ) : iaLast ( 3 ) ) => S_1 % Value ( :, iV )
        V_2 &
          ( iaFirst ( 1 ) : iaLast ( 1 ), &
            iaFirst ( 2 ) : iaLast ( 2 ), &  
            iaFirst ( 3 ) : iaLast ( 3 ) ) => S_2 % Value ( :, iV )
        V_O &
          ( iaFirst ( 1 ) : iaLast ( 1 ), &
            iaFirst ( 2 ) : iaLast ( 2 ), &  
            iaFirst ( 3 ) : iaLast ( 3 ) ) => S_O % Value ( :, iV )
        
        L1_Error ( iV ) &
          = L1_Norm ( V_1 ( 1 : nProperCells ( 1 ), &
                            1 : nProperCells ( 2 ), &
                            1 : nProperCells ( 3 ) ), &
                      V_2 ( 1 : nProperCells ( 1 ), &
                            1 : nProperCells ( 2 ), &
                            1 : nProperCells ( 3 ) ) )
        
        !-- Avoid computing fractional differences for vanishing variable
        if ( V_1 ( 1, 1, 1 ) == 0.0_KDR ) cycle
        
        associate &
          ( V_1_P => V_1 ( 1 : nProperCells ( 1 ), &
                           1 : nProperCells ( 2 ), &
                           1 : nProperCells ( 3 ) ), &
            V_2_P => V_1 ( 1 : nProperCells ( 1 ), &
                           1 : nProperCells ( 2 ), &
                           1 : nProperCells ( 3 ) ), &
            V_O_P => V_O ( 1 : nProperCells ( 1 ), &
                           1 : nProperCells ( 2 ), &
                           1 : nProperCells ( 3 ) ) ) 
                           
          V_O_P = ( V_2_P - V_1_P ) / V_1_P
        end associate 
        
      end do  !-- nVariables
                    
      call Show ( 'L1 Norm Error for ' // trim ( S_1 % Name ) )
      call Show ( L1_Error, S_1 % Variable )
      
      deallocate ( L1_Error )
      
      end associate
      
    end do    !-- nStorages
    
    call GI_O % Write ( )
    !GI_1 % Stream => GIS_O
    !call GI_1 % Write ( )
    
    call Show ('<<<< closing stream' )
    call GIS_O % Close ( )
    
    deallocate ( Storage )
    
    deallocate ( GI_O )
    deallocate ( GI_2 )
    deallocate ( GI_1 )
  
  end do
  
  deallocate ( GIS_O )
  deallocate ( GIS_1 ) 
  deallocate ( GIS_2 )
    
  
    
!  allocate ( GIS )
!  
!  call GIS % Initialize &
!         ( trim ( Filename ) // '_Difference', PH % Communicator, &
!           WorkingDirectoryOption = OutputDirectory_1 )
!  
!  call GIS % Open ( GIS % ACCESS_CREATE, SeriesOption = .false. )
!  
!  GI_1 % Stream => GIS
!  call GI_1 % Write ( )
!  
!  call GIS % Close ( )
!  
!  deallocate ( GI_2 )
!  deallocate ( GI_1 )
!  
!  deallocate ( GIS )
  
  deallocate ( PROGRAM_HEADER )
  
  end associate !-- PH

end program DirectoryDifference
