program FileDifference

  use Basics
  use L1_Norm_Function
  
  implicit none
  
  integer ( KDI ) :: &
    iV, &
    iS, &
    Filenumber_1, &
    Filenumber_2
  integer ( KDI ), dimension ( 3 ) :: &
    iaFirst, &
    iaLast, &
    nProperCells
  real ( KDR ), dimension ( : ), allocatable :: &
    L1_Error
  real ( KDR ), dimension ( :, :, : ), pointer :: &
    V_1, V_2
  character ( LDF ) :: &
    Filename, &
    OutputDirectory
  type ( GridImageStreamForm ), allocatable, target :: &
    GIS
  class ( GridImageSiloTemplate ), allocatable :: &
    GI_1, &
    GI_2
    
  allocate ( PROGRAM_HEADER )
  
  associate ( PH => PROGRAM_HEADER )
  call PH % Initialize &
         ( 'ProgramOutputDifference', AppendDimensionalityOption = .false. )
  
  Filename = 'FilenameNotSet'
  call PH % GetParameter ( Filename, 'Filename' )
  
  Filenumber_1 = 0
  call PH % GetParameter ( Filenumber_1, 'Filenumber_1' )
  
  Filenumber_2 = 0
  call PH % GetParameter ( Filenumber_2, 'Filenumber_2' )
  
  OutputDirectory = '../Output/'
  call PH % GetParameter ( OutputDirectory, 'OutputDirectory' )
  
  allocate ( GIS )
  call GIS % Initialize &
        ( Filename, PH % Communicator, &
          WorkingDirectoryOption = OutputDirectory )
  
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = Filenumber_1 )
  call GIS % ListContents ( ContentTypeOption = 'Directory' )
  
  if ( trim ( GIS % ContentList ( 1 ) ) == 'Curves' ) then
    allocate ( CurveImageForm :: GI_1 )
    allocate ( CurveImageForm :: GI_2 )
  else if ( trim ( GIS % ContentList ( 1 ) ) == 'Grid' ) then
    allocate ( StructuredGridImageForm :: GI_1 )
    allocate ( StructuredGridImageForm :: GI_2 )
  end if
  
  call GI_1 % Initialize ( GIS )
  call GI_1 % SetDirectory ( trim ( GIS % ContentList ( 1 ) ) )
  call GI_1 % Read ( )
  
  call GIS % Close ( )
 
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = Filenumber_2 )
  call GIS % ListContents ( ContentTypeOption = 'Directory' )
  
  call GI_2 % Initialize ( GIS )
  call GI_2 % SetDirectory ( trim ( GIS % ContentList ( 1 ) ) )
  call GI_2 % Read ( )
  
  call GIS % Close ( )
  
  deallocate ( GIS )
  
  !-- Calculate and write the differences
  
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
    
  do iS = 1, GI_1 % nStorages
    
    associate &
      ( S_1 => GI_1 % Storage ( iS ), &
        S_2 => GI_2 % Storage ( iS ) )
    
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
      
      L1_Error ( iV ) &
        = L1_Norm ( V_1 ( 1 : nProperCells ( 1 ), &
                          1 : nProperCells ( 2 ), &
                          1 : nProperCells ( 3 ) ), &
                    V_2 ( 1 : nProperCells ( 1 ), &
                          1 : nProperCells ( 2 ), &
                          1 : nProperCells ( 3 ) ) )
      
      !-- Avoid computing fractional differences for vanishing variable
      if ( V_1 ( 1, 1, 1 ) == 0.0_KDR ) cycle
      
      V_1 = ( V_2 - V_1 ) / V_1
      
    end do
                  
    call Show ( 'L1 Norm Error for ' // trim ( S_1 % Name ) )
    call Show ( L1_Error, S_1 % Variable )
    
    deallocate ( L1_Error )
    
    end associate
  end do
  
  allocate ( GIS )
  
  call GIS % Initialize &
         ( trim ( Filename ) // '_Difference', PH % Communicator, &
           WorkingDirectoryOption = OutputDirectory )
  
  call GIS % Open ( GIS % ACCESS_CREATE, SeriesOption = .false. )
  
  GI_1 % Stream => GIS
  call GI_1 % Write ( )
  
  call GIS % Close ( )
  
  deallocate ( GI_2 )
  deallocate ( GI_1 )
  
  deallocate ( GIS )
  
  deallocate ( PROGRAM_HEADER )
  
  end associate !-- PH

end program FileDifference
