program StructuredGridImage_Form_Test

  use VariableManagement
  use Display
  use MessagePassing
  use GridImageStream_Form
  use StructuredGridImage_Form
  
  implicit none
  
  integer ( KDI ) :: &
    oC_1   !-- oCoordinate_1
  real ( KDR ), dimension ( 20, 3 ) :: &
    NodeCoordinate
  character ( LDF ) :: &
    Name = 'StructuredGridImage_Form_Test'
  character ( LDL ), dimension ( 3 ) :: &
    VariableName
  character ( LDL ), dimension ( 1 ) :: &
    VectorName
  type ( Integer_1D_Form ), dimension ( 1 ) :: &
    VectorIndices
  type ( VariableGroupForm ) :: &
    VG
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( GridImageStreamForm ) :: &
    GIS
  type ( StructuredGridImageForm ) :: & 
    SGI_Rect, &
    SGI_Read
  
  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetVerbosity ( 'INFO_7' )
  
  NodeCoordinate ( :, 1 : 2  ) &
    = reshape &
        ( source = [ 0., 2.25, 1., 2.25, 2.5, 2.25, 5., 2.25, &
                     0., 2.55, 1., 2.55, 2.5, 2.55, 5., 2.55, &      
                     0., 0.,   1., 0.,   2.5, 0.,   5., 0., &
                     0., 5.,   1., 5.,   2.5, 5.,   5., 5., &
                     0., 2.,   1., 2.,   2.5, 2.,   5., 2. ], & 
          shape = [ 20, 2 ], order = [ 2, 1 ] )
  NodeCoordinate ( :, 3 ) = 0.0_KDR
  
  oC_1 = C % Rank * maxval ( NodeCoordinate ( :, 1 ) )
  NodeCoordinate ( :, 1 ) = NodeCoordinate ( :, 1 ) + oC_1
  
  call VectorIndices ( 1 ) % Initialize ( [ 1, 2, 3 ] )
  
  VariableName ( 1 ) = 'VariableName_1'
  VariableName ( 2 ) = 'VariableName_2'
  VariableName ( 3 ) = 'VariableName_3'
  VectorName ( 1 ) = 'VectorName'
  call VG % InitializeAllocate &
         ( [ 12, 3 ], VariableOption = VariableName, &
           NameOption = 'VariableGroup', VectorOption = VectorName, &
           VectorIndicesOption = VectorIndices )
  VG % Value ( :, 1 ) &
    = [ 1., 2., 3., 4., 5., 6., 2., 4., 6., 8., 10., 12. ]
  VG % Value ( :, 2 ) = 2.0_KDR * VG % Value ( :, 1 )
  VG % Value ( :, 3 ) = 3.0_KDR * VG % Value ( :, 1 )
    
  call GIS % Initialize ( Name, CommunicatorOption = C )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  
  call SGI_Rect % Initialize ( GIS )

  call SGI_Rect % AddVariableGroup ( VG )
  
  call SGI_Rect % SetGrid  &
         ( Directory = 'Rectilinear', &
           NodeCoordinate = NodeCoordinate, nDimensions = 2, &
           nProperCells = 12, nGhostCells = 0, oValue = 0, &
           nCells = [ 3, 4, 0 ] )

  call Show ( SGI_Rect % nDimensions, 'nDimensions' )  
  call Show ( SGI_Rect % nNodes, 'nNodes' )
  call Show ( SGI_Rect % NodeCoordinate_1, 'NodeCoordinate_1' )
  call Show ( SGI_Rect % NodeCoordinate_2, 'NodeCoordinate_2' )
  call Show ( SGI_Rect % NodeCoordinate_3, 'NodeCoordinate_3' )
  call Show &
         ( SGI_Rect % VariableGroup ( 1 ) % Value ( :, 1 ), &
           SGI_Rect % VariableGroup ( 1 ) % Variable ( 1 ) )
  call Show &
         ( SGI_Rect % VariableGroup ( 1 ) % Value ( :, 2 ), &
           SGI_Rect % VariableGroup ( 1 ) % Variable ( 2 ) )
  call Show &
         ( SGI_Rect % VariableGroup ( 1 ) % Value ( :, 3 ), &
           SGI_Rect % VariableGroup ( 1 ) % Variable ( 3 ) )

  call SGI_Rect % Write ( )
  
  call GIS % Close ( ) 
  
  call Clear ( VG % Value )

  call GIS % Open ( GIS % ACCESS_READ, NumberOption = 0 )
  
  call SGI_Read % Initialize ( GIS )
  
  call SGI_Read % SetReadAttributes ( Directory = 'Rectilinear', oValue = 0 )
  
  call SGI_Read % Read ( )
  
  call Show ( SGI_Read % nDimensions, 'nDimensions' )  
  call Show ( SGI_Read % nNodes, 'nNodes' )
  call Show ( SGI_Read % NodeCoordinate_1, 'NodeCoordinate_1' )
  call Show ( SGI_Read % NodeCoordinate_2, 'NodeCoordinate_2' )
  call Show ( SGI_Read % NodeCoordinate_3, 'NodeCoordinate_3' )
  call Show &
         ( SGI_Read % VariableGroup ( 1 ) % Value ( :, 1 ), &
           SGI_Read % VariableGroup ( 1 ) % Variable ( 1 ) )
  call Show &
         ( SGI_Read % VariableGroup ( 1 ) % Value ( :, 2 ), &
           SGI_Read % VariableGroup ( 1 ) % Variable ( 2 ) )
  call Show &
         ( SGI_Read % VariableGroup ( 1 ) % Value ( :, 3 ), &
           SGI_Read % VariableGroup ( 1 ) % Variable ( 3 ) )
  
  call GIS % Close ( )

  deallocate ( C )

end program StructuredGridImage_Form_Test
