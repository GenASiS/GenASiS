program UnstructuredGridImage_Form_Test

  use VariableManagement
  use Display
  use MessagePassing
  use GridImageStream_Form
  use UnstructuredGridImage_Form
  
  implicit none
  
  integer ( KDI ) :: &
    iN, &      !-- iNode
    oC_1, &    !-- oCoordinate_1
    nNodes
  integer ( KDI ), dimension ( :, : ), allocatable :: &
    NodeConnectivity
  real ( KDR ), dimension ( :, : ), allocatable :: &
    NodeCoordinate
  character ( LDF ) :: &
    Name = 'UnstructuredGridImage_Form_Test'
  character ( LDL ), dimension ( 2 ) :: &
    VariableName
  type ( VariableGroupForm ), dimension ( 2 ) :: &
    VG
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( GridImageStreamForm ) :: &
    GIS
  type ( UnstructuredGridImageForm ) :: & 
    UGI_Write, &
    UGI_Read
  
  allocate ( C )
  call C % Initialize ( )
  
  nNodes = 8
  allocate ( NodeCoordinate ( nNodes, 3 ) )
  allocate ( NodeConnectivity ( nNodes, 1 ) )
  
  NodeCoordinate ( :, 1 ) = [ 0., 2., 2., 0., 0., 2., 2., 0. ]
  NodeCoordinate ( :, 2 ) = [ 0., 0., 0., 0., 2., 2., 2., 2. ]
  NodeCoordinate ( :, 3 ) = [ 2., 2., 0., 0., 2., 2., 0., 0. ]
  NodeConnectivity ( :, 1 ) = [ ( iN, iN = 1, nNodes ) ]
    
  oC_1 = C % Rank * maxval ( NodeCoordinate ( :, 1 ) )
  NodeCoordinate ( :, 1 ) = NodeCoordinate ( :, 1 ) + oC_1
  
  VariableName = [ 'VG_1_Name_1', 'VG_1_Name_2' ]
  call VG ( 1 ) % Initialize &
         ( [ 1, 2 ], VariableOption = VariableName, &
           NameOption = 'VariableGroup_1' )
  VG ( 1 ) % Value ( :, : ) = NodeCoordinate ( 6 : 6, 1 : 2 )
  
  VariableName = [ 'VG_2_Name_1', 'VG_2_Name_2' ]
  call VG ( 2 ) % Initialize &
         ( [ 1, 2 ], VariableOption = VariableName, &
           NameOption = 'VariableGroup_2' )
  VG ( 2 ) % Value ( :, : ) = 2 * NodeCoordinate ( 6 : 6, 1 : 2 )
  
  call CONSOLE % Initialize ( C % Rank )

  call GIS % Initialize ( Name, CommunicatorOption = C )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  
  call UGI_Write % Initialize ( GIS )

  call UGI_Write % AddVariableGroup ( VG ( 1 ) )
  call UGI_Write % AddVariableGroup ( VG ( 2 ) )
  
  call UGI_Write % SetGrid  &
         ( Directory = 'Unstructured', &
           NodeCoordinate = NodeCoordinate, &
           NodeConnectivity = NodeConnectivity, &
           nDimensions = 3, nProperCells = 1, nGhostCells = 0, oValue = 0 )
  
  call UGI_Write % Write ()
  
  call GIS % Close ( ) 
  
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = 0 )
  
  call UGI_Read % Initialize ( GIS )
  
  call UGI_Read % SetReadAttributes ( Directory = 'Unstructured', oValue = 0 )
  
  call UGI_Read % Read ( )
  
  call Show ( UGI_Read % nDimensions, 'nDimensions' )  
  call Show ( UGI_Read % nNodes, 'nNodes' )
  call Show ( UGI_Read % nTotalCells, 'nTotalCells' )
  call Show ( UGI_Read % NodeCoordinate_1, 'NodeCoordinate_1' )
  call Show ( UGI_Read % NodeCoordinate_2, 'NodeCoordinate_2' )
  call Show ( UGI_Read % NodeCoordinate_3, 'NodeCoordinate_3' )
  
  call Show ( ( all ( UGI_Write % NodeCoordinate_1 &
                        == UGI_Read % NodeCoordinate_1 ) &
                .and. all ( UGI_Write % NodeCoordinate_2 &
                              == UGI_Read % NodeCoordinate_2 ) &
                .and. all ( UGI_Write % NodeCoordinate_3 &
                              == UGI_Read % NodeCoordinate_3 ) ), &
              'Every R/W NodeCoordinate matches' )
  
  call Show ( ( all ( UGI_Write % VariableGroup ( 1 ) % Value &
                        == UGI_Read % VariableGroup ( 1 ) % Value ) &
                .and. all ( UGI_Write % VariableGroup ( 2 ) % Value &
                              == UGI_Read % VariableGroup ( 2 ) % Value ) ), &
              'Every R/W Variable matches' )
  
  call GIS % Close ( ) 
  
  deallocate ( C )

end program UnstructuredGridImage_Form_Test
