program CurveImage_Form_Test

  use VariableManagement
  use Display
  use MessagePassing
  use GridImageStream_Form
  use CurveImage_Form
  
  implicit none
  
  integer ( KDI ) :: &
    iC, &
    oC_1, &
    DisplayRank = 0
  real ( KDR ), dimension ( 20 ) :: &
    NodeCoordinate
  character ( LDF ) :: &
    Name = 'CurveImage_Form_Test'
  character ( LDL ), dimension ( 1 ) :: &
    VariableName
  type ( VariableGroupForm ) :: &
    VG
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( GridImageStreamForm ) :: &
    GIS
  type ( CurveImageForm ) :: & 
    CI_Write, &
    CI_Read
  
  allocate ( C )
  call C % Initialize ( )
  
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetDisplayRank ( DisplayRank )

  NodeCoordinate =  [ ( iC * 0.5_KDR, iC = 1, 20 ) ]
  
  oC_1 = C % Rank * maxval ( NodeCoordinate )
  NodeCoordinate = NodeCoordinate  + oC_1
  
  VariableName ( 1 ) = 'VariableName'
  call VG % Initialize &
         ( [ 20, 1 ], VariableOption = VariableName, &
           NameOption = 'VariableGroup' )
  VG % Value ( :, 1 ) &
    = [ ( ( iC + 20 * C % Rank ) * 2.0_KDR, iC = 1, 20 ) ] 
    
  call Show ( NodeCoordinate, 'NodeCoordinate' )
  call Show ( VG % Value ( :, 1 ), 'Variable' )

  call GIS % Initialize ( Name, CommunicatorOption = C )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  
  call CI_Write % Initialize ( GIS )
  
  call CI_Write % AddVariableGroup ( VG )

  call CI_Write % SetGrid  &
         ( Directory = 'Curves', NodeCoordinate = NodeCoordinate, &
           nProperCells = 20, oValue = 0 )

  call CI_Write % Write ()

  call GIS % Close ( ) 
  
  call Clear ( VG % Value )
  
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = 0 )
  
  call CI_Read % Initialize ( GIS )
  
  call CI_Read % SetReadAttributes ( Directory = 'Curves', oValue = 0 )
  
  call CI_Read % Read ( )
  
  call Show ( CI_Read % NodeCoordinate_1, 'NodeCoordinate_1' )
  call Show &
         ( CI_Read % VariableGroup ( 1 ) % Value ( :, 1 ), &
           CI_Read % VariableGroup ( 1 ) % Variable ( 1 ) )
  
  call GIS % Close ( )

  deallocate ( C )
  
end program CurveImage_Form_Test
