program PointGridImage_Form_Test

  use VariableManagement
  use Display
  use MessagePassing
  use GridImageStream_Form
  use PointGridImage_Form
  
  implicit none
  
  integer ( KDI ), parameter :: &
    N_DIMENSIONS = 3, &
    N_CELLS  = 100
  integer ( KDI ) :: &
    iC          !-- iCell
  real ( KDR ) :: &
    T_0, &
    Angle
  real ( KDR ), dimension ( N_CELLS, 3 ) :: &
    Coordinate
  character ( LDF ) :: &
    Name = 'PointGridImage_Form_Test'
  character ( LDL ), dimension ( 3 ) :: &
    VariableName
  type ( VariableGroupForm ) :: &
    VG
    type ( CommunicatorForm ), allocatable :: &
    C
  type ( GridImageStreamForm ) :: &
    GIS
  type ( PointGridImageForm ) :: & 
    PGI_Write, &
    PGI_Read
  
  allocate ( C )
  call C % Initialize ( )
  
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetVerbosity ( 'INFO_7' )
  
  do iC = 1, N_CELLS
    T_0 = real ( iC - 1, kind = KDR ) / real ( N_CELLS - 1, kind = KDR )
    Angle = CONSTANT % PI * 10.0_KDR * T_0
    Coordinate ( iC, 1 ) = T_0 * cos ( Angle )
    Coordinate ( iC, 2 ) = T_0 * sin ( Angle )
    Coordinate ( iC, 3 ) = T_0 
  end do
  
  VariableName ( 1 ) = 'Center_1'
  VariableName ( 2 ) = 'Center_2'
  VariableName ( 3 ) = 'Center_3'
  
  call VG % Initialize &
         ( [ N_CELLS, 3 ], VariableOption = VariableName, &
           NameOption = 'Geometry' )
  VG % Value = Coordinate
  
  call GIS % Initialize ( Name, CommunicatorOption = C )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  
  call PGI_Write % Initialize ( GIS )

  call PGI_Write % AddVariableGroup ( VG )
  
  call PGI_Write % SetGrid  &
         ( Directory = 'Grid', &
           Coordinate = Coordinate, nCells = N_CELLS, &
           nDimensions = N_DIMENSIONS, oValue = 0 )

  call Show ( PGI_Write % nDimensions,      'nDimensions'  )  
  call Show ( PGI_Write % nCells,           'nCells'       )
  call Show ( PGI_Write % NodeCoordinate_1, 'Coordinate_1' )
  call Show ( PGI_Write % NodeCoordinate_2, 'Coordinate_2' )
  call Show ( PGI_Write % NodeCoordinate_3, 'Coordinate_3' )
  call Show ( VG % Value ( :, 1 ),          'Center_1' )
  call Show ( VG % Value ( :, 2 ),          'Center_1' )
  call Show ( VG % Value ( :, 3 ),          'Center_1' )
  
  call PGI_Write % Write ( )
  
  call GIS % Close ( ) 
  
  call Clear ( VG % Value )
  
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = 0 )
  
  call PGI_Read % Initialize ( GIS )
  
  call PGI_Read % SetReadAttributes ( Directory = 'Grid', oValue = 0 )
  
  call PGI_Read % Read ( )
  
  associate ( VG_Read => PGI_Read % VariableGroup ( 1 ) )
  
  call Show ( PGI_Read % nDimensions,      'nDimensions'  )  
  call Show ( PGI_Read % nCells,           'nCells'       )
  call Show ( PGI_Read % NodeCoordinate_1, 'Coordinate_1' )
  call Show ( PGI_Read % NodeCoordinate_2, 'Coordinate_2' )
  call Show ( PGI_Read % NodeCoordinate_3, 'Coordinate_3' )
  call Show ( VG_Read % Value ( :, 1 ),    'Center_1' )
  call Show ( VG_Read % Value ( :, 2 ),    'Center_1' )
  call Show ( VG_Read % Value ( :, 3 ),    'Center_1' )
  
  end associate

  call GIS % Close ( )
  
  deallocate ( C )

end program PointGridImage_Form_Test
