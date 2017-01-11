program GridImageStream_Form_Test

  use VariableManagement
  use Display
  use MessagePassing
  use GridImageStream_Form

  implicit none
  
  integer ( KDI ) :: &
    iStep
  character ( LDL ), dimension ( 6 ) :: &
    Variable &
      = [ 'Variable_1', 'Variable_2', 'Variable_3', 'Variable_4', &
          'Variable_5', 'Variable_6'  ]
  character ( LDF ) :: &
    Name = 'Stream_Form_Test'
  type ( VariableGroupForm ), dimension ( : ), allocatable :: &
    VG_W, &
    VG_R
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  
  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetVerbosity ( 'INFO_7' )

  allocate ( GIS )
  call GIS % Initialize ( Name, CommunicatorOption = C )
  
  do iStep = 1, 5

    call GIS % Open ( GIS % ACCESS_CREATE )
    
    call Show ( GIS % IsWritable ( ), 'Is Stream Writable ?' )
    call Show ( GIS % IsReadable ( ), 'Is Stream Readable ?' )
    
    call GIS % MakeDirectory ( 'SubDirectory_1' )
    call GIS % MakeDirectory ( 'SubSubDirectory_1_1' )
    call GIS % MakeDirectory ( 'SubSubSubDirectory_1_1_1' )
    
    call GIS % ChangeDirectory ( '../../' )
    call GIS % MakeDirectory ( 'SubSubDirectory_1_2/' )
    call GIS % ChangeDirectory ( '../' )
    call GIS % MakeDirectory ( 'SubSubDirectory_1_3/' )
    
    call GIS % ChangeDirectory ( '/' )
    call GIS % MakeDirectory ( 'SubDirectory_2' )
    call GIS % MakeDirectory ( 'SubSubDirectory_2_1' )
    call GIS % MakeDirectory ( 'SubSubSubDirectory_2_1_1' )
    
    call GIS % ChangeDirectory ( '.' )
    call GIS % ChangeDirectory ( '..' )
    call GIS % MakeDirectory ( 'SubSubSubDirectory_2_1_2' )
    
    call GIS % ChangeDirectory ( '../..' )
    call GIS % MakeDirectory ( 'SubSubDirectory_2_2' )
    
    call GIS % ChangeDirectory ( '/SubDirectory_1/SubSubDirectory_1_1')
    call GIS % MakeDirectory ( 'SubSubSubDirectory_1_1_2' )
    
    call GIS % ChangeDirectory ( '../..')
    call GIS % ChangeDirectory ( 'SubSubDirectory_1_3')
    call GIS % MakeDirectory ( 'SubSubSubDirectory_1_3_1' )
    
    call GIS % Close ( )
    
    call Show ( GIS % IsWritable ( ), 'Is Stream Writable ?' )
    call Show ( GIS % IsReadable ( ), 'Is Stream Readable ?' )
    
    call GIS % Open ( GIS % ACCESS_WRITE )
    
    call Show ( GIS % IsWritable ( ), 'Is Stream Writable ?' )
    call Show ( GIS % IsReadable ( ), 'Is Stream Readable ?' )
    
    call GIS % MakeDirectory ( 'SubDirectory_3' )
    call GIS % ChangeDirectory ( '/' )
    call GIS % MakeDirectory ( 'SubDirectory_4' )
    
    allocate ( VG_W ( 3 ) )
    call VG_W ( 1 ) % Initialize &
           ( ValueShape = [ 10, 6 ], &
             VariableOption = Variable, NameOption = 'VariableGroup_1' )
    call random_number ( VG_W ( 1 ) % Value )

    call VG_W ( 2 ) % Initialize &
           ( ValueShape = [ 15, 6 ], &
             VariableOption = Variable, NameOption = 'VariableGroup_2')
    call random_number ( VG_W ( 2 ) % Value )
    
    call VG_W ( 3 ) % Initialize &
           ( ValueShape = [ 26, 6 ], &
             VariableOption = Variable, NameOption = 'VariableGroup_3')
    call random_number ( VG_W ( 3 ) % Value )
    
    call GIS % ChangeDirectory ( '/' )
    call GIS % WriteVariableGroup &
           ( VG_W, DirectoryOption = 'SubDirectory_1' )
    
    call GIS % Close ( )
    
    call GIS % Open ( GIS % ACCESS_READ, NumberOption = iStep - 1 )
    
    call Show ( GIS % IsWritable ( ), 'Is Stream Writable ?' )
    call Show ( GIS % IsReadable ( ), 'Is Stream Readable ?' )
    
    call GIS % ListContents ( ContentTypeOption = 'Directory' )
    call Show ( size ( GIS % ContentList ), 'nDirectories' )
    
    allocate ( VG_R ( 3 ) )
    call GIS % ReadVariableGroup &
           ( VG_R, DirectoryOption = 'SubDirectory_1' )
           
    call Show &
           ( all ( VG_W ( 1 ) % Value == VG_R ( 1 ) % Value ) &
               .and. all ( VG_W ( 2 ) % Value == VG_R ( 2 ) % Value ) &
               .and. all ( VG_W ( 3 ) % Value == VG_R ( 3 ) % Value ), &
             'Ever R/W VG Value Matches' )
    
    call GIS % Close ( )
    
    deallocate ( VG_R )
    deallocate ( VG_W )
    
  end do
  
  call GIS % Open ( GIS % ACCESS_WRITE, NumberOption = 0 )
  call GIS % MakeDirectory ( 'SubDirectory_5' )
  
  call GIS % ChangeDirectory ( '/')
  call GIS % ListContents ( ContentTypeOption = 'Directory' )
  call Show ( size ( GIS % ContentList ), 'nDirectories' )
  
  call GIS % Close ( ) 
  
  call GIS % Open ( GIS % ACCESS_CREATE ) 
  call GIS % MakeDirectory ( 'SubDirectory_1' )
  
  call GIS % Close ( ) 

  deallocate ( GIS )
  deallocate ( C )

end program GridImageStream_Form_Test
