program SplitJoin_Command_Test
  
  use KIND_DEFAULT_Singleton
  use LEN_DEFAULT_Singleton
  use Split_Command
  use Join_Command

  implicit none
  
  integer ( KDI ) :: &
    iW  !-- iWord
  character ( LDL ), dimension ( : ), allocatable :: &
    Word
  character ( LDB ) :: &
    Buffer
  
  print *
  Buffer = 'To boldly go where no man has gone before'
  call Split ( Buffer, ' ', Word )
  print *, 'Buffer: ', trim ( Buffer )
  print *, 'Words : '
  do iW = 1, size ( Word )
    print *, iW, Word ( iW )
  end do
  Buffer = '' 
  call Join ( Word, ' ', Buffer )
  print *, 'Joined: ', trim ( Buffer )

  print *  
  Buffer = 'To  ,  boldly ,  go ,  where ,  no ,  man ,  has ,  gone ,  before'
  call Split ( Buffer, ',', Word )
  print *, 'Buffer: ', trim ( Buffer )
  print *, 'Words : '
  do iW = 1, size ( Word )
    print *, iW, Word ( iW )
  end do
  Buffer = '' 
  call Join ( Word, ' ,  ', Buffer )
  print *, 'Joined: ', trim ( Buffer )

  print *  
  Buffer = 'To@@@boldly@@@go@@@where@@@no@@@man@@@has@@@gone@@@before'
  call Split ( Buffer, '@@@', Word )
  print *, 'Buffer: ', trim ( Buffer )
  print *, 'Words : '
  do iW = 1, size ( Word )
    print *, iW, Word ( iW )
  end do
  Buffer = '' 
  call Join ( Word, '@@@', Buffer )
  print *, 'Joined: ', trim ( Buffer )
  
end program SplitJoin_Command_Test
