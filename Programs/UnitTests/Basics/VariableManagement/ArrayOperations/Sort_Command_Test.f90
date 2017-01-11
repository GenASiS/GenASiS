program Sort_Command_Test
  
  use Specifiers
  use Sort_Command

  implicit none
  
  integer ( KDI ) :: &
    iA
  integer ( KDI ), dimension ( 10 ) :: &
    A_I
  real ( KDR ), dimension ( 10 ) :: &
    A_R
  real ( KDR ) :: &
    Number
    
  print *
  print *, 'Test sorting reals'
  print *, 'Unsorted:'
  do iA = 1, size ( A_R )
    call random_number ( Number )
    A_R ( iA ) =  Number * 10
    print *, A_R ( iA )
  end do
  
  print *
  
  call Sort ( A_R )
  
  print * , 'Sorted: '
  do iA = 1, size ( A_R )
    print *, A_R ( iA )
  end do
  
  print *
  print *, 'Test sorting integers'  
  print * , 'Unsorted:'
  do iA = 1, size ( A_I )
    call random_number ( Number )
    A_I ( iA ) =  floor ( Number * 10 )
    print *, A_I ( iA )
  end do
  
  print*
  
  call Sort ( A_I )
  
  print* , 'Sorted: '
  do iA = 1, size ( A_I )
    print*, A_I ( iA )
  end do
  
end program Sort_Command_Test
