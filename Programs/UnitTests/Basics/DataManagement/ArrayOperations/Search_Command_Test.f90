program Search_Command_Test
  
  use Specifiers
  use Sort_Command
  use Search_Command

  implicit none
  
  integer ( KDI ) :: &
    iA
  integer ( KDI ), dimension ( 10 ) :: &
    A_I
  integer ( KDI ) :: &
    V_I
  real ( KDR ) :: &
    Number, &
    V_R
  real ( KDR ), dimension ( 10 ) :: &
    A_R
    
  print*
  print*, 'Test searching integers'
  
  print* , 'Array:'
  do iA = 1, size ( A_I )
    call random_number ( Number )
    A_I ( iA ) =  floor ( Number * 100 )
  end do
  
  call Sort ( A_I )
  do iA = 1, size ( A_I )
    print*, A_I ( iA )
  end do
  
  print*
  
  V_I = A_I ( 4 )
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  V_I = A_I ( 1 )
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  V_I = A_I ( 10 )
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  V_I = A_I ( 1 ) + 1
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  V_I = A_I ( 1 ) - 1
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  V_I = A_I ( 5 ) + 1
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  V_I = A_I ( 5 ) - 1
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  V_I = A_I ( 10 ) + 1
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  V_I = A_I ( 10 ) - 1
  call Search ( A_I, V_I, iA )
  print* , 'Searched : ', V_I
  print* , 'Location : ', iA
print*
  
  print*
  print*, 'Test searching reals'
  
  print* , 'Array:'
  do iA = 1, size ( A_R )
    call random_number ( Number )
    A_R ( iA ) =  Number * 10
  end do
  
  call Sort ( A_R )
  do iA = 1, size ( A_R )
    print*, A_R ( iA )
  end do
  
  print*
  
  V_R = A_R ( 4 )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
  V_R = A_R ( 1 )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
  V_R = A_R ( 10 )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
  V_R = A_R ( 1 ) + 10.0_KDR * epsilon ( 1.0_KDR )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
  V_R = A_R ( 1 ) - 10.0_KDR * epsilon ( 1.0_KDR )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
  V_R = A_R ( 5 ) + 10.0_KDR * epsilon ( 1.0_KDR )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
  V_R = A_R ( 5 ) - 10.0_KDR * epsilon ( 1.0_KDR )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
  V_R = A_R ( 10 ) + 10.0_KDR * epsilon ( 1.0_KDR )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
  V_R = A_R ( 10 ) - 10.0_KDR * epsilon ( 1.0_KDR )
  call Search ( A_R, V_R, iA )
  print* , 'Searched : ', V_R
  print* , 'Location : ', iA
  print*
  
end program Search_Command_Test
