program REDUCTION_Singleton_Test

  use VariableManagement
  use Display
  use MessagePassingBasics
  use REDUCTION_Singleton

  implicit none

  type ( CommunicatorForm ), allocatable :: &
    C

  allocate ( C )

  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetDisplayRank ( 0 )
  
  call Show & 
         ( REDUCTION % MAX, 'REDUCTION % MAX' )
  call Show & 
         ( REDUCTION % MIN, 'REDUCTION % MIN' )
  call Show & 
         ( REDUCTION % SUM, 'REDUCTION % SUM' )
  call Show & 
         ( REDUCTION % PRODUCT, 'REDUCTION % PRODUCT' )
  call Show & 
         ( REDUCTION % LOGICAL_AND, 'REDUCTION % LOGICAL_AND' )
  call Show & 
         ( REDUCTION % BITWISE_AND, 'REDUCTION % BITWISE_AND' )
  call Show & 
         ( REDUCTION % LOGICAL_OR, 'REDUCTION % LOGICAL_OR' )
  call Show & 
         ( REDUCTION % BITWISE_OR, 'REDUCTION % BITWISE_OR' )
  call Show & 
         ( REDUCTION % LOGICAL_EXCLUSIVE_OR, &
           'REDUCTION % LOGICAL_EXCLUSIVE_OR' )
  call Show & 
         ( REDUCTION % BITWISE_EXCLUSIVE_OR, &
           'REDUCTION % BITWISE_EXCLUSIVE_OR' )
  call Show & 
         ( REDUCTION % MIN_LOCATION, 'REDUCTION % MIN_LOCATION' )
  call Show & 
         ( REDUCTION % MAX_LOCATION, 'REDUCTION % MAX_LOCATION' )

  deallocate ( C )

end program REDUCTION_Singleton_Test
