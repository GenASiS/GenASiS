program KIND_TINY_Singleton_Test

  use KIND_TINY_Singleton

  implicit none

  print *
  print *, 'KIND_TINY % INTEGER = ', KIND_TINY % INTEGER
  print *, 'KIND_TINY % LOGICAL = ', KIND_TINY % LOGICAL

  print *
  print *, 'KTI = ', KTI
  print *, 'KTL = ', KTL

  print *
  print *, 'huge ( 0_KTI ) = ', huge ( 0_KTI )

end program KIND_TINY_Singleton_Test
