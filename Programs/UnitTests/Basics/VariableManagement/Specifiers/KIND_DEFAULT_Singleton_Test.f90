program KIND_DEFAULT_Singleton_Test

  use KIND_DEFAULT_Singleton

  implicit none

  print *
  print *, 'KIND_DEFAULT % INTEGER = ', KIND_DEFAULT % INTEGER
  print *, 'KIND_DEFAULT % REAL    = ', KIND_DEFAULT % REAL
  print *, 'KIND_DEFAULT % COMPLEX = ', KIND_DEFAULT % COMPLEX
  print *, 'KIND_DEFAULT % LOGICAL = ', KIND_DEFAULT % LOGICAL

  print *
  print *, 'KDI = ', KDI
  print *, 'KDR = ', KDR
  print *, 'KDC = ', KDC
  print *, 'KDL = ', KDL

  print *
  print *, 'huge ( 0_KDI ) = ', huge ( 0_KDI )

  print *
  print *, 'tiny ( 0.0_KDR ) = ', tiny ( 0.0_KDR )
  print *, 'huge ( 0.0_KDR ) = ', huge ( 0.0_KDR )

end program KIND_DEFAULT_Singleton_Test
