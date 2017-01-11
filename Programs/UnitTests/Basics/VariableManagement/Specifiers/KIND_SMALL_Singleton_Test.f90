program KIND_SMALL_Singleton_Test

  use KIND_SMALL_Singleton

  implicit none

  print *
  print *, 'KIND_SMALL % INTEGER = ', KIND_SMALL % INTEGER

  print *
  print *, 'KSI = ', KSI

  print *
  print *, 'huge ( 0_KSI ) = ', huge ( 0_KSI )

end program KIND_SMALL_Singleton_Test
