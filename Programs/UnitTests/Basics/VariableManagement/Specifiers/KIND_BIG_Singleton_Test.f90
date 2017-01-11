program KIND_BIG_Singleton_Test

  use KIND_BIG_Singleton

  implicit none

  print *
  print *, 'KIND_BIG % INTEGER = ', KIND_BIG % INTEGER

  print *
  print *, 'KBI = ', KBI

  print *
  print *, 'huge ( 0_KBI ) = ', huge ( 0_KBI )

end program KIND_BIG_Singleton_Test
