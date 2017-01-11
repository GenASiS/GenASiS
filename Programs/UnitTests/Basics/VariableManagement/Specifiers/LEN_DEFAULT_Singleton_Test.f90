program LEN_DEFAULT_Singleton_Test

  use LEN_DEFAULT_Singleton

  implicit none

  print *
  print *, 'LEN_DEFAULT % NUMBER   = ', LEN_DEFAULT % NUMBER
  print *, 'LEN_DEFAULT % LABEL    = ', LEN_DEFAULT % LABEL
  print *, 'LEN_DEFAULT % FILENAME = ', LEN_DEFAULT % FILENAME
  print *, 'LEN_DEFAULT % BUFFER   = ', LEN_DEFAULT % BUFFER

  print *
  print *, 'LDN = ', LDN
  print *, 'LDL = ', LDL
  print *, 'LDF = ', LDF
  print *, 'LDB = ', LDB

end program LEN_DEFAULT_Singleton_Test
