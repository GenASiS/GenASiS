program KIND_DEFAULT_Singleton_Test

  use KIND_DEFAULT_Singleton

  implicit none

  integer ( KDI ), parameter :: &
    PRINTABLE_MIN = 32, &
    PRINTABLE_MAX = 126, &
    CHARACTERS_PER_ROW  = 32

  integer ( KDI ) :: &
    iV, &  !-- iValue
    iC     !-- iCharacter
  character ( 80 ) :: &
    Row

  print *
  print *, 'KIND_DEFAULT % INTEGER   = ', KIND_DEFAULT % INTEGER
  print *, 'KIND_DEFAULT % REAL      = ', KIND_DEFAULT % REAL
  print *, 'KIND_DEFAULT % COMPLEX   = ', KIND_DEFAULT % COMPLEX
  print *, 'KIND_DEFAULT % LOGICAL   = ', KIND_DEFAULT % LOGICAL
  print *, 'KIND_DEFAULT % CHARACTER = ', KIND_DEFAULT % CHARACTER

  print *
  print *, 'KDI  =  ', KDI
  print *, 'KDR  =  ', KDR
  print *, 'KDC  =  ', KDC
  print *, 'KDL  =  ', KDL
  print *, 'KDCH =  ', KDCH

  print *
  print *, 'huge ( 0_KDI ) = ', huge ( 0_KDI )

  print *
  print *, 'tiny ( 0.0_KDR ) = ', tiny ( 0.0_KDR )
  print *, 'huge ( 0.0_KDR ) = ', huge ( 0.0_KDR )

  print *
  print *, 'ASCII Printable'

  iC  = 0
  Row = ''
  do iV = PRINTABLE_MIN, PRINTABLE_MAX

    iC = iC + 1
    Row ( 2 * iC : 2 * iC ) = char ( iV, kind = KDCH )

    if ( iC == CHARACTERS_PER_ROW .or. iV == PRINTABLE_MAX ) then
      print *, Row
      iC  = 0
      Row = ''
    end if

  end do !-- iV

end program KIND_DEFAULT_Singleton_Test
