program KIND_BIG_Singleton_Test

  use iso_fortran_env
  use KIND_DEFAULT_Singleton
  use KIND_BIG_Singleton

  implicit none

  integer ( KDI ), parameter :: &
    N_BLOCKS = 7, &
    CHARACTERS_PER_ROW  = 32
  integer ( KDI ), dimension ( N_BLOCKS ), parameter :: &
    BLOCK_MIN = [ int ( z'0020' ), &
                  int ( z'00A0' ), &
                  int ( z'0100' ), &
                  int ( z'2070' ), &
                  int ( z'2200' ), &
                  int ( z'2600' ), &
                  int ( z'1D400' ) ], &
    BLOCK_MAX = [ int ( z'007E' ), &
                  int ( z'00FF' ), &
                  int ( z'017F' ), &
                  int ( z'209F' ), &
                  int ( z'22FF' ), &
                  int ( z'26FF' ), &
                  int ( z'1D7FF' ) ]
  character ( 64 ), dimension ( N_BLOCKS ), parameter :: &
    BLOCK_NAME = [ 'Basic Latin                      ', &
                   'Latin-1 Supplement               ', &
                   'Latin Extended-A                 ', &
                   'Superscripts and Subscripts      ', &
                   'Mathematical Operators           ', &
                   'Miscellaneous Symbols            ', &
                   'Mathematical Alphanumeric Symbols' ]

  integer ( KDI ) :: &
    iB, &  !-- iBlock
    iV, &  !-- iValue
    iC     !-- iCharacter
  character ( 80, kind = KBCH ) :: &
    Row
  character ( 5 ) :: &
    Encoding

!-- Runtime error with CCE
!  if ( KBCH == selected_char_kind ( 'ASCII' ) ) then
!print*, '>>> 1'
!    open ( OUTPUT_UNIT, encoding = 'DEFAULT' )
!print*, '>>> 2'
!  else if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
  if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
    Encoding = 'UTF-8'
    open ( OUTPUT_UNIT, encoding = Encoding )
  end if

  print *
  print *, 'KIND_BIG % INTEGER   = ', KIND_BIG % INTEGER
  print *, 'KIND_BIG % CHARACTER = ', KIND_BIG % CHARACTER

  print *
  print *, 'KBI  = ', KBI
  print *, 'KBCH = ', KBCH

  print *
  print *, 'huge ( 0_KBI ) = ', huge ( 0_KBI )

  do iB = 1, N_BLOCKS

    print *
    print *, 'Unicode ' // BLOCK_NAME ( iB )

    iC  = 0
    Row = ''
    do iV = BLOCK_MIN ( iB ), BLOCK_MAX ( iB )

      iC = iC + 1
      Row ( 2 * iC : 2 * iC ) = char ( iV, kind = KBCH )

      if ( iC == CHARACTERS_PER_ROW .or. iV == BLOCK_MAX ( iB ) ) then
        print *, Row
        iC  = 0
        Row = ''
      end if

    end do !-- iV
  end do !-- iB

end program KIND_BIG_Singleton_Test
