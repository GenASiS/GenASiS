!-- LEN_DEFAULT_Singleton defines some standard len parameters for the
!   character data type. 

module LEN_DEFAULT_Singleton

  use KIND_DEFAULT_Singleton

  implicit none
  private

  type, public :: LenDefaultSingleton
    integer ( kind ( KDI ) ) :: &
      NUMBER   = 7, &
      LABEL    = 31, &
      FILENAME = 255, &
      BUFFER   = 1023
  end type LenDefaultSingleton

  type ( LenDefaultSingleton ), public, parameter :: &
    LEN_DEFAULT = LenDefaultSingleton ( )

  integer ( KDI ), public, parameter :: &
    LDN = LEN_DEFAULT % NUMBER, &
    LDL = LEN_DEFAULT % LABEL, &
    LDF = LEN_DEFAULT % FILENAME, &
    LDB = LEN_DEFAULT % BUFFER

end module LEN_DEFAULT_Singleton
