!-- ATLAS contains flags used in connection with Atlas classes.

module ATLAS_Singleton

  use Basics

  implicit none
  private

  type, public :: AtlasSingleton
    integer ( KDI ) :: &
      MAX_DIMENSIONS = 3, &
      MAX_CHARTS     = 8, &
      MAX_FIELDS     = 96, &
      MAX_STREAMS    = 8
  end type AtlasSingleton

  type ( AtlasSingleton ), public, parameter :: &
    ATLAS = AtlasSingleton ( )

end module ATLAS_Singleton
