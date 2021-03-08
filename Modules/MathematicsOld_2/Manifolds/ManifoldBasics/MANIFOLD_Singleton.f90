!-- MANIFOLD contains flags used in connection with Manifold classes.

module MANIFOLD_Singleton

  use Basics

  implicit none
  private

  type, public :: ManifoldSingleton
    integer ( KDI ) :: &
      MAX_DIMENSIONS = 3, &
      MAX_CHARTS     = 8, &
      MAX_FIELDS     = 96, &
      MAX_STREAMS    = 8
  end type ManifoldSingleton

  type ( ManifoldSingleton ), public, parameter :: &
    MANIFOLD = ManifoldSingleton ( )

end module MANIFOLD_Singleton
