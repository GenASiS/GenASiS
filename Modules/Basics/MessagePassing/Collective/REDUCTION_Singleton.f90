!-- ReductionSingleton abstracts specified constants and flags for collective
!   operation.

module REDUCTION_Singleton

  use MPI
  use Specifiers

  implicit none
  private
  
  type, public :: ReductionSingleton
    integer ( KDI ), public :: &
      MAX                  = MPI_MAX, &
      MIN                  = MPI_MIN, &
      SUM                  = MPI_SUM, &
      PRODUCT              = MPI_PROD, &
      LOGICAL_AND          = MPI_LAND, &
      BITWISE_AND          = MPI_BAND, &
      LOGICAL_OR           = MPI_LOR, &
      BITWISE_OR           = MPI_BOR, &
      LOGICAL_EXCLUSIVE_OR = MPI_LXOR, &
      BITWISE_EXCLUSIVE_OR = MPI_BXOR, &
      MIN_LOCATION         = MPI_MINLOC, &
      MAX_LOCATION         = MPI_MAXLOC
  end type ReductionSingleton
  
  type ( ReductionSingleton ), public, parameter :: &
    REDUCTION = ReductionSingleton ( )

end module REDUCTION_Singleton
