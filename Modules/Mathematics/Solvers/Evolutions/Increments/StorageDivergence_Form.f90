!-- StorageDivergence contains storage for an IncrementDivergence.

module StorageDivergence_Form

  use Basics
  use Operations

  implicit none
  private

  type, public :: StorageDivergenceForm
    type ( VariableGroupForm ), allocatable :: &
      Geometry_I, &              !-- Geometry_Inner
      Current_IL, Current_IR, &  !-- Current_InnerLeft, Current_InnerRight
      ModifiedSpeeds_I, &        !-- ModifiedSpeeds_Inner
      DiffusionFactor_I, &       !-- DiffusionFactor_Inner
      Flux_IL, Flux_IR, &        !-- Flux_InnerLeft, Flux_InnerRight
      Flux_I                     !-- Flux_Inner
    type ( GradientForm ), allocatable :: &
      GradientPrimitive
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type StorageDivergenceForm

contains


  subroutine Initialize &
              ( SD, nCurrent, nConserved, nPrimitive, nModifiedSpeeds, &
                nGeometry, nValues )

    class ( StorageDivergenceForm ), intent ( inout ) :: &
      SD
    integer ( KDI ), intent ( in ) :: &
      nCurrent, &
      nConserved, &
      nPrimitive, &
      nModifiedSpeeds, &
      nGeometry, &
      nValues

    allocate ( SD % Geometry_I )
    call SD % Geometry_I % Initialize &
           ( [ nValues, nGeometry ], ClearOption = .true. )

    allocate ( SD % Current_IL, SD % Current_IR )
    call SD % Current_IL % Initialize &
           ( [ nValues, nCurrent ], ClearOption = .true. )
    call SD % Current_IR % Initialize &
           ( [ nValues, nCurrent ], ClearOption = .true. )

    allocate ( SD % ModifiedSpeeds_I )
    call SD % ModifiedSpeeds_I % Initialize &
           ( [ nValues, nModifiedSpeeds ], ClearOption = .true. )

    allocate ( SD % DiffusionFactor_I )
    call SD % DiffusionFactor_I % Initialize &
           ( [ nValues, nConserved ], ClearOption = .true. )

    allocate ( SD % Flux_IL, SD % Flux_IR )
    call SD % Flux_IL % Initialize &
           ( [ nValues, nConserved ], ClearOption = .true. )
    call SD % Flux_IR % Initialize &
           ( [ nValues, nConserved ], ClearOption = .true. )

    allocate ( SD % Flux_I )
    call SD % Flux_I % Initialize &
           ( [ nValues, nConserved ], ClearOption = .true. )

    allocate ( SD % GradientPrimitive )
    call SD % GradientPrimitive % Initialize &
           ( 'Primitive', [ nValues, nPrimitive ] )

  end subroutine Initialize


  subroutine Finalize ( SD )

    type ( StorageDivergenceForm ), intent ( inout ) :: &
      SD

    if ( allocated ( SD % GradientPrimitive ) ) &
      deallocate ( SD % GradientPrimitive )
    if ( allocated ( SD % DiffusionFactor_I ) ) &
      deallocate ( SD % DiffusionFactor_I )
    if ( allocated ( SD % ModifiedSpeeds_I ) ) &
      deallocate ( SD % ModifiedSpeeds_I )
    if ( allocated ( SD % Flux_I ) ) &
      deallocate ( SD % Flux_I )
    if ( allocated ( SD % Current_IR ) ) &
      deallocate ( SD % Current_IR )
    if ( allocated ( SD % Current_IL ) ) &
      deallocate ( SD % Current_IL )
    if ( allocated ( SD % Geometry_I ) ) &
      deallocate ( SD % Geometry_I )

  end subroutine Finalize


end module StorageDivergence_Form
