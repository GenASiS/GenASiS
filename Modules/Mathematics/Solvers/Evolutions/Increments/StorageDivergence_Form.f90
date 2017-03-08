!-- StorageDivergence contains storage for an IncrementDivergence.

module StorageDivergence_Form

  use Basics
  use Operations

  implicit none
  private

  type, public :: StorageDivergenceForm
    integer ( KDI ) :: &
      iTimerDivergence     = 0, &
      iTimerReconstruction = 0, &
      iTimerFluxes         = 0, &
      iTimerIncrement      = 0
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
      InitializeTimers
    procedure, public, pass :: &
      Allocate
    procedure, public, pass :: &
      Deallocate
    final :: &
      Finalize
  end type StorageDivergenceForm

contains


  subroutine InitializeTimers ( SD, BaseLevel )

    class ( StorageDivergenceForm ), intent ( inout ) :: &
      SD
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    if ( SD % iTimerDivergence > 0  &
         .or.  BaseLevel > PROGRAM_HEADER % TimerLevel ) &
      return

    call PROGRAM_HEADER % AddTimer &
           ( 'Divergence', SD % iTimerDivergence, &
             Level = BaseLevel )
      call PROGRAM_HEADER % AddTimer &
             ( 'Reconstruction', SD % iTimerReconstruction, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'Fluxes', SD % iTimerFluxes, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'Increment', SD % iTimerIncrement, &
               Level = BaseLevel + 1 )

  end subroutine InitializeTimers


  subroutine Allocate &
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

  end subroutine Allocate


  subroutine Deallocate ( SD )

    class ( StorageDivergenceForm ), intent ( inout ) :: &
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

  end subroutine Deallocate


  subroutine Finalize ( SD )

    type ( StorageDivergenceForm ), intent ( inout ) :: &
      SD

    call SD % Deallocate ( )

  end subroutine Finalize


end module StorageDivergence_Form
