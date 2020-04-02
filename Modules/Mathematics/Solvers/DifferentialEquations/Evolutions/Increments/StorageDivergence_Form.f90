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
      iTimerIncrement      = 0, &
      iTimerDataToDevice   = 0, &
      iTimerDataToHost     = 0
    type ( StorageForm ), allocatable :: &
      Geometry_I, &              !-- Geometry_Inner
      Current_IL, Current_IR, &  !-- Current_InnerLeft, Current_InnerRight
      SolverSpeeds_I, &          !-- SolverSpeeds_Inner
      DiffusionFactor_I, &       !-- DiffusionFactor_Inner
      Flux_IL, Flux_IR, &        !-- Flux_InnerLeft, Flux_InnerRight
      Flux_I, &                  !-- Flux_Inner
      Current_ICL, Current_ICR   !-- Current_InnerCenterLeft, _InnerCenterRight
    type ( GradientForm ), allocatable :: &
      GradientReconstructed
  contains
    procedure, public, pass :: &
      InitializeTimers
    procedure, public, pass :: &
      Allocate
    procedure, public, pass :: &
      Allocate_HLLC
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
      call PROGRAM_HEADER % AddTimer &
             ( 'SD_DataToDevice', SD % iTimerDataToDevice, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'SD_DataToHost', SD % iTimerDataToHost, &
               Level = BaseLevel + 1 )

  end subroutine InitializeTimers


  subroutine Allocate &
              ( SD, AllocateDevice, nCurrent, nConserved, nReconstructed, &
                nSolverSpeeds, nGeometry, nValues )

    class ( StorageDivergenceForm ), intent ( inout ) :: &
      SD
    logical ( KDL ), intent ( in ) :: &
      AllocateDevice
    integer ( KDI ), intent ( in ) :: &
      nCurrent, &
      nConserved, &
      nReconstructed, &
      nSolverSpeeds, &
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

    allocate ( SD % SolverSpeeds_I )
    call SD % SolverSpeeds_I % Initialize &
           ( [ nValues, nSolverSpeeds ], ClearOption = .true. )

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

    allocate ( SD % GradientReconstructed )
    call SD % GradientReconstructed % Initialize &
           ( 'Reconstructed', [ nValues, nReconstructed ] )
           
    if ( AllocateDevice ) then
      call SD % Geometry_I % AllocateDevice ( )
      call SD % Current_IL % AllocateDevice ( )
      call SD % Current_IR % AllocateDevice ( )
      call SD % SolverSpeeds_I % AllocateDevice ( )
      call SD % DiffusionFactor_I % AllocateDevice ( )
      call SD % Flux_IL % AllocateDevice ( )
      call SD % Flux_IR % AllocateDevice ( )
      call SD % Flux_I % AllocateDevice ( )
      call SD % GradientReconstructed % AllocateDevice ( )
    end if

  end subroutine Allocate


  subroutine Allocate_HLLC ( SD, AllocateDevice, nCurrent, nValues )

    class ( StorageDivergenceForm ), intent ( inout ) :: &
      SD
    logical ( KDL ), intent ( in ) :: &
      AllocateDevice
    integer ( KDI ), intent ( in ) :: &
      nCurrent, &
      nValues

    allocate ( SD % Current_ICL, SD % Current_ICR )
    call SD % Current_ICL % Initialize &
           ( [ nValues, nCurrent ], ClearOption = .true. )
    call SD % Current_ICR % Initialize &
           ( [ nValues, nCurrent ], ClearOption = .true. )
    
    if ( AllocateDevice ) then
      call SD % Current_ICL % AllocateDevice ( )
      call SD % Current_ICR % AllocateDevice ( )
    end if

  end subroutine Allocate_HLLC


  subroutine Deallocate ( SD )

    class ( StorageDivergenceForm ), intent ( inout ) :: &
      SD

    if ( allocated ( SD % GradientReconstructed ) ) &
      deallocate ( SD % GradientReconstructed )
    if ( allocated ( SD % DiffusionFactor_I ) ) &
      deallocate ( SD % DiffusionFactor_I )
    if ( allocated ( SD % SolverSpeeds_I ) ) &
      deallocate ( SD % SolverSpeeds_I )
    if ( allocated ( SD % Current_ICR ) ) &
      deallocate ( SD % Current_ICR )
    if ( allocated ( SD % Current_ICL ) ) &
      deallocate ( SD % Current_ICL )
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
