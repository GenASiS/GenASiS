module VolumeIntegral_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none
  private

  type, public :: VolumeIntegralForm
    real ( KDR ), dimension ( : ), allocatable :: &
      Output
  contains
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type VolumeIntegralForm

    private :: &
      ComputeIntegral_CGS


contains


  subroutine Compute ( VI, I, G, ReduceOption, IgnorabilityOption )

    class ( VolumeIntegralForm ), intent ( inout ) :: &
      VI
    class ( FieldSetForm ), intent ( in ) :: &
      I  !-- Integrand
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    logical ( KDL ), intent ( in ), optional :: &
      ReduceOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iI, &  !-- iIntegral
      iF, &  !-- iF
      Ignorability
    real ( KDR ), dimension ( : ), allocatable :: &
      MyIntegral
    logical ( KDR ) :: &
      Reduce
    type ( CollectiveOperation_R_Form ) :: &
      CO

    Reduce = .true.
    if ( present ( ReduceOption ) ) &
      Reduce = ReduceOption

    Ignorability = CONSOLE % INFO_5
    if ( present ( IgnorabilityOption ) ) &
      Ignorability = IgnorabilityOption

    call Show ( 'Computing an integral', Ignorability )
    call Show ( I % Name, 'Integrand', Ignorability )
    call Show ( I % Atlas % Name, 'Atlas', Ignorability )

    select type ( A  =>  I % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      (  C   =>  A % Chart_GS, &
         IV  =>  I % Storage ( 1 ) % Value, &
         GV  =>  G % Storage ( 1 ) % Value, &
        nI   =>  I % nFields )

    allocate ( VI % Output ( nI ) )
    allocate ( MyIntegral ( nI ) )

    if ( C % Distributed .and. Reduce ) then
      call CO % Initialize &
             ( C % Communicator, &
               nOutgoing = [ nI ], nIncoming = [ nI ] )
    end if

    do iI = 1, nI
      iF  =  I % iaSelected ( iI )
      call ComputeIntegral_CGS &
             ( C % ProperCell, IV ( :, iF ), GV ( :, G % VOLUME ), &
               MyIntegral ( iI ) )
    end do !-- iI
    call Show ( MyIntegral, 'MyIntegral', Ignorability )

    if ( C % Distributed .and. Reduce ) then
      CO % Outgoing % Value  =  MyIntegral
      call CO % Reduce ( REDUCTION % SUM )
      VI % Output  =  CO % Incoming % Value
    else
      VI % Output  =  MyIntegral
    end if

    call Show ( VI % Output, 'Integral', Ignorability )

    end associate !-- C, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'VolumeIntegral_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Compute


  subroutine Finalize ( VI )

    type ( VolumeIntegralForm ), intent ( inout ) :: &
      VI

    if ( allocated ( VI % Output ) ) &
      deallocate ( VI % Output )

  end subroutine Finalize


  subroutine ComputeIntegral_CGS ( ProperCell, dIdV, dV, I )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      dIdV, &
      dV
    real ( KDR ), intent ( out ) :: &
      I

    integer ( KDI ) :: &
      iV, &
      nV

    nV  =  size ( dIdV )
    
    I  =  0.0_KDR

    !$OMP parallel do reduction ( + : I )
    do iV  =  1, nV
      if ( ProperCell ( iV ) ) &
        I  =  I  +  dIdV ( iV ) * dV ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeIntegral_CGS


end module VolumeIntegral_Form
