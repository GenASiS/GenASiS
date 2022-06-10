module VolumeIntegral_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none
  private

  type, public :: VolumeIntegralForm
    integer ( KDI ) :: &
      IGNORABILITY, &
      nIntegrals
    character ( LDL ) :: &
      Name
    type ( FieldSetForm ), allocatable :: &
      Integrand
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
    real ( KDR ), dimension ( : ), allocatable :: &
      Output
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type VolumeIntegralForm

    private :: &
      ComputeIntegral_SCG


contains


  subroutine Initialize ( VI, G, nIntegrals, NameOption, IgnorabilityOption )

    class ( VolumeIntegralForm ), intent ( inout ) :: &
      VI
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nIntegrals
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    VI % IGNORABILITY  =  CONSOLE % INFO_2
    if ( present ( IgnorabilityOption ) ) &
      VI % IGNORABILITY  =  IgnorabilityOption

    VI % Name  =  'VolumeIntegral'
    if ( present ( NameOption ) ) &
      VI % Name  =  NameOption

    call Show ( 'Initializing a VolumeIntegral', VI % IGNORABILITY )
    call Show ( VI % Name, 'Name', VI % IGNORABILITY )

    VI % nIntegrals  =   nIntegrals
    VI % Geometry    =>  G

    allocate ( VI % Output ( nIntegrals ) )

    allocate ( VI % Integrand )
    associate ( I  =>  VI % Integrand )
    call I % Initialize &
           ( G % Atlas, NameOption = 'Integrand_' // trim ( VI % Name ), &
             nFieldsOption = nIntegrals )
    end associate !-- I

  end subroutine Initialize


  subroutine Compute ( VI, ReduceOption )

    class ( VolumeIntegralForm ), intent ( inout ) :: &
      VI
    logical ( KDL ), intent ( in ), optional :: &
      ReduceOption

    logical ( KDR ) :: &
      Reduce
    type ( CollectiveOperation_R_Form ) :: &
      CO

    call Show ( 'Computing a VolumeIntegral', VI % IGNORABILITY + 1 )
    call Show ( VI % Name, 'Name', VI % IGNORABILITY + 1 )

    Reduce = .true.
    if ( present ( ReduceOption ) ) &
      Reduce = ReduceOption

    associate &
      ( G  =>  VI % Geometry, &
        I  =>  VI % Integrand )

    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      (  C   =>   A % Chart_GS, &
         IV  =>   I % Storage ( 1 ) % Value, &
         GV  =>   G % Storage ( 1 ) % Value, &
        nI   =>  VI % nIntegrals )

    call ComputeIntegral_SCG &
           ( C % ProperCell, IV, GV ( :, G % VOLUME ), VI % Output )
    call Show ( VI % Output, 'MyIntegral', VI % IGNORABILITY + 1 )

    if ( C % Distributed .and. Reduce ) then
      call CO % Initialize &
             ( C % Communicator, nOutgoing = [ nI ], nIncoming = [ nI ] )
      CO % Outgoing % Value  =  VI % Output
      call CO % Reduce ( REDUCTION % SUM )
      VI % Output  =  CO % Incoming % Value
    end if

    call Show ( VI % Output, 'Integral', VI % IGNORABILITY + 1 )

    end associate !-- C, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'VolumeIntegral_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- G, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( VI )

    type ( VolumeIntegralForm ), intent ( inout ) :: &
      VI

    nullify ( VI % Geometry )

    if ( allocated ( VI % Integrand ) ) &
      deallocate ( VI % Integrand )
    if ( allocated ( VI % Output ) ) &
      deallocate ( VI % Output )

    call Show ( 'Finalizing a VolumeIntegral', VI % IGNORABILITY )
    call Show ( VI % Name, 'Name', VI % IGNORABILITY )

  end subroutine Finalize


  subroutine ComputeIntegral_SCG ( ProperCell, dIdV, dV, I )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      dIdV
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      dV
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I

    integer ( KDI ) :: &
      iV, &
      nV

    nV  =  size ( ProperCell )
    
    I  =  0.0_KDR

    !$OMP parallel do reduction ( + : I )
    do iV  =  1, nV
      if ( ProperCell ( iV ) ) &
        I  =  I  +  dIdV ( iV, : ) * dV ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeIntegral_SCG


end module VolumeIntegral_Form
