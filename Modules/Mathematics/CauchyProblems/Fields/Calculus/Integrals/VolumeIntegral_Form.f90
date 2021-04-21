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
    procedure, private, pass :: &
      Compute_ASG
    generic, public :: &
      Compute => Compute_ASG
    final :: &
      Finalize
  end type VolumeIntegralForm

    private :: &
      ComputeIntegral_GS


contains


  subroutine Compute_ASG ( VI, IA, GA, ReduceOption, IgnorabilityOption )

    class ( VolumeIntegralForm ), intent ( inout ) :: &
      VI
    class ( FieldSet_ASG_Form ), intent ( in ) :: &
      IA  !-- Integrand
    class ( Geometry_F_ASG_Form ), intent ( in ) :: &
      GA
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
    call Show ( IA % Name, 'Integrand', Ignorability )
    call Show ( IA % Atlas % Name, 'Atlas', Ignorability )

    !-- Integrand

    associate ( IG  =>  IA % FieldSet_G )
    associate ( IV  =>  IG % Storage % Value )

    !-- Geometry

    associate ( GG  =>  GA % Geometry_G )
    associate ( GV  =>  GA % FieldSet_G % Storage % Value )

    !-- Grid

    select type ( G  =>  IG % Chart ) 
    type is ( Grid_S_Form )

    !-- Integrals

    associate ( nI  =>  IG % nFields )

    allocate ( VI % Output ( nI ) )
    allocate ( MyIntegral ( nI ) )

    if ( G % Distributed .and. Reduce ) then
      call CO % Initialize &
             ( G % Communicator, &
               nOutgoing = [ nI ], nIncoming = [ nI ] )
    end if


    do iI = 1, nI
      iF  =  IG % iaSelected ( iI )
      call ComputeIntegral_GS &
             ( G % ProperCell, IV ( :, iF ), GV ( :, GG % VOLUME ), &
               MyIntegral ( iI ) )
    end do !-- iI
    call Show ( MyIntegral, 'MyIntegral', Ignorability )

    if ( G % Distributed .and. Reduce ) then
      CO % Outgoing % Value  =  MyIntegral
      call CO % Reduce ( REDUCTION % SUM )
      VI % Output  =  CO % Incoming % Value
    else
      VI % Output  =  MyIntegral
    end if

    call Show ( VI % Output, 'Integral', Ignorability )

    end associate !-- nI

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'VolumeIntegral_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_ASG', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G

    end associate !-- GV
    end associate !-- GG

    end associate !-- FSV
    end associate !-- FSG

  end subroutine Compute_ASG


  subroutine Finalize ( VI )

    type ( VolumeIntegralForm ), intent ( inout ) :: &
      VI

    if ( allocated ( VI % Output ) ) &
      deallocate ( VI % Output )

  end subroutine Finalize


  subroutine ComputeIntegral_GS ( ProperCell, dIdV, dV, I )

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

  end subroutine ComputeIntegral_GS


end module VolumeIntegral_Form
