module SurfaceIntegral_Form

  use Basics
  use Manifolds
  use Geometries

  implicit none
  private

  type, public :: SurfaceIntegralForm
    integer ( KDI ) :: &
      IGNORABILITY, &
      nIntegrals
    character ( LDL ) :: &
      Name
    type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
      Integrand  !-- surface slab, coordinate basis
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
  end type SurfaceIntegralForm

    private :: &
      ComputeIntegral_SCG

contains


  subroutine Initialize ( SI, G, nIntegrals, NameOption, IgnorabilityOption )

    class ( SurfaceIntegralForm ), intent ( inout ) :: &
      SI
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nIntegrals
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension, etc.
      iI  !-- iIntegrand
    integer ( KDI ), dimension ( 3 ) :: &
      nS  !-- nSurface

    SI % IGNORABILITY  =  CONSOLE % INFO_2
    if ( present ( IgnorabilityOption ) ) &
      SI % IGNORABILITY  =  IgnorabilityOption

    SI % Name  =  'SurfaceIntegral'
    if ( present ( NameOption ) ) &
      SI % Name  =  NameOption

    call Show ( 'Initializing a SurfaceIntegral', SI % IGNORABILITY )
    call Show ( SI % Name, 'Name', SI % IGNORABILITY )

    SI % nIntegrals  =   nIntegrals
    SI % Geometry    =>  G

    allocate ( SI % Output ( nIntegrals ) )

    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    allocate ( SI % Integrand ( nIntegrals, 2 * C % nDimensions ) )
    associate ( I  =>  SI % Integrand )
    do iD  =  1,  C % nDimensions
      jD  =  mod ( iD, 3 ) + 1
      kD  =  mod ( jD, 3 ) + 1
      nS ( iD )  =  1
      nS ( jD )  =  C % nCellsBrick ( jD )
      nS ( kD )  =  C % nCellsBrick ( kD )
      do iI  =  1,  nIntegrals
        associate &
          ( I_I  =>  SI % Integrand ( iI, C % Connectivity % iaInner ( iD ) ), &
            I_O  =>  SI % Integrand ( iI, C % Connectivity % iaOuter ( iD ) ) )
          call I_I % Initialize ( nS, ClearOption = .true. )
          call I_O % Initialize ( nS, ClearOption = .true. )
        end associate !-- I_I, I_O
      end do !-- iI
    end do !-- iD
    end associate !-- I

    end associate !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'SurfaceIntegral_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Initialize


  subroutine Compute ( SI, ReduceOption )

    class ( SurfaceIntegralForm ), intent ( inout ) :: &
      SI
    logical ( KDL ), intent ( in ), optional :: &
      ReduceOption

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension
      iI  !-- iIntegral
    integer ( KDI ), dimension ( 3 ) :: &
      nS, &  !-- nSurface
      LB, UB  !-- LowerBound, UpperBound
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Area_I_D, &
      dA_I, dA_O
    logical ( KDR ) :: &
      Reduce
    type ( CollectiveOperation_R_Form ) :: &
      CO

    call Show ( 'Computing a SurfaceIntegral', SI % IGNORABILITY + 1 )
    call Show ( SI % Name, 'Name', SI % IGNORABILITY + 1 )

    Reduce = .true.
    if ( present ( ReduceOption ) ) &
      Reduce = ReduceOption

    SI % Output  =  0.0_KDR

    associate &
      ( G  =>  SI % Geometry, &
        I  =>  SI % Integrand )

    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      (  C   =>   A % Chart_GS, &
         GV  =>   G % Storage ( 1 ) % Value, &
        nI   =>  SI % nIntegrals )
    associate &
      ( Cy  =>  C % Connectivity )

    do iD  =  1, C % nDimensions

      jD  =  mod ( iD, 3 ) + 1
      kD  =  mod ( jD, 3 ) + 1
    
      call C % SetFieldPointer &
             ( GV ( :, G % AREA_I_D ( iD ) ), Area_I_D )

      nS ( iD )  =  1
      nS ( jD )  =  C % nCellsBrick ( jD )
      nS ( kD )  =  C % nCellsBrick ( kD )

      !-- Inner boundary
      if ( C % iaBrick ( iD )  ==  1 ) then

        LB = 1
        UB = nS

        dA_I  =>  Area_I_D ( LB ( 1 ) : UB ( 1 ), &
                             LB ( 2 ) : UB ( 2 ), &
                             LB ( 3 ) : UB ( 3 ) )

        do iI  =  1,  nI
          associate ( dIdA  =>  I ( iI, Cy % iaInner ( iD ) ) % Value )
          !-- Outward normal points left
!          MyIntegral ( iI ) = MyIntegral ( iI ) - sum ( dIdA * dA )
          call ComputeIntegral_SCG &
                 ( SI % Output ( iI ), dIdA, dA_I, Direction = -1.0_KDR )
          end associate !-- dIdA
        end do !-- iI

      end if !-- Inner boundary

      !-- Outer boundary
      if ( C % iaBrick ( iD )  ==  C % nBricks ( iD ) ) then

        LB = 1
        UB = nS
        LB ( iD )  =  LB ( iD )  +  C % nCellsBrick ( iD )
        UB ( iD )  =  UB ( iD )  +  C % nCellsBrick ( iD )

        dA_O  =>  Area_I_D ( LB ( 1 ) : UB ( 1 ), &
                             LB ( 2 ) : UB ( 2 ), &
                             LB ( 3 ) : UB ( 3 ) )

        do iI = 1, nI
          associate ( dIdA  =>  I ( iI, Cy % iaOuter ( iD ) ) % Value )
          !-- Outward normal points right
!          MyIntegral ( iI ) = MyIntegral ( iI ) + sum ( dIdA * dA )
          call ComputeIntegral_SCG &
                 ( SI % Output ( iI ), dIdA, dA_O, Direction = +1.0_KDR )
          end associate !-- dIdA
        end do !-- iI

      end if !-- Outer boundary

    end do !-- iD

    call Show ( SI % Output, 'MyIntegral', SI % IGNORABILITY + 1 )

    if ( C % Distributed .and. Reduce ) then
      call CO % Initialize &
             ( C % Communicator, nOutgoing = [ nI ], nIncoming = [ nI ] )
      CO % Outgoing % Value  =  SI % Output
      call CO % Reduce ( REDUCTION % SUM )
      SI % Output  =  CO % Incoming % Value
    end if !-- Reduce

    call Show ( SI % Output, 'Integral', SI % IGNORABILITY + 1 )

    end associate !-- Cy
    end associate !-- C, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'SurfaceIntegral_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- G, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( SI )

    type ( SurfaceIntegralForm ), intent ( inout ) :: &
      SI

    nullify ( SI % Geometry )

    if ( allocated ( SI % Integrand ) ) &
      deallocate ( SI % Integrand )
    if ( allocated ( SI % Output ) ) &
      deallocate ( SI % Output )

    call Show ( 'Finalizing a SurfaceIntegral', SI % IGNORABILITY )
    call Show ( SI % Name, 'Name', SI % IGNORABILITY )

  end subroutine Finalize


  subroutine ComputeIntegral_SCG ( I, dIdA, dA, Direction )

    real ( KDR ), intent ( inout ) :: &
      I
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      dIdA, &
      dA
    real ( KDR ), intent ( in ) :: &
      Direction

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV

    nV = shape ( dIdA )

    !$OMP parallel do reduction ( + : I ) collapse ( 3 )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          I  =  I  +  Direction * dIdA ( iV, jV, kV ) * dA ( iV, jV, kV )
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ComputeIntegral_SCG


end module SurfaceIntegral_Form
