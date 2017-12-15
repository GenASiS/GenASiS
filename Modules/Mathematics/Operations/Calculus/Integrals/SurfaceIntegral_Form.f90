!-- SurfaceIntegral computes a surface integral on a chart.

module SurfaceIntegral_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: SurfaceIntegralForm
  contains
    procedure, private, nopass :: &
      Compute_CSL
    generic :: &
      Compute => Compute_CSL
  end type SurfaceIntegralForm

    private :: &
      ComputeIntegral_CSL

contains


  subroutine Compute_CSL &
               ( CSL, Integrand, Integral, ReduceOption, IgnorabilityOption )

    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      Integrand  !-- surface slab, coordinate basis
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Integral
    logical ( KDL ), intent ( in ), optional :: &
      ReduceOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension
      iI, &  !-- iIntegral
      Ignorability
    integer ( KDI ), dimension ( 3 ) :: &
      nS, &  !-- nSurface
      LB, UB  !-- LowerBound, UpperBound
    real ( KDR ), dimension ( size ( Integral ) ) :: &
      MyIntegral
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      AreaInner_D, &
      dA_I, dA_O
    logical ( KDR ) :: &
      Reduce
    type ( CollectiveOperation_R_Form ) :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G

    Reduce = .true.
    if ( present ( ReduceOption ) ) Reduce = ReduceOption

    Ignorability = CONSOLE % INFO_5
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    associate &
      ( nI => size ( Integral ), &
        C  => CSL % Atlas % Connectivity )

    G => CSL % Geometry ( )

    call Clear ( MyIntegral )

    do iD = 1, CSL % nDimensions

      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1
    
      call CSL % SetVariablePointer &
             ( G % Value ( :, G % AREA_INNER_D ( iD ) ), AreaInner_D )

      nS ( iD ) = 1
      nS ( jD ) = CSL % nCellsBrick ( jD )
      nS ( kD ) = CSL % nCellsBrick ( kD )

      !-- Inner boundary
      if ( CSL % iaBrick ( iD ) == 1 ) then

        LB = 1
        UB = nS

        dA_I => AreaInner_D ( LB ( 1 ) : UB ( 1 ), &
                              LB ( 2 ) : UB ( 2 ), &
                              LB ( 3 ) : UB ( 3 ) )

        do iI = 1, nI
          associate ( dIdA => Integrand ( iI, C % iaInner ( iD ) ) % Value )
          !-- Outward normal points left
!          MyIntegral ( iI ) = MyIntegral ( iI ) - sum ( dIdA * dA )
          call ComputeIntegral_CSL &
                 ( MyIntegral ( iI ), dIdA, dA_I, Direction = -1.0_KDR )
          end associate !-- dIdA
        end do !-- iI

      end if !-- iaBrick ( iD ) == 1

      !-- Outer boundary
      if ( CSL % iaBrick ( iD ) == CSL % nBricks ( iD ) ) then

        LB = 1
        UB = nS
        LB ( iD ) = LB ( iD ) + CSL % nCellsBrick ( iD )
        UB ( iD ) = UB ( iD ) + CSL % nCellsBrick ( iD )

        dA_O => AreaInner_D ( LB ( 1 ) : UB ( 1 ), &
                              LB ( 2 ) : UB ( 2 ), &
                              LB ( 3 ) : UB ( 3 ) )

        do iI = 1, nI
          associate ( dIdA => Integrand ( iI, C % iaOuter ( iD ) ) % Value )
          !-- Outward normal points right
!          MyIntegral ( iI ) = MyIntegral ( iI ) + sum ( dIdA * dA )
          call ComputeIntegral_CSL &
                 ( MyIntegral ( iI ), dIdA, dA_O, Direction = +1.0_KDR )
          end associate !-- dIdA
        end do !-- iI

      end if !-- iaBrick ( iD ) == nBricks ( iD )

    end do !-- iD

    call Show ( 'Local contribution to SurfaceIntegral', Ignorability )
    call Show ( MyIntegral, 'MyIntegral', Ignorability )

    if ( Reduce ) then
      call CO % Initialize &
             ( CSL % Atlas % Communicator, &
               nOutgoing = [ nI ], nIncoming = [ nI ] )
      CO % Outgoing % Value = MyIntegral
      call CO % Reduce ( REDUCTION % SUM )
      Integral = CO % Incoming % Value
    else !-- don't reduce
      Integral = MyIntegral
    end if !-- Reduce

    nullify ( AreaInner_D, dA_I, dA_O, G )
    end associate !-- nI, etc.

  end subroutine Compute_CSL


  subroutine ComputeIntegral_CSL ( I, dIdA, dA, Direction )

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

    !$OMP parallel do private ( iV, jV, kV ) reduction ( + : I ) collapse ( 3 )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          I  =  I  +  Direction * dIdA ( iV, jV, kV ) * dA ( iV, jV, kV )
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ComputeIntegral_CSL


end module SurfaceIntegral_Form
