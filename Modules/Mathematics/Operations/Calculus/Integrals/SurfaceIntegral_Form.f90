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
      ComputeArea_CSL, &
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
      oS, &  !-- oSurface
      nS     !-- nSurface
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
    
!      nS ( iD ) = 1
!      nS ( jD ) = CSL % nCellsBrick ( jD )
!      nS ( kD ) = CSL % nCellsBrick ( kD )
!      allocate ( dA ( nS ( 1 ), nS ( 2 ), nS ( 3 ) ) )
!      call Clear ( dA )
      
!      !-- assume VolumeJacobian already included in integrand
!      oS = CSL % nGhostLayers
!      call CSL % SetVariablePointer &
!             ( G % Value ( :, G % WIDTH ( jD ) ), dX_J )
!      call CSL % SetVariablePointer &
!             ( G % Value ( :, G % WIDTH ( kD ) ), dX_K )
!      call ComputeArea_CSL &
!             ( dX_J, dX_K, dA, nS, oS, CSL % nDimensions, jD, kD )

      call CSL % SetVariablePointer &
             ( G % Value ( :, G % AREA_INNER_D ( iD ) ), AreaInner_D )

      nS ( iD ) = 1
      nS ( jD ) = CSL % nCellsBrick ( jD )
      nS ( kD ) = CSL % nCellsBrick ( kD )

      !-- Inner boundary
      if ( CSL % iaBrick ( iD ) == 1 ) then

        oS = CSL % nGhostLayers

        dA_I => AreaInner_D ( oS ( iD ) + 1 : oS ( iD ) + nS ( iD ), &
                              oS ( jD ) + 1 : oS ( jD ) + nS ( jD ), &
                              oS ( kD ) + 1 : oS ( kD ) + nS ( kD ) )

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

        oS         =  CSL % nGhostLayers
        oS ( iD )  =  oS ( iD ) + CSL % nCellsBrick ( iD )

        dA_O => AreaInner_D ( oS ( iD ) + 1 : oS ( iD ) + nS ( iD ), &
                              oS ( jD ) + 1 : oS ( jD ) + nS ( jD ), &
                              oS ( kD ) + 1 : oS ( kD ) + nS ( kD ) )

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


  subroutine ComputeArea_CSL ( dX_J, dX_K, dA, nS, oS, nDimensions, jD, kD )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      dX_J, &
      dX_K
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      dA
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nS, &
      oS
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      jD, kD

    integer ( KDI ) :: &
      iV, jV, kV

    if ( jD <= nDimensions .and. kD <= nDimensions ) then
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nS ( 3 )
        do jV = 1, nS ( 2 )
          do iV = 1, nS ( 1 )
            dA ( iV, jV, kV ) &
              =    dX_J ( oS ( 1 ) + iV, oS ( 2 ) + jV, oS ( 3 ) + kV ) &
                 * dX_K ( oS ( 1 ) + iV, oS ( 2 ) + jV, oS ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
    else if ( jD <= nDimensions ) then
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nS ( 3 )
        do jV = 1, nS ( 2 )
          do iV = 1, nS ( 1 )
            dA ( iV, jV, kV ) &
              = dX_J ( oS ( 1 ) + iV, oS ( 2 ) + jV, oS ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
    else if ( kD <= nDimensions ) then
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nS ( 3 )
        do jV = 1, nS ( 2 )
          do iV = 1, nS ( 1 )
            dA ( iV, jV, kV ) &
              = dX_K ( oS ( 1 ) + iV, oS ( 2 ) + jV, oS ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
    else
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nS ( 3 )
        do jV = 1, nS ( 2 )
          do iV = 1, nS ( 1 )
            dA ( iV, jV, kV ) = 1.0_KDR
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
    end if

  end subroutine ComputeArea_CSL


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
