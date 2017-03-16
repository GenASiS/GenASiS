!-- VolumeIntegral computes a volume integral on a chart.

module VolumeIntegral_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: VolumeIntegralForm
  contains
    procedure, private, nopass :: &
      Compute_CSL
    generic :: &
      Compute => Compute_CSL
  end type VolumeIntegralForm

    private :: &
      ComputeVolume, &
      ComputeIntegral_CSL

contains


  subroutine Compute_CSL ( CSL, Integrand, Integral, IgnorabilityOption )

    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    type ( Real_1D_Form ), dimension ( : ), intent ( in ) :: &
      Integrand
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Integral
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iI, &  !-- iIntegral
      Ignorability
    real ( KDR ), dimension ( size ( Integral ) ) :: &
      MyIntegral
    real ( KDR ), dimension ( : ), allocatable :: &
      dVolume
    type ( CollectiveOperation_R_Form ) :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G

    Ignorability = CONSOLE % INFO_5
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    associate ( nI => size ( Integral ) )

    G => CSL % Geometry ( )

    allocate ( dVolume ( G % nValues ) )
    call ComputeVolume &
           ( G % Value ( :, G % VOLUME_JACOBIAN ), &
             G % Value ( :, G % WIDTH ( 1 ) ), &
             G % Value ( :, G % WIDTH ( 2 ) ), &
             G % Value ( :, G % WIDTH ( 3 ) ), &
             CSL % nDimensions, dVolume )

    associate ( A => CSL % Atlas )

    if ( A % IsDistributed ) then
      call CO % Initialize &
             ( CSL % Atlas % Communicator, &
               nOutgoing = [ nI ], nIncoming = [ nI ] )
    end if

    do iI = 1, nI
      call ComputeIntegral_CSL &
             ( CSL % IsProperCell, Integrand ( iI ) % Value, dVolume, &
               MyIntegral ( iI ) )
    end do !-- iI
    call Show ( 'Local contribution to VolumeIntegral', Ignorability )
    call Show ( MyIntegral, 'MyIntegral', Ignorability )

    if ( A % IsDistributed ) then
      CO % Outgoing % Value = MyIntegral
      call CO % Reduce ( REDUCTION % SUM )
      Integral = CO % Incoming % Value
    else
      Integral = MyIntegral
    end if

    end associate !-- A
    end associate !-- nI, etc.
    nullify ( G )

  end subroutine Compute_CSL


  subroutine ComputeVolume ( VJ, dX_1, dX_2, dX_3, nDimensions, dV )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      VJ, &
      dX_1, dX_2, dX_3
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      dV

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( VJ )

    select case ( nDimensions )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        dV ( iV )  =  VJ ( iV ) * dX_1 ( iV )
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        dV ( iV )  =  VJ ( iV ) * dX_1 ( iV ) * dX_2 ( iV )
      end do
      !$OMP end parallel do
    case ( 3 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        dV ( iV )  =  VJ ( iV ) * dX_1 ( iV ) * dX_2 ( iV ) * dX_3 ( iV )
      end do
      !$OMP end parallel do
    end select !-- nDimensions

  end subroutine ComputeVolume


  subroutine ComputeIntegral_CSL ( IsProperCell, dIdV, dV, I )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      dIdV, &
      dV
    real ( KDR ), intent ( out ) :: &
      I

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( dIdV )
    
    I = 0.0_KDR

    !$OMP parallel do private ( iV ) reduction ( + : I )
    do iV = 1, nV
      if ( IsProperCell ( iV ) ) &
        I  =  I  +  dIdV ( iV ) * dV ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeIntegral_CSL


end module VolumeIntegral_Form
