!-- SphericalAverage computes a spherical average on a chart.

module SphericalAverage_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: SphericalAverageForm
  contains
    procedure, private, nopass :: &
      Compute_CSL
    generic :: &
      Compute => Compute_CSL
  end type SphericalAverageForm

    private :: &
      ComputeSolidAngleKernel, &
      ComputeAverageKernel

contains


  subroutine Compute_CSL ( Source, Target, CSL, CSL_SA, IgnorabilityOption )

    type ( StorageForm ), intent ( inout ) :: &
      Source, &
      Target
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL, &
      CSL_SA
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      oB, &  !-- oBuffer
      iA, &  !-- iAverage
      Ignorability
    integer ( KDI ), dimension ( : ), allocatable :: &
      iRadius
    real ( KDR ), dimension ( : ), allocatable :: &
      SolidAngle
    ! real ( KDR ), dimension ( size ( Integral ) ) :: &
    !   MyIntegral
    type ( CollectiveOperation_R_Form ) :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G

    Ignorability = CONSOLE % INFO_5
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    G => CSL % Geometry ( )

    allocate ( iRadius ( G % nValues ) )
    allocate ( SolidAngle ( G % nValues ) )

    call ComputeSolidAngleKernel &
           ( SolidAngle, iRadius, &
             G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), &
             G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
             G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), &
             G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
             G % Value ( :, G % CENTER_U ( 1 ) ), &
             G % Value ( :, G % CENTER_U ( 2 ) ), &
             CSL_SA % Edge ( 1 ) % Value, CSL % nDimensions )

    call CO % Initialize &
           ( CSL % Atlas % Communicator, &
             nOutgoing = [ Target % nValues  *  Target % nVariables ], &
             nIncoming = [ Target % nValues  *  Target % nVariables ], &
             RootOption = CONSOLE % DisplayRank )

    oB = 0
    do iA = 1, Target % nVariables
      associate &
        ( MyAverage => CO % Outgoing &
                          % Value ( oB + 1 : oB + Target % nValues ), &
          Variable  => Source % Value ( :, Source % iaSelected ( iA ) ) )
      call Clear ( MyAverage )
      call ComputeAverageKernel &
             ( MyAverage, CSL % IsProperCell, Variable, SolidAngle, iRadius )
      oB  =  oB  +  Target % nValues
      end associate !-- MyAverage, etc. 
    end do !-- iA

    call CO % Reduce ( REDUCTION % SUM )

    if ( CSL % Atlas % Communicator % Rank == CONSOLE % DisplayRank ) then
      oB = 0
      do iA = 1, Target % nVariables
        associate &
          ( Average => CO % Incoming &
                          % Value ( oB + 1 : oB + Target % nValues ), &
            Variable_SA => Target % Value ( :, Target % iaSelected ( iA ) ) )
        call Copy ( Average, Variable_SA )
        oB  =  oB  +  Target % nValues
        end associate !-- MyIntegral, etc. 
      end do !-- iA
    end if

    nullify ( G )

  end subroutine Compute_CSL


  subroutine ComputeSolidAngleKernel &
               ( dOmega, iR, W_L_2, W_L_3, W_R_2, W_R_3, R_C, Th_C, &
                 EdgeValue, nDimensions )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dOmega
    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      iR
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      W_L_2, W_L_3, &
      W_R_2, W_R_3, &
      R_C, &
      Th_C, &
      EdgeValue
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      TwoPi, FourPi, &
      Th_I, Th_O, &
      dPh
    
    TwoPi   =  2.0_KDR  *  CONSTANT % PI
    FourPi  =  4.0_KDR  *  CONSTANT % PI

    !$OMP parallel do private ( iV, Th_I, Th_O, dPh )
    do iV = 1, size ( dOmega )
      
      Th_I  =  Th_C ( iV )  -  W_L_2 ( iV )
      Th_O  =  Th_C ( iV )  +  W_R_2 ( iV )

      dPh  =  W_L_3 ( iV )  +  W_R_3 ( iV ) 

      select case ( nDimensions )
      case ( 1 )
        dOmega ( iV )  =  FourPi
      case ( 2 )
        dOmega ( iV )  =  TwoPi  *  ( cos ( Th_I )  -  cos ( Th_O ) )
      case ( 3 ) 
        dOmega ( iV )  =  dPh  *  ( cos ( Th_I )  -  cos ( Th_O ) )
      end select !-- nDimensions

      call Search ( EdgeValue, R_C ( iV ), iR ( iV ) )

    end do
    !$OMP end parallel do

  end subroutine ComputeSolidAngleKernel


  subroutine ComputeAverageKernel ( MyA, IsProperCell, V, dOmega, iR )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      MyA
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V, &
      dOmega
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iR

    integer ( KDI ) :: &
      iV

    !$OMP parallel do private ( iV )
    do iV = 1, size ( V )
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      MyA ( iR ( iV ) )  =  MyA ( iR ( iV ) )  +  V ( iV ) * dOmega ( iV )
    end do
    !$OMP end parallel do
 
    MyA  =  MyA  /  ( 4.0_KDR  *  CONSTANT % PI )

  end subroutine ComputeAverageKernel


end module SphericalAverage_Form
