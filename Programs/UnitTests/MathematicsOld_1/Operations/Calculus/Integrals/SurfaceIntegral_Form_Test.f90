module SurfaceIntegral_Form_Test__Form

  use Basics
  use Manifolds
  use SurfaceIntegral_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'SurfaceIntegral_Form_Test'

  type, public :: SurfaceIntegral_Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( Atlas_SC_CC_Form ), allocatable :: &
      Atlas
    type ( SurfaceIntegralForm ), allocatable :: &
      SurfaceIntegral
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type SurfaceIntegral_Form_Test_Form

contains


  subroutine Initialize ( SIFT, Name )

    class ( SurfaceIntegral_Form_Test_Form ), intent ( inout ), target :: &
      SIFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iD, jD, kD  !-- iDimension, etc.
    integer ( KDI ), dimension ( 3 ) :: &
      nS  !-- nSurface
    real ( KDR ), dimension ( 1 ) :: &
      Integral
    type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
      Integrand

    !-- Stream

    allocate ( SIFT % GridImageStream )
    associate ( GIS => SIFT % GridImageStream )
    call GIS % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )
    
    !-- Atlas, etc.

    allocate ( SIFT % Atlas )
    associate ( A => SIFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( )  
    call A % SetGeometry ( )

    associate ( Cnnct => A % Connectivity )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
       
    !-- Integrand

    allocate ( Integrand ( 1, 2 * C % nDimensions ) )

    do iD = 1, C % nDimensions
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1
      nS ( iD ) = 1
      nS ( jD ) = C % nCellsBrick ( jD )
      nS ( kD ) = C % nCellsBrick ( kD )
      associate &
        ( I_I => Integrand ( 1, Cnnct % iaInner ( iD ) ), &
          I_O => Integrand ( 1, Cnnct % iaOuter ( iD ) ) )
        call I_I % Initialize ( nS )
        call I_O % Initialize ( nS )
        select case ( iD )
        case ( 1, 2 ) !-- outward normal
          I_I % Value  =  -1.0_KDR
          I_O % Value  =  +1.0_KDR
        case ( 3 )  !-- must be continuous at periodic boundary
          I_I % Value  =  +1.0_KDR
          I_O % Value  =  +1.0_KDR
        end select !-- iD
      end associate !-- I_I, I_O
    end do !-- iD

    !-- SurfaceIntegral

    allocate ( SIFT % SurfaceIntegral )
    associate ( SI => SIFT % SurfaceIntegral )

    call SI % Compute ( C, Integrand, Integral )
    call Show ( Integral ( 1 ), '*** SurfaceIntegral', &
                nLeadingLinesOption = 1 )
    call Show ( 4.0_KDR  *  CONSTANT % PI  *  C % MaxCoordinate ( 1 ) ** 2, &
                '*** Expected', nTrailingLinesOption = 1 )

    !-- Write

    call A % OpenStream ( GIS, '1', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    !-- Cleanup

    end associate !-- SI
    end select !-- C
    end associate !-- Cnnct
    end associate !-- A
    end associate !-- GIS
    
  end subroutine Initialize


  subroutine Finalize ( SIFT )

    type ( SurfaceIntegral_Form_Test_Form ), intent ( inout ) :: &
      SIFT

    deallocate ( SIFT % SurfaceIntegral )
    deallocate ( SIFT % Atlas )
    deallocate ( SIFT % GridImageStream )

  end subroutine Finalize


end module SurfaceIntegral_Form_Test__Form



program SurfaceIntegral_Form_Test

  use Basics
  use SurfaceIntegral_Form_Test__Form

  implicit none
  
  type ( SurfaceIntegral_Form_Test_Form ), allocatable :: &
    SIFT
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )
    
  allocate ( SIFT )
  call SIFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( SIFT )

  deallocate ( PROGRAM_HEADER )

end program SurfaceIntegral_Form_Test
