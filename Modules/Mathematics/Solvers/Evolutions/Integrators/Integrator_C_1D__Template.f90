!-- Integrator_C_1D is a template for time evolution of multiple conserved 
!   currents.

module Integrator_C_1D__Template

  !-- Integrator_Current_1D__Template

  use Basics
  use Manifolds
  use Fields
  use Integrator_C__Template

  implicit none
  private

  type, public, extends ( Integrator_C_Template ), abstract :: &
    Integrator_C_1D_Template
      integer ( KDI ) :: &
        N_CURRENTS
      type ( Current_ASC_ElementForm ), dimension ( : ), allocatable :: &
        Current_ASC_1D
  contains
    procedure, public, pass :: &
      FinalizeTemplate_C_1D
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &
      ComputeTimeStepLocal
  end type Integrator_C_1D_Template

contains


  subroutine FinalizeTemplate_C_1D ( I )

    class ( Integrator_C_1D_Template ), intent ( inout ) :: &
      I

   if ( allocated ( I % Current_ASC_1D ) ) &
     deallocate ( I % Current_ASC_1D )

    call I % FinalizeTemplate_C ( )

  end subroutine FinalizeTemplate_C_1D


  subroutine ComputeTally ( I, ComputeChangeOption )

    class ( Integrator_C_1D_Template ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption

    integer ( KDI ) :: &
      iC  !-- iCurrent

    if ( .not. allocated ( I % Current_ASC_1D ) ) &
      return

    associate ( Timer => PROGRAM_HEADER % Timer ( I % iTimerComputeTally ) )
    call Timer % Start ( )

    do iC = 1, I % N_CURRENTS
      associate ( CA => I % Current_ASC_1D ( iC ) % Element )
      call CA % ComputeTally ( ComputeChangeOption = ComputeChangeOption )
      end associate !-- CA
    end do !-- iC

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeTally


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( Integrator_C_1D_Template ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    integer ( KDI ) :: &
      iC  !-- iCurrent
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( CurrentTemplate ), pointer :: &
      C

    if ( .not. allocated ( I % Current_ASC_1D ) ) &
      return

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( CSL => PS % Chart )
    class is ( Chart_SL_Template )

    G => CSL % Geometry ( )

    do iC = 1, I % N_CURRENTS
      associate ( CA => I % Current_ASC_1D ( iC ) % Element )
      C => CA % Current ( )

      call I % ComputeTimeStepKernel_CSL &
             ( CSL % IsProperCell, &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
               G % Value ( :, G % WIDTH ( 1 ) ), &
               G % Value ( :, G % WIDTH ( 2 ) ), & 
               G % Value ( :, G % WIDTH ( 3 ) ), &
               CSL % nDimensions, TimeStepCandidate ( iC ) )

      end associate !-- CA
    end do !-- iC

    end select !-- CSL
    end select !-- PS

    TimeStepCandidate = I % CourantFactor * TimeStepCandidate

    nullify ( C, G )

  end subroutine ComputeTimeStepLocal


end module Integrator_C_1D__Template
