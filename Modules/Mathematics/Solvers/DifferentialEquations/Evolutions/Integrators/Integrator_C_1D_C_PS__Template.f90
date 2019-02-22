!-- Integrator_C_1D_MS_C_PS is a template for time evolution of multiple
!   similar conserved currents, and an additional conserved
!   current on position space.

module Integrator_C_1D_C_PS__Template

  !-- Integrator_Current_1D_Current_PositionSpace__Template

  use Basics
  use Manifolds
  use Fields
  use Steps
  use Integrator_C_PS__Template

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    Integrator_C_1D_C_PS_Template
      integer ( KDI ) :: &
        N_CURRENTS_1D = 0
      class ( Step_RK_C_ASC_Template ), allocatable :: &
        Step_1D
  contains
    procedure ( I_1D ), private, pass, deferred :: &
      Initialize_1D
    generic, public :: &
      Initialize => Initialize_1D
    procedure, public, pass :: &  !-- 1
      InitializeTemplate_C_1D_C_PS
    procedure, public, pass :: &  !-- 1
      FinalizeTemplate_C_1D_C_PS
    procedure, private, pass :: &  !-- 2
      ComputeCycle
    procedure, private, pass :: &  !-- 3
      InitializeStepTimers
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &
      ComputeCycle_C_1D_C_ASC
    procedure ( CT_1D ), private, pass, deferred :: &
      ComputeTally_1D
    procedure, private, pass :: &
      PrepareStep_1D
    procedure, private, pass :: &
      PrepareStep
    procedure, private, pass :: &
      ComputeTimeStepLocal
    procedure, public, pass :: &
      ComputeTimeStepLocalTemplate
    procedure ( CAP ), private, pass, deferred :: &
      Current_ASC_Pointer   
  end type Integrator_C_1D_C_PS_Template

  abstract interface

    subroutine I_1D ( I, Name, TimeUnitOption, FinishTimeOption, &
                      CourantFactorOption, nWriteOption )
      use Basics
      import Integrator_C_1D_C_PS_Template
      class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
        I
      character ( * ), intent ( in )  :: &
        Name
      type ( MeasuredValueForm ), intent ( in ), optional :: &
        TimeUnitOption
      real ( KDR ), intent ( in ), optional :: &
        FinishTimeOption, &
        CourantFactorOption
      integer ( KDI ), intent ( in ), optional :: &
        nWriteOption
    end subroutine I_1D

    subroutine CT_1D ( I, ComputeChangeOption, IgnorabilityOption )
      use Basics
      import Integrator_C_1D_C_PS_Template
      class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
        I
      logical ( KDL ), intent ( in ), optional :: &
        ComputeChangeOption
      integer ( KDI ), intent ( in ), optional :: &
        IgnorabilityOption
    end subroutine CT_1D

    function CAP ( I, iC ) result ( CA )
      use Basics
      use Fields
      import Integrator_C_1D_C_PS_Template
      class ( Integrator_C_1D_C_PS_Template ), intent ( inout ), target :: &
        I
      integer ( KDI ), intent ( in ) :: &
        iC
      class ( Current_ASC_Template ), pointer :: &
        CA
    end function CAP

  end interface

contains


  subroutine InitializeTemplate_C_1D_C_PS &
               ( I, Name, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nWriteOption )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I
    character ( * ), intent ( in )  :: &
      Name
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_C_1D_C_PS'

    if ( .not. allocated ( I % Step_1D ) ) then
      call Show ( 'Step_1D not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_1D_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_1D_C_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    call I % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

  end subroutine InitializeTemplate_C_1D_C_PS


  subroutine FinalizeTemplate_C_1D_C_PS ( I )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I

    if ( allocated ( I % Step_1D ) ) &
      deallocate ( I % Step_1D )

    call I % FinalizeTemplate_C_PS ( )

  end subroutine FinalizeTemplate_C_1D_C_PS


  subroutine ComputeCycle ( I )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I

    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerCycle )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
      call I % ComputeCycle_C_1D_C_ASC ( PS )
    end select

!    select type ( MS => I % MomentumSpace )
!    class is ( Bundle_SLL_ASC_CSLD_Form )
!      call I % ComputeCycle_BSLL_ASC_CSLD_1D ( MS )
!    end select

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeCycle


  subroutine InitializeStepTimers ( I, BaseLevel )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    if ( allocated ( I % Step_1D ) ) &
      call I % Step_1D % InitializeTimers ( BaseLevel )
    if ( allocated ( I % Step ) ) &
      call I % Step % InitializeTimers ( BaseLevel )

  end subroutine InitializeStepTimers


  subroutine ComputeTally ( I, ComputeChangeOption, IgnorabilityOption )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call I % ComputeTally_1D &
           ( ComputeChangeOption = ComputeChangeOption, &
             IgnorabilityOption  = IgnorabilityOption )

    if ( allocated ( I % Current_ASC ) ) then 
      associate ( CA => I % Current_ASC )
      call CA % ComputeTally &
             ( ComputeChangeOption = ComputeChangeOption, &
               IgnorabilityOption  = IgnorabilityOption )
      end associate !-- CA
    end if

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeTally


  subroutine ComputeCycle_C_1D_C_ASC ( I, PS )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I
    class ( Atlas_SC_Form ), intent ( inout ) :: &
      PS

    real ( KDR ) :: &
      TimeNew

    call I % PrepareCycle ( )
    call I % ComputeNewTime ( TimeNew )

    associate ( TimeStep => TimeNew - I % Time )    

    call I % PrepareStep_1D ( )
    call I % Step_1D % Compute ( I % Time, TimeStep )

    if ( allocated ( I % Step ) ) then
      call I % PrepareStep ( )
      call I % Step % Compute ( I % Time, TimeStep )
    end if
      
    I % iCycle = I % iCycle + 1
    I % Time = I % Time + TimeStep

    if ( I % WriteTimeExact ) then
      if ( I % Time == I % WriteTime ) &
        I % IsCheckpointTime = .true.
    else 
      if ( I % Time  >  I % WriteTime &
           .or. abs ( I % Time - I % WriteTime )  <  0.5_KDR * TimeStep ) &
        I % IsCheckpointTime = .true.
    end if

    end associate !-- TimeStep

  end subroutine ComputeCycle_C_1D_C_ASC


  subroutine PrepareStep_1D ( I )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I

  end subroutine PrepareStep_1D


  subroutine PrepareStep ( I )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I

  end subroutine PrepareStep


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    call I % ComputeTimeStepLocalTemplate ( TimeStepCandidate )

  end subroutine ComputeTimeStepLocal


  subroutine ComputeTimeStepLocalTemplate ( I, TimeStepCandidate )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    integer ( KDI ) :: &
      iC, &  !-- iCurrent
      N_PS
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( CurrentTemplate ), pointer :: &
      C
    class ( Current_ASC_Template ), pointer :: &
      CA

    N_PS  =  I % N_CURRENTS_1D + 1

    select type ( PS  =>  I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( CSL  =>  PS % Chart )
    class is ( Chart_SL_Template )

    G  =>  CSL % Geometry ( )

    do iC = 1, N_PS

      if ( iC == N_PS ) then
        if ( allocated ( I % Current_ASC ) ) then
          CA  =>  I % Current_ASC
        else
          CA  =>  null ( )
        end if
      else
        CA  =>  I % Current_ASC_Pointer ( iC )
      end if

      if ( associated ( CA ) ) then

        C => CA % Current ( )
        call I % ComputeTimeStepKernel_CSL &
               ( CSL % IsProperCell, &
                 C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
                 C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
                 C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
                 C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
                 C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
                 C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
                 G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
                 G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), & 
                 G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
                 G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
                 G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), & 
                 G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
                 G % Value ( :, G % COARSENING ( 2 ) ), &
                 G % Value ( :, G % COARSENING ( 3 ) ), &
                 CSL % nDimensions, TimeStepCandidate ( iC ) )
      else
        TimeStepCandidate ( iC )  =  huge ( 1.0_KDR )
      end if

    end do !-- iC

    end select !-- CSL
    end select !-- PS

    TimeStepCandidate = I % CourantFactor * TimeStepCandidate

    nullify ( G, C, CA )

  end subroutine ComputeTimeStepLocalTemplate


end module Integrator_C_1D_C_PS__Template