!-- Step_RK is a template for a Runge-Kutta time step.

module Step_RK__Template

  !-- Step_RungeKutta_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics

  implicit none
  private

  type, public, abstract :: Step_RK_Template
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerStep = 0, &
      iTimerLoadInitial = 0, &
      iTimerTemplate = 0, &
      iTimerInitializeIntermediate = 0, &
      iTimerIncrementIntermediate = 0, &
      iTimerStage = 0, &
      iTimerIncrementSolution = 0, &
      nStages
    real ( KDR ), dimension ( : ), allocatable :: &
      C, & 
      B
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      A
    character ( LDF ) :: &
      Type = '', &
      Name = ''
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      InitializeTimers
    procedure, public, pass :: &
      ComputeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
    procedure ( II ), private, pass, deferred :: &
      InitializeIntermediate
    procedure ( II_A_iK ), private, pass, deferred :: &
      IncrementIntermediate
    procedure ( CS ), private, pass, deferred :: &
      ComputeStage
    procedure ( IS_B_iS ), private, pass, deferred :: &
      IncrementSolution
  end type Step_RK_Template

  abstract interface

    subroutine IT ( S, BaseLevel )
      use Basics
      import Step_RK_Template
      class ( Step_RK_Template ), intent ( inout ) :: &
        S
      integer ( KDI ), intent ( in ) :: &
        BaseLevel
    end subroutine IT

    subroutine II ( S )
      import Step_RK_Template
      class ( Step_RK_Template ), intent ( inout ) :: &
        S
    end subroutine II

    subroutine II_A_iK ( S, A, iK )
      use Basics
      import Step_RK_Template
      class ( Step_RK_Template ), intent ( inout ) :: &
        S
      real ( KDR ), intent ( in ) :: &
        A
      integer ( KDI ), intent ( in ) :: &
        iK
    end subroutine II_A_iK

    subroutine CS ( S, Time, TimeStep, iStage )
      use Basics
      import Step_RK_Template
      class ( Step_RK_Template ), intent ( inout ) :: &
        S
      real ( KDR ), intent ( in ) :: &
        Time, &
        TimeStep
      integer ( KDI ), intent ( in ) :: &
        iStage
    end subroutine CS

    subroutine IS_B_iS ( S, B, iS )
      use Basics
      import Step_RK_Template
      class ( Step_RK_Template ), intent ( inout ) :: &
        S
      real ( KDR ), intent ( in ) :: &
        B
      integer ( KDI ), intent ( in ) :: &
        iS
    end subroutine IS_B_iS

  end interface


contains


  subroutine InitializeTemplate ( S, NameSuffix, A, B, C )

    class ( Step_RK_Template ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ) :: &
      NameSuffix
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iS  !-- iStage

    S % IGNORABILITY = CONSOLE % INFO_1

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK' 

    S % Name = 'Step_' // trim ( NameSuffix ) 

    call Show ( 'Initializing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    S % nStages = size ( B )
    associate ( nS => S % nStages )

    allocate ( S % A ( 2 : nS ) )
    do iS = 2, nS
      call S % A ( iS ) % Initialize ( iS - 1 )
      S % A ( iS ) % Value = A ( iS, 1 : iS - 1 )
    end do !-- iS

    allocate ( S % B ( nS ) )
    S % B = B

    allocate ( S % C ( 2 : nS ) )
    S % C = C

    end associate !-- nS

  end subroutine InitializeTemplate


  subroutine InitializeTimers ( S, BaseLevel )

    class ( Step_RK_Template ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    if ( S % iTimerStep > 0  &
         .or.  BaseLevel > PROGRAM_HEADER % TimerLevel ) &
      return

    call PROGRAM_HEADER % AddTimer &
           ( S % Name, S % iTimerStep, &
             Level = BaseLevel )
      call PROGRAM_HEADER % AddTimer &
             ( 'LoadInitial', S % iTimerLoadInitial, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'Template', S % iTimerTemplate, &
               Level = BaseLevel + 1 )
        call PROGRAM_HEADER % AddTimer &
               ( 'InitializeIntermediate', S % iTimerInitializeIntermediate, &
                 Level = BaseLevel + 2 )
        call PROGRAM_HEADER % AddTimer &
               ( 'IncrementIntermediate', S % iTimerIncrementIntermediate, &
                 Level = BaseLevel + 2 )
        call PROGRAM_HEADER % AddTimer &
               ( 'Stage', S % iTimerStage, &
                 Level = BaseLevel + 2 )
        call PROGRAM_HEADER % AddTimer &
               ( 'IncrementSolution', S % iTimerIncrementSolution, &
                 Level = BaseLevel + 2 )

  end subroutine InitializeTimers


!   subroutine ComputeTemplate ( S, Solution, Time, TimeStep )

!     class ( Step_RK_Template ), intent ( inout ) :: &
!       S
!     type ( VariableGroupForm ), dimension ( : ), intent ( inout ) :: &
!       Solution
!     real ( KDR ), intent ( in ) :: &
!       Time, &
!       TimeStep

!     integer ( KDI ) :: &
!       iS, &  !-- iStage
!       iG, &  !-- iGroup
!       iK     !-- iIncrement
!     integer ( KDI ), dimension ( size ( Solution ) ) :: &
!       nValues, &
!       nEquations
!     type ( VariableGroupForm ), dimension ( size ( Solution ) ) :: &
!       Y  !-- Argument of right-hand side
!     type ( VariableGroupForm ), &
!       dimension ( size ( Solution ), S % nStages ) :: &
!         K  !-- Increments

!     associate ( nGroups => size ( Solution ) ) !-- nGroups

!     do iG = 1, nGroups
!       nEquations ( iG )  =  Solution ( iG ) % nVariables
!       nValues    ( iG )  =  Solution ( iG ) % nValues
!       call Y ( iG ) % Initialize ( [ nValues ( iG ), nEquations ( iG ) ] )
!     end do !-- iG

!     !-- Compute increments

!     do iS = 1, S % nStages
!       do iG = 1, nGroups

!         call K ( iG, iS ) % Initialize &
!                ( [ nValues ( iG ), nEquations ( iG ) ] )

!         associate &
!           ( YV => Y ( iG ) % Value, &
!             SV => Solution ( iG ) % Value )

!         YV = SV
!         do iK = 1, iS - 1
!           associate &
!             ( KV => K ( iG, iK ) % Value, &
!               A  => S % A ( iS ) % Value ( iK ) )

! !         YV  =  YV  +  A * KV  
!           call MultiplyAdd ( YV, KV, A )

!           end associate !-- KV, etc.
!         end do !-- iK

!         call S % ComputeIncrement ( K, Y, Time, TimeStep, iS, iG ) 

!         end associate !-- YV, etc.

!       end do !-- iG
!     end do !-- iS

!     !-- Assemble increments

!     do iS = 1, S % nStages
!       do iG = 1, nGroups
!         associate &
!           ( SV => Solution ( iG ) % Value, &
!             KV => K ( iG, iS ) % Value, &
!             B  => S % B ( iS ) )
! !       SV  =  SV  +  B * KV
!         call MultiplyAdd ( SV, KV, B )
!         end associate !-- SV, etc.
!       end do !-- iG
!     end do !-- iS

!     end associate !-- nGroups

!   end subroutine ComputeTemplate


  subroutine ComputeTemplate ( S, Time, TimeStep )

    class ( Step_RK_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep

    integer ( KDI ) :: &
      iS, &  !-- iStage
      iK     !-- iIncrement
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerTemplate )
    if ( associated ( Timer ) ) call Timer % Start ( )

    !-- On entry, Solution = Y_N (old value)

    !-- Compute stages
    do iS = 1, S % nStages

      call Show ( 'Computing a stage', S % IGNORABILITY + 2 )
      call Show ( iS, 'iStage', S % IGNORABILITY + 2 )

      !-- Set Y  =  Solution
      call S % InitializeIntermediate ( )

      do iK = 1, iS - 1
        associate ( A  =>  S % A ( iS ) % Value ( iK ) )
        !-- Set Y  =  Y  +  A * K ( iK )
        call S % IncrementIntermediate ( A, iK )
        end associate !-- A
      end do !-- iK

      call S % ComputeStage ( Time, TimeStep, iS )

    end do !-- iS

    !-- Assemble stages
    do iS = 1, S % nStages
      associate ( B  =>  S % B ( iS ) )
      !-- Set Solution  =  Solution  +  B * K ( iS )
      call S % IncrementSolution ( B, iS )
      end associate !-- B
    end do !-- iS

    !-- On exit, Solution  =  Y_(N+1) (new value)

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeTemplate


  impure elemental subroutine FinalizeTemplate ( S )

    class ( Step_RK_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % A ) ) &
      deallocate ( S % A )
    if ( allocated ( S % B ) ) &
      deallocate ( S % B )
    if ( allocated ( S % C ) ) &
      deallocate ( S % C )

    if ( S % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

  end subroutine FinalizeTemplate


end module Step_RK__Template
