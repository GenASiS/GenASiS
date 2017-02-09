!-- Step_RK is a template for a Runge-Kutta time step.

module Step_RK__Template

  !-- Step_RungeKutta_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Operations

  implicit none
  private

  type, public, abstract :: Step_RK_Template
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerComputeStep = 0, &
      iTimerComputeIncrement = 0, &
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
      ComputeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
    procedure ( CI ), private, pass, deferred :: &
      ComputeIncrement
  end type Step_RK_Template

  abstract interface

    subroutine CI ( S, K, Y, Time, TimeStep, iStage, iGroup )
      use Basics
      import Step_RK_Template
      class ( Step_RK_Template ), intent ( inout ) :: &
        S
      type ( VariableGroupForm ), dimension ( :, : ), intent ( inout ) :: &
        K
      type ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
        Y
      real ( KDR ), intent ( in ) :: &
        Time, &
        TimeStep
      integer ( KDI ), intent ( in ) :: &
        iStage, &
        iGroup
    end subroutine CI

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
    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    S % IGNORABILITY = CONSOLE % INFO_3

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK' 

    call Split ( S % Type, ' ', TypeWord )
    S % Name = trim ( TypeWord ( 2 ) ) // '_' // trim ( NameSuffix ) 

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

    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeStep', S % iTimerComputeStep )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeIncrement', S % iTimerComputeIncrement )

  end subroutine InitializeTemplate


  subroutine ComputeTemplate ( S, Solution, Time, TimeStep )

    class ( Step_RK_Template ), intent ( inout ) :: &
      S
    type ( VariableGroupForm ), dimension ( : ), intent ( inout ) :: &
      Solution
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep

    integer ( KDI ) :: &
      iS, &  !-- iStage
      iG, &  !-- iGroup
      iK     !-- iIncrement
    integer ( KDI ), dimension ( size ( Solution ) ) :: &
      nValues, &
      nEquations
    type ( VariableGroupForm ), dimension ( size ( Solution ) ) :: &
      Y  !-- Argument of right-hand side
    type ( VariableGroupForm ), &
      dimension ( size ( Solution ), S % nStages ) :: &
        K  !-- Increments

    associate ( nGroups => size ( Solution ) ) !-- nGroups

    do iG = 1, nGroups
      nEquations ( iG )  =  Solution ( iG ) % nVariables
      nValues    ( iG )  =  Solution ( iG ) % nValues
      call Y ( iG ) % Initialize ( [ nValues ( iG ), nEquations ( iG ) ] )
    end do !-- iG

    !-- Compute increments

    do iS = 1, S % nStages
      do iG = 1, nGroups

        call K ( iG, iS ) % Initialize &
               ( [ nValues ( iG ), nEquations ( iG ) ] )

        associate &
          ( YV => Y ( iG ) % Value, &
            SV => Solution ( iG ) % Value )

        YV = SV
        do iK = 1, iS - 1
          associate &
            ( KV => K ( iG, iK ) % Value, &
              A  => S % A ( iS ) % Value ( iK ) )

!         YV  =  YV  +  A * KV  
          call MultiplyAdd ( YV, KV, A )

          end associate !-- KV, etc.
        end do !-- iK

        call S % ComputeIncrement ( K, Y, Time, TimeStep, iS, iG ) 

        end associate !-- YV, etc.

      end do !-- iG
    end do !-- iS

    !-- Assemble increments

    do iS = 1, S % nStages
      do iG = 1, nGroups
        associate &
          ( SV => Solution ( iG ) % Value, &
            KV => K ( iG, iS ) % Value, &
            B  => S % B ( iS ) )
!       SV  =  SV  +  B * KV
        call MultiplyAdd ( SV, KV, B )
        end associate !-- SV, etc.
      end do !-- iG
    end do !-- iS

    end associate !-- nGroups

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
