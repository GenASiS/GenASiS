module Step_RK_CS__Form

  !-- Step_RungeKutta_CurrentSet_Form

  use Basics
  use Algebra
  use Fields
  use Slopes
  use Step_RK_H__Form

  implicit none
  private

  type, public, extends ( Step_RK_H_Form ) :: Step_RK_CS_Form
    type ( FieldSetForm ), allocatable :: &
      Balanced, &
      Intermediate, &
      Solution
    type ( FieldSetElement ), dimension ( : ), allocatable :: &
      SolutionStage
    class ( CurrentSetForm ), pointer :: &
      CurrentSet
    class ( RiemannSolver_HLL_Form ), allocatable :: &
      RiemannSolver
  contains
    procedure, private, pass :: &
      Initialize_CS
    generic, public :: &
      Initialize => Initialize_CS
    procedure, public, pass :: &
      Show => Show_S
    final :: &
      Finalize
    procedure, private, pass :: &
      LoadSolution
    procedure, private, pass :: &
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, private, pass :: &
      ComputeStage
    procedure, private, pass :: &
      IncrementSolution
    procedure, private, pass :: &
      StoreSolution
  end type Step_RK_CS_Form


contains


  subroutine Initialize_CS ( S, CS, NameOption, A_Option, B_Option, C_Option )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    real ( KDR ), dimension ( 2 : , : ), intent ( in ), optional :: &
      A_Option
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      B_Option
    real ( KDR ), dimension ( 2 : ), intent ( in ), optional :: &
      C_Option

    integer ( KDI ) :: &
      iS  !-- iStage
    character ( 1 ) :: &
      StageNumber
    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_CS'

    Name  =  'Step_' // trim ( CS % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    S % CurrentSet  =>  CS

    call S % Step_RK_H_Form % Initialize &
           ( CS % Atlas, &
             NameOption = Name, &
             A_Option = A_Option, &
             B_Option = B_Option, &
             C_Option = C_Option )

    !-- Balanced

    allocate ( S % Balanced )
    associate ( Blncd  =>  S % Balanced )
    call Blncd % Initialize &
           ( CS, iaSelected = CS % iaBalanced, &
             NameOption = 'Balanced', &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- Blncd

    !-- Intermediate storage

    allocate ( S % Intermediate )
    associate ( Intmdt  =>  S % Intermediate )
    call Intmdt % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = 'Intermediate', &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- Intmdt

    !-- Solution storage

    allocate ( S % Solution )
    associate ( Sltn  =>  S % Solution )
    call Sltn % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = 'Solution', &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- Sltn

    !-- RiemannSolver

    if ( .not. allocated ( S % RiemannSolver ) ) then
      allocate ( S % RiemannSolver )
      associate ( RS  =>  S % RiemannSolver )
      call RS % Initialize ( CS )
      end associate !-- RSA
    end if !-- allocated RiemannSolver

    !-- Slopes

    if ( .not. allocated ( S % SlopeStage ) ) then
      associate ( nS  =>  S % nStages )
      allocate ( S % SlopeStage ( nS ) )
      do iS  =  1,  nS
        write ( StageNumber, fmt = '(i1.1)' ) iS
        allocate ( Slope_DFV_F_Form :: S % SlopeStage ( iS ) % Element )
        select type ( SA  =>  S % SlopeStage ( iS ) % Element )
          class is ( Slope_DFV_F_Form )
        call SA % Initialize &
               ( S % RiemannSolver, SuffixOption = StageNumber )
        end select !-- SA
      end do !-- iS
      end associate !-- nS
    end if !-- allocated SlopeStage

  end subroutine Initialize_CS


  subroutine Show_S ( S )

    class ( Step_RK_CS_Form ), intent ( in ) :: &
      S

    call S % Step_RK_H_Form % Show ( )

    call S % Solution % Show ( )
    call S % Intermediate % Show ( )
    call S % RiemannSolver % Show ( )

  end subroutine Show_S


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_CS_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % RiemannSolver ) ) &
      deallocate ( S % RiemannSolver )
    if ( allocated ( S % SolutionStage ) ) &
      deallocate ( S % SolutionStage )
    if ( allocated ( S % Solution ) ) &
      deallocate ( S % Solution )
    if ( allocated ( S % Intermediate ) ) &
      deallocate ( S % Intermediate )
    if ( allocated ( S % Balanced ) ) &
      deallocate ( S % Balanced )

    nullify ( S % CurrentSet )

  end subroutine Finalize


  subroutine LoadSolution ( S )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S

    associate &
      ( CS_B  =>  S % Balanced, &
        Y     =>  S % Solution )

    call CS_B % Copy ( Y )

    !-- For diagnostic I/O
    if ( allocated ( S % SolutionStage ) ) then
      associate ( Y_Stg  =>  S % SolutionStage ( 1 ) % Element )
      call CS_B % Copy ( Y_Stg )
      end associate !-- SltnStg
    end if

    end associate !-- CS_B, etc.

  end subroutine LoadSolution


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage

    if ( iS  >  1 ) then
      associate &
        ( Y    =>  S % Solution, &
          Y_I  =>  S % Intermediate )

      call Y % Copy ( Y_I )

      end associate !-- Y, etc.
    end if !-- iS > 1

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, dT, iK )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       A, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iK

    associate &
      ( Y_I  =>  S % Intermediate, &
        K    =>  S % SlopeStage ( iK ) % Element )

    call Y_I % MultiplyAdd ( K, dT * A )

    end associate !-- Y_I, etc.

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, T, iS, T_Option )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      T
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC  !-- iChart
    type ( TimerForm ), pointer :: &
      T_C, &
      T_EG

    associate ( K  =>  S % SlopeStage ( iS ) % Element )

    if ( iS  >  1 ) then

      associate ( K_1  =>  S % SlopeStage ( 1 ) % Element )
      call K % CloneTimers ( K_1 )
      end associate !-- K_1

      associate &
        (  Y_I  =>  S % Intermediate, &
          CS_B  =>  S % Balanced, &
          CS    =>  S % CurrentSet )

      call Y_I % Copy ( CS_B )
      call CS % ComputeFromBalanced ( )
      call CS % ApplyBoundaryConditions ( )
    
      !-- For diagnostic I/O
      if ( allocated ( S % SolutionStage ) ) then
        associate ( Y_Stg  =>  S % SolutionStage ( iS ) % Element )
        call Y_I % Copy ( Y_Stg )
        end associate !-- Y_Stg
      end if

      end associate !-- Y_I, etc.
 
    end if !-- iStage > 1

    !-- Compute slope

    if ( present ( T_Option ) ) then
      T_C  =>  K % Timer ( LevelOption = T_Option % Level + 1 )
      call T_C % Start ( )
    else
      T_C   =>  null ( )
    end if
    call K % Compute ( T_Option = T_C, iS_Option = iS )
    if ( associated ( T_C ) ) call T_C % Stop ( )

    !-- Slope ghost exchange

    if ( present ( T_Option ) ) then
      T_EG  =>  K % TimerGhost &
                  ( NameRootOption = K % TimerName, &
                    LevelOption = T_Option % Level + 1 )
      call T_EG % Start ( )
    else
      T_EG  =>  null ( )
    end if
    call K % ExchangeGhostData ( )
    if ( associated ( T_EG ) ) call T_EG % Stop ( )

    end associate !-- K

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, dT, iS )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       B, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iC  !-- iChart

    associate &
      ( Y  =>  S % Solution, &
        K  =>  S % SlopeStage ( iS ) % Element )

    call Y % MultiplyAdd ( K, dT * B )

    end associate !-- Y, etc.

  end subroutine IncrementSolution


  subroutine StoreSolution ( S )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S

    associate &
      (  Y    =>  S % Solution, &
        CS_B  =>  S % Balanced, &
        CS    =>  S % CurrentSet )

    call  Y % Copy ( CS_B )
    call CS % ComputeFromBalanced ( )
    call CS % ApplyBoundaryConditions ( )
  
    end associate !-- Y, etc.
 
  end subroutine StoreSolution


end module Step_RK_CS__Form
