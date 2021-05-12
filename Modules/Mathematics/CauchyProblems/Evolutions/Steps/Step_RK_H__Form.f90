module Step_RK_H__Form
  
  !-- Step_RungeKutta_Header_Form

  use Basics
  
  implicit none
  private

  type, public :: Step_RK_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimer       = 0, &
      iTimer_LI    = 0, &  !-- LoadInitial
      iTimer_II    = 0, &  !-- InitializeIntermediate
      iTimer_II_A  = 0, &  !-- IncrementIntermediate
      iTimer_CS    = 0, &  !-- ComputeStage
      iTimer_IS_B  = 0, &  !-- IncrementSolution
      iTimer_SF    = 0, &  !-- StoreFinal
      nEquations, &
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
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    procedure, public, pass :: &
      Compute
    procedure, private, pass :: &
      Show_S
    generic, public :: &
      Show => Show_S
    final :: &
      Finalize
    procedure, private, pass :: &
      LoadSolution
    procedure, private, pass :: &
      StoreSolution
    procedure, private, pass :: &
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, private, pass :: &
      ComputeStage
    procedure, private, pass :: &
      IncrementSolution
  end type Step_RK_H_Form


contains


  subroutine Initialize_H ( S, A, B, C, nEquations, NameOption )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nEquations
    character ( * ), intent ( in ), optional :: &
      NameOption

    integer ( KDI ) :: &
      iS

    S % IGNORABILITY  =  CONSOLE % INFO_1

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK' 

    S % Name = 'Step'
    if ( present ( NameOption ) ) &
      S % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    S % nEquations  =  nEquations

    S % nStages  =  size ( B )
    associate ( nS  =>  S % nStages )

    allocate ( S % A ( 2 : nS ) )
    do iS  =  2,  nS
      call S % A ( iS ) % Initialize ( iS - 1 )
      S % A ( iS ) % Value  =  A ( iS, 1 : iS - 1 )
    end do !-- iS

    allocate ( S % B ( nS ) )
    S % B  =  B

    allocate ( S % C ( 2 : nS ) )
    S % C  =  C

    end associate !-- nS

  end subroutine Initialize_H


  subroutine Compute ( S, T, dT, TimerLevelOption )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      T, &
      dT
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    integer ( KDI ) :: &
      iS, &  !-- iStage
      iK     !-- iIncrement
    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_LI, &
      Timer_II, &
      Timer_II_A, &
      Timer_CS, &
      Timer_IS_B, &
      Timer_SF

    associate ( iT  =>  S % iTimer )
    if ( iT == 0 ) then
      TimerName  =  S % Name
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    Timer  =>  PROGRAM_HEADER % TimerPointer ( S % iTimer )
    call Timer % Start ( )

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    !-- Set  Solution  =  Y_N  (old value)

    associate ( iT  =>  S % iTimer_LI )
    if ( iT == 0 ) then
      TimerName  =  trim ( Timer % Name ) // '_LI' 
      call PROGRAM_HEADER % AddTimer &
             ( TimerName, iT, Level = Timer % Level + 1 )
    end if
    end associate !-- iT
    Timer_LI  =>  PROGRAM_HEADER % TimerPointer ( S % iTimer_LI )

    call Timer_LI % Start ( )
    call S % LoadSolution ( )
    call Timer_LI % Stop ( )

    !-- Compute stages

    do iS = 1, S % nStages

      call Show ( 'Computing a stage', S % IGNORABILITY + 3 )
      call Show ( iS, 'iStage', S % IGNORABILITY + 3 )

      !-- Set  Y  =  Solution

      associate ( iT  =>  S % iTimer_II )
      if ( iT == 0 ) then
        TimerName  =  trim ( Timer % Name ) // '_II' 
        call PROGRAM_HEADER % AddTimer &
               ( TimerName, iT, Level = Timer % Level + 1 )
      end if
      end associate !-- iT
      Timer_II  =>  PROGRAM_HEADER % TimerPointer ( S % iTimer_II )

      call Timer_II % Start ( )
      call S % InitializeIntermediate ( iS )
      call Timer_II % Stop ( )

      !-- Loop: Set  Y  =   Y  +  A * K ( iK )

      associate ( iT  =>  S % iTimer_II_A )
      if ( iT == 0 ) then
        TimerName  =  trim ( Timer % Name ) // '_II_A' 
        call PROGRAM_HEADER % AddTimer &
               ( TimerName, iT, Level = Timer % Level + 1 )
      end if
      end associate !-- iT
      Timer_II_A  =>  PROGRAM_HEADER % TimerPointer ( S % iTimer_II_A )

      call Timer_II_A % Start ( )
      do iK = 1, iS - 1
        associate ( A  =>  S % A ( iS ) % Value ( iK ) )
        !-- Set Y  =  Y  +  dT * A * K ( iK )
        call S % IncrementIntermediate ( A, dT, iK )
        end associate !-- A
      end do !-- iK
      call Timer_II_A % Stop ( )

      associate ( iT  =>  S % iTimer_CS )
      if ( iT == 0 ) then
        TimerName  =  trim ( Timer % Name ) // '_CS' 
        call PROGRAM_HEADER % AddTimer &
               ( TimerName, iT, Level = Timer % Level + 1 )
      end if
      end associate !-- iT
      Timer_CS  =>  PROGRAM_HEADER % TimerPointer ( S % iTimer_CS )

      call Timer_CS % Start ( )
      call S % ComputeStage ( T, iS )
      call Timer_CS % Stop ( )

    end do !-- iS

    associate ( iT  =>  S % iTimer_IS_B )
    if ( iT == 0 ) then
      TimerName  =  trim ( Timer % Name ) // '_IS_B' 
      call PROGRAM_HEADER % AddTimer &
             ( TimerName, iT, Level = Timer % Level + 1 )
    end if
    end associate !-- iT
    Timer_IS_B  =>  PROGRAM_HEADER % TimerPointer ( S % iTimer_IS_B )

    call Timer_IS_B % Start ( )
    !-- Assemble stages
    do iS = 1, S % nStages
      associate ( B  =>  S % B ( iS ) )
      !-- Set Solution  =  Solution  +  dT * B * K ( iS )
      call S % IncrementSolution ( B, dT, iS )
      end associate !-- B
    end do !-- iS
    call Timer_IS_B % Stop ( )

    associate ( iT  =>  S % iTimer_SF )
    if ( iT == 0 ) then
      TimerName  =  trim ( Timer % Name ) // '_IS_B' 
      call PROGRAM_HEADER % AddTimer &
             ( TimerName, iT, Level = Timer % Level + 1 )
    end if
    end associate !-- iT
    Timer_SF  =>  PROGRAM_HEADER % TimerPointer ( S % iTimer_SF )

    !-- On exit, Solution  =  Y_(N+1) (new value)
    call Timer_SF % Start ( )
    call S % StoreSolution ( )
    call Timer_SF % Stop ( )

    call Timer % Stop ( )

  end subroutine Compute


  subroutine Show_S ( S )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

   integer ( KDI ) :: &
     iA
   character ( 1 ) :: &
     Index

    call Split ( S % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    call Show ( S % nEquations, 'nEquations', S % IGNORABILITY )
    call Show ( S % nStages, 'nStages', S % IGNORABILITY )

    do iA  =  2, S % nStages
      write ( Index, fmt = '(i1.1)' ) iA
      call Show ( S % A ( iA ) % Value, 'A ( ' // Index // ' )' )
    end do !-- iA

    call Show ( S % B, 'B' )
    call Show ( S % C, 'C', lRealOption = 2 )

  end subroutine Show_S


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_H_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % A ) ) &
      deallocate ( S % A )
    if ( allocated ( S % B ) ) &
      deallocate ( S % B )
    if ( allocated ( S % C ) ) &
      deallocate ( S % C )

    if ( S % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

  end subroutine Finalize


  subroutine LoadSolution ( S )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S

    call Show ( 'LoadSolution must be overridden', CONSOLE % ERROR )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % ERROR )
    call Show ( 'LoadSolution', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine LoadSolution


  subroutine StoreSolution ( S )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S

    call Show ( 'StoreSolution must be overridden', CONSOLE % ERROR )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % ERROR )
    call Show ( 'StoreSolution', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine StoreSolution


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS

    call Show ( 'InitializeIntermediate must be overridden', CONSOLE % ERROR )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % ERROR )
    call Show ( 'InitializeIntermediate', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, dT, iK )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       A, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iK

    call Show ( 'IncrementIntermediate must be overridden', CONSOLE % ERROR )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % ERROR )
    call Show ( 'IncrementIntermediate', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, T, iS )

      class ( Step_RK_H_Form ), intent ( inout ) :: &
        S
      real ( KDR ), intent ( in ) :: &
        T
      integer ( KDI ), intent ( in ) :: &
        iS  !-- iStage

    call Show ( 'ComputeStage must be overridden', CONSOLE % ERROR )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % ERROR )
    call Show ( 'ComputeStage', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, dT, iS )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       B, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    call Show ( 'IncrementSolution must be overridden', CONSOLE % ERROR )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % ERROR )
    call Show ( 'IncrementSolution', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine IncrementSolution


end module Step_RK_H__Form
