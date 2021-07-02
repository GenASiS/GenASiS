module Step_RK_H__Form
  
  !-- Step_RungeKutta_Header_Form

  use Basics
  
  implicit none
  private

  type, public :: Step_RK_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nStages = 0
    integer ( KDI ) :: &
      iTimer      = 0, &
      iTimer_LI   = 0, &  !-- LoadInitial
      iTimer_II   = 0, &  !-- InitializeIntermediate
      iTimer_II_A = 0, &  !-- IncrementIntermediate
      iTimer_CS   = 0, &  !-- ComputeStage
      iTimer_IS_B = 0, &  !-- IncrementSolution
      iTimer_SF   = 0     !-- StoreFinal
    real ( KDR ), dimension ( : ), allocatable :: &
      C, &  !-- RungeKutta nodes
      B     !-- RungeKutta weights
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      A  !-- RungeKutta matrix
    character ( LDF ) :: &
      Type = '', &
      Name = ''
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    procedure, public, pass :: &
      Show => Show_S
    procedure, public, pass :: &
      Timer
    procedure, private, pass :: &
      Timer_LI
    procedure, private, pass :: &
      Timer_II
    procedure, private, pass :: &
      Timer_II_A
    procedure, private, pass :: &
      Timer_CS
    procedure, private, pass :: &
      Timer_IS_B
    procedure, private, pass :: &
      Timer_SF
    procedure, public, pass :: &
      Compute
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


  subroutine Initialize_H ( S, NameOption, A_Option, B_Option, C_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ), optional :: &
      NameOption
    real ( KDR ), dimension ( 2 : , : ), intent ( in ), optional :: &
      A_Option  !-- RungeKutta matrix
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      B_Option  !-- RungeKutta weights
    real ( KDR ), dimension ( 2 : ), intent ( in ), optional :: &
      C_Option  !-- RungeKutta nodes

    integer ( KDI ) :: &
      iS
    real ( KDR ), dimension ( :, : ), allocatable :: &
      A  !-- RungeKutta matrix
    real ( KDR ), dimension ( : ), allocatable :: &
      B  !-- RungeKutta weights
    real ( KDR ), dimension ( : ), allocatable :: &
      C  !-- RungeKutta nodes

    S % IGNORABILITY  =  CONSOLE % INFO_1

    if ( present ( A_Option ) ) then
      allocate ( A ( 2 : ubound ( A_Option, dim = 1 ), &
                     size ( A_Option, dim = 2 ) ) )
      A  =  A_Option
    else
      allocate ( A ( 2 : 2, 1 : 1 ) )
      A           =  0.0_KDR
      A ( 2, 1 )  =  1.0_KDR
    end if

    if ( present ( B_Option ) ) then
      allocate ( B, source = B_Option )
    else
      allocate ( B ( 1 : 2 ) )
      B ( 1 )  =  0.5_KDR
      B ( 2 )  =  0.5_KDR
    end if

    if ( present ( C_Option ) ) then
      allocate ( C ( 2 : ubound ( C_Option, dim = 1 ) ) )
      C  =  C_Option
    else
      allocate ( C ( 2 : 2 ) )
      C ( 2 )  =  1.0_KDR
    end if

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK' 

    S % Name = 'Step'
    if ( present ( NameOption ) ) &
      S % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

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

    call Show ( S % nStages, 'nStages', S % IGNORABILITY )

    do iA  =  2, S % nStages
      write ( Index, fmt = '(i1.1)' ) iA
      call Show ( S % A ( iA ) % Value, 'A ( ' // Index // ' )' )
    end do !-- iA

    call Show ( S % B, 'B' )
    call Show ( S % C, 'C', lRealOption = 2 )

  end subroutine Show_S


  function Timer ( S, LevelOption ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimer )

    if ( iT == 0 ) then
      TimerName  =  S % Name
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer


  function Timer_LI ( S, LevelOption ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimer_LI )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % Name ) // '_LdIntl'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_LI


  function Timer_II ( S, LevelOption ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimer_II )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % Name ) // '_IntlzIntmdt'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_II


  function Timer_II_A ( S, LevelOption ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimer_II_A )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % Name ) // '_IncrmntIntmdt'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_II_A


  function Timer_CS ( S, LevelOption ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimer_CS )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % Name ) // '_CmptStg'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_CS


  function Timer_IS_B ( S, LevelOption ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimer_IS_B )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % Name ) // '_IncrmntSltn'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_IS_B


  function Timer_SF ( S, LevelOption ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimer_SF )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % Name ) // '_StrFnl'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_SF


  subroutine Compute ( S, T, dT, T_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      T, &
      dT
    type ( TimerForm ), intent ( in ), pointer, optional :: &
      T_Option

    integer ( KDI ) :: &
      iS, &  !-- iStage
      iK     !-- iIncrement
    type ( TimerForm ), pointer :: &
      T_LI, &
      T_II, &
      T_II_A, &
      T_CS, &
      T_IS_B, &
      T_SF

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    if ( present ( T_Option ) ) then
      T_LI    =>  S % Timer_LI   ( LevelOption = T_Option % Level + 1 )
      T_II    =>  S % Timer_II   ( LevelOption = T_Option % Level + 1 )
      T_II_A  =>  S % Timer_II_A ( LevelOption = T_Option % Level + 1 )
      T_CS    =>  S % Timer_CS   ( LevelOption = T_Option % Level + 1 )
      T_IS_B  =>  S % Timer_IS_B ( LevelOption = T_Option % Level + 1 )
      T_SF    =>  S % Timer_SF   ( LevelOption = T_Option % Level + 1 )
    else
      T_LI    =>  null ( )
      T_II    =>  null ( )
      T_II_A  =>  null ( )
      T_CS    =>  null ( )
      T_IS_B  =>  null ( )
      T_SF    =>  null ( )
    end if

    !-- Set  Solution  =  Y_N  (old value)

    if ( associated ( T_LI ) ) call T_LI % Start ( )
    call S % LoadSolution ( )
    if ( associated ( T_LI ) ) call T_LI % Stop ( )

    !-- Compute stages

    do iS = 1, S % nStages

      call Show ( 'Computing a stage', S % IGNORABILITY + 3 )
      call Show ( iS, 'iStage', S % IGNORABILITY + 3 )

      !-- Set  Y  =  Solution

      if ( associated ( T_II ) ) call T_II % Start ( )
      call S % InitializeIntermediate ( iS )
      if ( associated ( T_II ) ) call T_II % Stop ( )

      !-- Loop: Set  Y  =   Y  +  A * K ( iK )

      if ( associated ( T_II_A ) ) call T_II_A % Start ( )
      do iK = 1, iS - 1
        associate ( A  =>  S % A ( iS ) % Value ( iK ) )
        !-- Set Y  =  Y  +  dT * A * K ( iK )
        call S % IncrementIntermediate ( A, dT, iK )
        end associate !-- A
      end do !-- iK
      if ( associated ( T_II_A ) ) call T_II_A % Stop ( )

      if ( associated ( T_CS ) ) then
        call T_CS % Start ( )
        call S % ComputeStage ( T, iS, T_Option = T_CS )
        call T_CS % Stop ( )
      else
        call S % ComputeStage ( T, iS )
      end if

    end do !-- iS

    if ( associated ( T_IS_B ) ) call T_IS_B % Start ( )
    !-- Assemble stages
    do iS = 1, S % nStages
      associate ( B  =>  S % B ( iS ) )
      !-- Set Solution  =  Solution  +  dT * B * K ( iS )
      call S % IncrementSolution ( B, dT, iS )
      end associate !-- B
    end do !-- iS
    if ( associated ( T_IS_B ) ) call T_IS_B % Stop ( )

    !-- On exit, Solution  =  Y_(N+1) (new value)
    if ( associated ( T_SF ) ) call T_SF % Start ( )
    call S % StoreSolution ( )
    if ( associated ( T_SF ) ) call T_SF % Stop ( )

  end subroutine Compute


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

    call Show ( 'LoadSolution should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'LoadSolution', 'subroutine', CONSOLE % WARNING )

  end subroutine LoadSolution


  subroutine StoreSolution ( S )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S

    call Show ( 'StoreSolution should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'StoreSolution', 'subroutine', CONSOLE % WARNING )

  end subroutine StoreSolution


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS

    call Show ( 'InitializeIntermediate should be overridden', &
                CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'InitializeIntermediate', 'subroutine', CONSOLE % WARNING )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, dT, iK )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       A, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iK

    call Show ( 'IncrementIntermediate should be overridden', &
                CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'IncrementIntermediate', 'subroutine', CONSOLE % WARNING )

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, T, iS, T_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      T
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( in ), pointer, optional :: &
      T_Option

    call Show ( 'ComputeStage should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'ComputeStage', 'subroutine', CONSOLE % WARNING )

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, dT, iS )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       B, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    call Show ( 'IncrementSolution should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'IncrementSolution', 'subroutine', CONSOLE % WARNING )

  end subroutine IncrementSolution


end module Step_RK_H__Form
