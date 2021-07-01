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
    procedure, private, pass :: &
      Show_S
    generic, public :: &
      Show => Show_S
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
      call S % ComputeStage ( T, iS, TimerLevelOption = Timer_CS % Level + 1 )
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
      TimerName  =  trim ( Timer % Name ) // '_SF' 
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


  subroutine ComputeStage ( S, T, iS, TimerLevelOption )

      class ( Step_RK_H_Form ), intent ( inout ) :: &
        S
      real ( KDR ), intent ( in ) :: &
        T
      integer ( KDI ), intent ( in ) :: &
        iS  !-- iStage
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

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
