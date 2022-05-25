module Poisson_H__Form

  !-- Poisson_Header__Form

  use Basics
  use Fields
  use Laplacian_M_H__Form
  
  implicit none
  private

  type, public :: Poisson_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimer = 0, &
      iTimer_CM = 0, &  !-- CombineMoments
      iTimer_CS = 0, &  !-- ClearSolution
      iTimer_LS = 0, &  !-- LocalSolution
      iTimer_C  = 0, &  !-- Coarsening
      iTimer_ES = 0, &  !-- ExchangeSolution
      iTimer_BS = 0, &  !-- BoundarySolution
      nEquations = 0, &
      MaxDegree = 0
    logical ( KDL ) :: &
      Coarsen
    character ( LDF ) :: &
      Type = '', &
      Name = '', &
      SolverType = ''
    type ( Coarsening_C_Form ), allocatable :: &
      Coarsening
    class ( Laplacian_M_H_Form ), allocatable :: &
      Laplacian_M
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, public, pass :: &
      Show => Show_P
    procedure, public, pass :: &
      Timer
    procedure, public, pass :: &
      Solve
    final :: &
      Finalize
    procedure, private, pass :: &
      Solve_M
    procedure, private, pass :: &
      CombineMoments
    procedure, private, pass :: &
      CombineMomentsLocal
    procedure, private, pass :: &
      ExchangeSolution
    procedure, private, pass :: &
      ApplyBoundarySolution
  end type Poisson_H_Form


contains


  subroutine Initialize_H &
               ( P, G, SolverType, MaxDegreeOption, nEquationsOption )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    character ( * ), intent ( in ) :: &
      SolverType
    integer ( KDI ), intent ( in ), optional :: &
      MaxDegreeOption, &
      nEquationsOption

    P % IGNORABILITY  =  G % IGNORABILITY

    if ( P % Type == '' ) &
      P % Type = 'a Poisson' 

    P % Name = 'Poisson'

    call Show ( 'Initializing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )

    P % nEquations  =  1
    if ( present ( nEquationsOption ) ) &
      P % nEquations = nEquationsOption

    P % SolverType  =  SolverType

    P % MaxDegree  =  0
    if ( present ( MaxDegreeOption ) ) &
      P % MaxDegree  =  MaxDegreeOption

    P % Coarsen  =  .true.
    call PROGRAM_HEADER % GetParameter ( P % Coarsen, 'Coarsen' )
    if ( P % Coarsen ) then
      allocate ( P % Coarsening )
      associate ( C  =>  P % Coarsening )
      call C % Initialize ( G )
      end associate !-- C
    end if

  end subroutine Initialize_H


  subroutine Show_P ( P )

    class ( Poisson_H_Form ), intent ( in ) :: &
      P

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( P % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )
    call Show ( P % SolverType, 'SolverType', P % IGNORABILITY )

    if ( allocated ( P % Laplacian_M ) ) &
      call P % Laplacian_M % Show ( )

  end subroutine Show_P


  function Timer ( P, Level ) result ( T )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = P % iTimer, &
               Name = trim ( P % Name ) // '_Slv', &
               Level = Level )

  end function Timer


  subroutine Solve ( P, Solution, Source, T_Option )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution
    class ( FieldSetForm ), intent ( inout ) :: &
      Source
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select case ( trim ( P % SolverType ) )
    case ( 'MULTIPOLE' )

      call P % Solve_M ( Solution, Source, T_Option = T_Option )
   
    case default
      call Show ( 'Solver type not supported', CONSOLE % ERROR )
      call Show ( P % SolverType, 'Type', CONSOLE % ERROR )
      call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- SolverType

  end subroutine Solve


  impure elemental subroutine Finalize ( P )

    type ( Poisson_H_Form ), intent ( inout ) :: &
      P

    if ( allocated ( P % Laplacian_M ) ) &
      deallocate ( P % Laplacian_M )
    if ( allocated ( P % Coarsening ) ) &
      deallocate ( P % Coarsening )

    if ( P % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )
    
  end subroutine Finalize


  subroutine Solve_M ( P, Solution, Source, T_Option )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution
    class ( FieldSetForm ), intent ( inout ) :: &
      Source
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( TimerForm ), pointer :: &
      T_L, &
      T_CM

    call Show ( 'Poisson solve, multipole', P % IGNORABILITY + 3 )
    call Show ( P % Name, 'Name', P % IGNORABILITY + 3 )

    if ( allocated ( P % Laplacian_M ) ) then
      associate ( L  =>  P % Laplacian_M )
      if ( present ( T_Option ) ) then

        T_L  =>  L % Timer ( Level = T_Option % Level + 1 )
        call T_L % Start ( )
        call L % ComputeMoments ( Source, T_Option = T_L )
        call T_L % Stop ( )

        T_CM  =>  PROGRAM_HEADER % Timer &
                    ( Handle = P % iTimer_CM, &
                      Name = trim ( P % Name ) // '_CmbnMmnts', &
                      Level = T_Option % Level + 1 )
        call T_CM % Start ( )
        call P % CombineMoments ( Solution, T_Option = T_CM )
        call T_CM % Stop ( )

      else
        call L % ComputeMoments ( Source )
        call P % CombineMoments ( Solution )
      end if
      end associate !-- LA
    else
      call Show ( 'Laplacian_M not allocated', CONSOLE % ERROR )
      call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Solve_M', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine Solve_M


  subroutine CombineMoments ( P, Solution, T_Option )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( TimerForm ), pointer :: &
      T_CS, &
      T_LS, &
      T_C, &
      T_ES, &
      T_BS

    if ( .not. allocated ( P % Laplacian_M ) ) then
      call Show ( 'Laplacian_M not allocated', CONSOLE % ERROR )
      call Show ( 'Poisson_H_Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMoments', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( present ( T_Option ) ) then
      T_CS  =>  PROGRAM_HEADER % Timer &
                  ( Handle = P % iTimer_CS, &
                    Name = trim ( P % Name ) // '_ClrSltn', &
                    Level = T_Option % Level + 1 )
      T_LS  =>  PROGRAM_HEADER % Timer &
                  ( Handle = P % iTimer_LS, &
                    Name = trim ( P % Name ) // '_LclSltn', &
                    Level = T_Option % Level + 1 )
      T_C   =>  PROGRAM_HEADER % Timer &
                  ( Handle = P % iTimer_C, &
                    Name = trim ( P % Name ) // '_Crsng', &
                    Level = T_Option % Level + 1 )
      T_ES  =>  PROGRAM_HEADER % Timer &
                  ( Handle = P % iTimer_ES, &
                    Name = trim ( P % Name ) // '_ExchngSltn', &
                    Level = T_Option % Level + 1 )
      T_BS  =>  PROGRAM_HEADER % Timer &
                  ( Handle = P % iTimer_BS, &
                    Name = trim ( P % Name ) // '_BndrySltn', &
                    Level = T_Option % Level + 1 )
    else
      T_CS  =>  null ( )
      T_LS  =>  null ( )
      T_C   =>  null ( )
      T_ES  =>  null ( )
      T_BS  =>  null ( )
    end if

    call Show ( 'Combining Moments', P % IGNORABILITY + 4 )

    if ( associated ( T_CS ) ) call T_CS % Start ( )
    call Solution % Clear ( )
    if ( associated ( T_CS ) ) call T_CS % Stop ( )

    if ( associated ( T_LS ) ) call T_LS % Start ( )
    if ( allocated ( P % Laplacian_M ) ) &
      call P % CombineMomentsLocal ( Solution )
    if ( associated ( T_LS ) ) call T_LS % Stop ( )

    if ( associated ( T_C ) ) call T_C % Start ( )
    if ( allocated ( P % Coarsening ) ) &
      call P % Coarsening % Compute ( Solution )
    if ( associated ( T_C ) ) call T_C % Stop ( )

    if ( associated ( T_ES ) ) call T_ES % Start ( )
    call P % ExchangeSolution ( Solution )
    if ( associated ( T_ES ) ) call T_ES % Stop ( )

    if ( associated ( T_BS ) ) call T_BS % Start ( )
    call P % ApplyBoundarySolution ( Solution )
    if ( associated ( T_BS ) ) call T_BS % Stop ( )

  end subroutine CombineMoments


  subroutine CombineMomentsLocal ( P, Solution )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'CombineMomentsLocal', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine CombineMomentsLocal


  subroutine ExchangeSolution ( P, Solution )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'ExchangeSolution', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ExchangeSolution


  subroutine ApplyBoundarySolution ( P, Solution )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'ApplyBoundarySolution', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ApplyBoundarySolution


end module Poisson_H__Form
