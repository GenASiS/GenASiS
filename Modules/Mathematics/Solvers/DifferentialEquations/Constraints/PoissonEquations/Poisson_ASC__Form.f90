module Poisson_ASC__Form

  !-- Poisson_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipoleOld_Template
  use LaplacianMultipoleOld_ASC__Form
  use LaplacianMultipole_ASC__Form
  use Poisson_Template

  implicit none
  private

  type, public, extends ( PoissonTemplate ) :: Poisson_ASC_Form
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      Solution => null ( )
    class ( Atlas_SC_Form ), pointer :: &
      Atlas => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SolveOld
    procedure, private, pass :: &
      CombineMomentAtlas
    procedure, private, pass :: &
      ExchangeSolution
    procedure, private, pass :: &
      ApplyBoundarySolution
  end type Poisson_ASC_Form

    private :: &
      SolveMultipole_CSL_Old

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!    private :: &
    public :: &
      CombineMoment_CSL_S_Kernel, &
      SolveCells_CSL_Kernel, &
      AssembleSolutionKernel

    interface

      module subroutine CombineMoment_CSL_S_Kernel &
                          ( S, SH_RC, SH_IC, SH_RS, SH_IS, R_C, &
                            M_RC,  M_IC,  M_RS,  M_IS, R_I, &
                            IsFirstShell, IsLastShell, Delta_M_FourPi, &
                            nC, oC, nE, oR, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          S  !-- Solution
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          SH_RC, SH_IC, SH_RS, SH_IS, &  !-- SolidHarmonics
          R_C
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          M_RC, M_IC, M_RS, M_IS  !-- MyMoments
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          R_I
        logical ( KDL ), intent ( in ) :: &
          IsFirstShell, IsLastShell
        real ( KDR ), intent ( in ) :: &
          Delta_M_FourPi
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, oC
        integer ( KDI ), intent ( in ) :: &
          nE, &  !-- nEquations
          oR     !-- oRadius
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine

      module subroutine SolveCells_CSL_Kernel &
                          ( Solution, CoordinateSystem, IsProperCell, &
                            M_RC, M_IC, M_RS, M_IS, Center, Origin, Delta, &
                            RadialEdge, FourPi, iaSolution, MaxDegree, &
                            MaxOrder, nDimensions, nCells, nEquations, &
                            nRadialCells, nAngularMomentCells, GridError, &
                            SH_RC, SH_IC, SH_RS, SH_IS )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          Solution
        character ( * ), intent ( in ) :: &
          CoordinateSystem
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          IsProperCell
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          M_RC, M_IC, &
          M_RS, M_IS
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          Center
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Origin, &
          Delta, &
          RadialEdge
        real ( KDR ), intent ( in ) :: &
          FourPi
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSolution
        integer ( KDI ), intent ( in ) :: &
          MaxDegree, &
          MaxOrder, &
          nDimensions, &
          nCells, &
          nEquations, &
          nRadialCells, &
          nAngularMomentCells
        logical ( KDL ), intent ( out ) :: &
          GridError
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          SH_RC, SH_IC, &  !-- SolidHarmonics
          SH_RS, SH_IS  
      end subroutine SolveCells_CSL_Kernel

      module subroutine AssembleSolutionKernel &
                          ( S, M_R_In, M_I_In, M_R_Out, M_I_Out, SH_R, SH_I, &
                            Delta, R_In, R_Out, R, nA )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          S
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_R_In, M_I_In, &
          M_R_Out, M_I_Out, &
          SH_R, SH_I, &
          Delta
        real ( KDR ), intent ( in ) :: &
          R_In, R_Out, R
        integer ( KDI ), intent ( in ) :: &
          nA
      end subroutine AssembleSolutionKernel
      
    end interface


    private :: &
      AssignSolutionPointer

contains


  subroutine Initialize &
               ( P, A, SolverType, MaxDegreeOption, nEquationsOption )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      SolverType
    integer ( KDI ), intent ( in ), optional :: &
      MaxDegreeOption, &
      nEquationsOption

    if ( P % Type == '' ) &
      P % Type = 'a Poisson_ASC' 

    call P % InitializeTemplate &
           ( A, SolverType, MaxDegreeOption, nEquationsOption )

    P % Atlas => A

    select case ( trim ( P % SolverType ) )
    case ( 'MULTIPOLE_OLD' )
      allocate ( LaplacianMultipoleOld_ASC_Form :: P % LaplacianMultipoleOld )
      select type ( L => P % LaplacianMultipoleOld )
      class is ( LaplacianMultipoleOld_ASC_Form )
        call L % Initialize ( A, P % MaxDegree, P % nEquations )
      end select !-- L
    case ( 'MULTIPOLE' )
      allocate ( LaplacianMultipole_ASC_Form :: P % LaplacianMultipole )
      select type ( L => P % LaplacianMultipole )
      class is ( LaplacianMultipole_ASC_Form )
        call L % Initialize ( A, P % MaxDegree, P % nEquations )
      end select !-- L
    case default
      call Show ( 'Solver type not supported', CONSOLE % ERROR )
      call Show ( SolverType, 'Type', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine Initialize


  impure elemental subroutine Finalize ( P )

    type ( Poisson_ASC_Form ), intent ( inout ) :: &
      P

    nullify ( P % Atlas )
    nullify ( P % Solution )

    call P % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SolveOld ( P, Solution, Source )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution
    class ( FieldAtlasTemplate ), intent ( in ) :: &
      Source

    class ( Storage_CSL_Form ), pointer :: &
      Source_CSL, &
      Solution_CSL

    select type ( Source )
    class is ( Storage_ASC_Form )
    Source_CSL => Source % Storage_CSL ( )

    select type ( Solution )
    class is ( Storage_ASC_Form )
    Solution_CSL => Solution % Storage_CSL ( )

    select type ( C => P % Atlas % Chart )
    class is ( Chart_SLD_Form )
      call SolveMultipole_CSL_Old ( P, C, Solution_CSL, Source_CSL )
    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end select !-- Solution
    end select !-- Source

  end subroutine SolveOld


  subroutine SolveMultipole_CSL_Old ( P, C, Solution_CSL, Source_CSL )
 
    type ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( Storage_CSL_Form ), intent ( inout ), target :: &
      Solution_CSL
    class ( Storage_CSL_Form ), intent ( in ), target :: &
      Source_CSL

    logical ( KDL ) :: &
      GridError
    type ( TimerForm ), pointer :: &
      Timer_CM, &
      Timer_ES, &
      Timer_BS
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( StorageForm ), pointer :: &
      Source, &
      Solution
    
    Source    => Source_CSL % Storage ( )
    Solution  => Source_CSL % Storage ( )

    Timer_CM  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerCombineMoments )
    Timer_ES  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerExchangeSolution )
    Timer_BS  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerBoundarySolution )

    call Show ( 'Poisson solve, multipole old', P % IGNORABILITY + 2 )
    call Show ( P % Name, 'Name', P % IGNORABILITY + 2 )

    associate ( L  =>  P % LaplacianMultipoleOld )

    call L % ComputeMoments ( Source )

    G => C % Geometry ( )

    if ( associated ( Timer_CM ) ) call Timer_CM % Start ( )
     call SolveCells_CSL_Kernel &
            ( Solution % Value, C % CoordinateSystem, C % IsProperCell, &
              L % M_RC, L % M_IC, L % M_RS, L % M_IS, &
              G % Value ( :, G % CENTER_U ( 1 ) : G % CENTER_U ( 3 ) ), &
              L % Origin, L % Delta, L % RadialEdge, 4.0_KDR * CONSTANT % PI, &
              Solution % iaSelected, L % MaxDegree, L % MaxOrder, &
              C % nDimensions, G % nValues, L % nEquations, L % nRadialCells, &
              L % nAngularMomentCells, GridError, &
              L % SolidHarmonic_RC, L % SolidHarmonic_IC, &
              L % SolidHarmonic_RS, L % SolidHarmonic_IS )
     if ( GridError ) &
       call PROGRAM_HEADER % Abort ( )
    if ( associated ( Timer_CM ) ) call Timer_CM % Stop ( )

    if ( associated ( Timer_ES ) ) call Timer_ES % Start ( )
    call C % ExchangeGhostData ( Solution_CSL )
    if ( associated ( Timer_ES ) ) call Timer_ES % Stop ( )

    if ( associated ( Timer_BS ) ) call Timer_BS % Start ( )
    call P % Atlas % ApplyBoundaryConditionsFaces ( Solution )
    if ( associated ( Timer_BS ) ) call Timer_BS % Stop ( )

    end associate !-- L

    nullify ( G )

  end subroutine SolveMultipole_CSL_Old


  subroutine CombineMomentAtlas &
               ( P, Solution, Delta_M_FourPi, iA, iSH_0 )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution
    real ( KDR ), intent ( in ) :: &
      Delta_M_FourPi
    integer ( KDI ), intent ( in ) :: &
      iA, &  
      iSH_0  

    logical ( KDL ) :: &
      IsFirstShell, IsLastShell
    class ( StorageForm ), pointer :: &
      Solution_S

    select type ( L => P % LaplacianMultipole )
    class is ( LaplacianMultipole_ASC_Form )

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

    select type ( Solution )
    class is ( Storage_ASC_Form )
    Solution_S => Solution % Storage ( )

    associate &
      (  nV => Solution_S % nVariables, &
        iaS => Solution_S % iaSelected )
 
    if ( nV /= L % nEquations ) then
      call Show ( 'Wrong number of variables in Solution', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iaS ( nV ) - iaS ( 1 ) + 1  /=  nV ) then
      call Show ( 'Solution variables must be contiguous', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call AssignSolutionPointer &
           ( Solution_S % Value ( :, iaS ( 1 ) : iaS ( nV ) ), &
             C % nCellsBrick, C % nGhostLayers, L % nEquations, P % Solution )

    end associate !-- nV, etc.

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
      IsFirstShell  =  ( C % iaBrick ( 1 )  ==  1 )
      IsLastShell  =  ( C % iaBrick ( 1 )  ==  C % nBricks ( 1 ) )
      call CombineMoment_CSL_S_Kernel &
             ( P % Solution, &
               L % SolidHarmonic_RC ( :, :, :, iSH_0 ), &
               L % SolidHarmonic_IC ( :, :, :, iSH_0 ), &
               L % SolidHarmonic_RS ( :, :, :, iSH_0 ), &
               L % SolidHarmonic_IS ( :, :, :, iSH_0 ), &
               L % Radius, &
               L % Moment_RC ( :, :, iA ), L % Moment_IC ( :, :, iA ), &
               L % Moment_RS ( :, :, iA ), L % Moment_IS ( :, :, iA ), &
               L % RadialEdges % Value ( :, 1 ), IsFirstShell, IsLastShell, &
               Delta_M_FourPi, C % nCellsBrick, C % nGhostLayers, &
               L % nEquations, &
               oR = ( C % iaBrick ( 1 ) - 1 ) * C % nCellsBrick ( 1 ), &
               UseDeviceOption = L % UseDevice )
    case default
      call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    class default
      call Show ( 'Solution type not supported', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Source

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    class default 
      call Show ( 'Laplacian type not supported', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- L

    nullify ( Solution_S )

  end subroutine CombineMomentAtlas


  subroutine ExchangeSolution ( P, Solution )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution

    class ( StorageForm ), pointer :: &
      Solution_S
    class ( Storage_CSL_Form ), pointer :: &
      Solution_CSL

    select type ( Solution )
    class is ( Storage_ASC_Form )
    Solution_S    => Solution % Storage ( )
    Solution_CSL  => Solution % Storage_CSL ( )

    select type ( C => P % Atlas % Chart )
    class is ( Chart_SLD_Form )
      call Solution_S % UpdateHost ( )
      call C % ExchangeGhostData ( Solution_CSL )
      call Solution_S % UpdateDevice ( )
    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ExchangeSolution', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    class default
      call Show ( 'Solution type not supported', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ExchangeSolution', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Source

    nullify ( Solution_S )
    
  end subroutine ExchangeSolution


  subroutine ApplyBoundarySolution ( P, Solution )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution

    class ( StorageForm ), pointer :: &
      Solution_S

    select type ( Solution )
    class is ( Storage_ASC_Form )
    Solution_S => Solution % Storage ( )

    call P % Atlas % ApplyBoundaryConditionsFaces ( Solution_S )

    class default
      call Show ( 'Solution type not supported', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ExchangeSolution', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Source

    nullify ( Solution_S )
    
  end subroutine ApplyBoundarySolution


  subroutine AssignSolutionPointer ( S_2D, nC, nG, nE, S_4D )

    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
      S_2D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nC, &  !-- nCellsBrick
      nG     !-- nGhostLayers
    integer ( KDI ), intent ( in ) :: &
      nE  !-- nEquations
    real ( KDR ), dimension ( :, :, :, : ), intent ( out ), pointer :: &
      S_4D

    associate &
      ( n1  =>  nC ( 1 )  +  2 * nG ( 1 ), &
        n2  =>  nC ( 2 )  +  2 * nG ( 2 ), &
        n3  =>  nC ( 3 )  +  2 * nG ( 3 ) )

    S_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : nE )  =>  S_2D

    end associate !-- n1, etc.

  end subroutine AssignSolutionPointer


end module Poisson_ASC__Form
