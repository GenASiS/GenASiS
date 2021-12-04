module Poisson_ASCG__Form

  !-- Poisson_AtlasSingleChartGrid__Form

  use Basics
  use Manifolds
  use Fields
  use Laplacian_M_ASCG__Form
  use Poisson_H__Form

  implicit none
  private

  type, public, extends ( Poisson_H_Form ) :: Poisson_ASCG_Form
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      Solution_4D => null ( )
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      CombineMomentsLocal
    procedure, private, pass :: &
      ExchangeSolution
    procedure, private, pass :: &
      ApplyBoundarySolution
  end type Poisson_ASCG_Form

    private :: &
      CombineMoments_CGS_S_Kernel

    interface

      module subroutine CombineMoments_CGS_S_Kernel &
                          ( S, RM_R, RM_I, AF, RF_R, RF_I, DF, &
                            nC, oC, nE, nAM, oR, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          S  !-- Solution
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          RM_R, RM_I, &  !-- RadialMoment_Regular, _Irregular
          AF             !-- AngularFunction
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          RF_R, RF_I  !-- RadialFunction_Regular, _Irregular
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          DF  !-- DeltaFactor
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, oC  !-- nCells, oCell
        integer ( KDI ), intent ( in ) :: &
          nE, &   !-- nEquations
          nAM, &  !-- nAngularMoments
          oR      !-- oRadius
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine CombineMoments_CGS_S_Kernel

    end interface


    private :: &
      AssignSolutionPointer


contains


  subroutine Initialize &
               ( P, G, SolverType, MaxDegreeOption, nEquationsOption )

    class ( Poisson_ASCG_Form ), intent ( inout ) :: &
      P
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    character ( * ), intent ( in ) :: &
      SolverType
    integer ( KDI ), intent ( in ), optional :: &
      MaxDegreeOption, &
      nEquationsOption

    if ( P % Type == '' ) &
      P % Type = 'a Poisson_ASCG' 

     call P % Initialize_H &
           ( G, SolverType, MaxDegreeOption, nEquationsOption )

    P % Geometry  =>  G

    select case ( trim ( P % SolverType ) )
    case ( 'MULTIPOLE' )
      allocate ( Laplacian_M_ASCG_Form :: P % Laplacian_M )
      select type ( L => P % Laplacian_M )
      class is ( Laplacian_M_ASCG_Form )
        call L % Initialize ( G, P % MaxDegree, P % nEquations )
      end select !-- L
    case default
      call Show ( 'Solver type not recognized', CONSOLE % ERROR )
      call Show ( SolverType, 'Type', CONSOLE % ERROR )
      call Show ( 'Poisson_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine Initialize


  impure elemental subroutine Finalize ( P )

    type ( Poisson_ASCG_Form ), intent ( inout ) :: &
      P

    nullify ( P % Geometry )
    nullify ( P % Solution_4D )

  end subroutine Finalize


  subroutine CombineMomentsLocal ( P, Solution )

    class ( Poisson_ASCG_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution

    call Show ( 'Combining Moments Local', P % IGNORABILITY + 5 )

    select type ( L  =>  P % Laplacian_M )
      class is ( Laplacian_M_ASCG_Form )
    select type ( A  =>  P % Geometry % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )
    associate &
      ( SS  =>  Solution % Storage ( 1 ) )
    associate &
      (  nV => SS % nVariables, &
        iaS => SS % iaSelected )
 
    if ( nV  /=  L % nEquations ) then
      call Show ( 'Wrong number of variables in Solution', CONSOLE % ERROR )
      call Show ( 'Poisson_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iaS ( nV ) - iaS ( 1 ) + 1  /=  nV ) then
      call Show ( 'Solution variables must be contiguous', CONSOLE % ERROR )
      call Show ( 'Poisson_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call AssignSolutionPointer &
           ( SS % Value ( :, iaS ( 1 ) : iaS ( nV ) ), &
             C % nCellsBrick, C % nGhostLayers, L % nEquations, &
             P % Solution_4D )
    
    call SS % ReassociateHost ( AssociateVariablesOption = .false. )

    end associate !-- nV, etc.
    

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
      call CombineMoments_CGS_S_Kernel &
             ( P % Solution_4D, L % RadialMoment_R_3D, L % RadialMoment_I_3D, &
               L % AngularFunction_3D, L % RadialFunctions_R % Value, &
               L % RadialFunctions_I % Value, &
               L % DeltaFactor % Value ( :, 1 ), &
               C % nCellsBrick, C % nGhostLayers, &
               L % nEquations, L % nAngularMoments, &
               oR = ( C % iaBrick ( 1 ) - 1 ) * C % nCellsBrick ( 1 ), &
               UseDeviceOption = L % DeviceMemory )
    case default
      call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'Poisson_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select
    
    call SS % ReassociateHost ( AssociateVariablesOption = .true. )
    
    end associate !-- SS
    end associate !-- C

    class default 
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Poisson_ASCG_Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    class default 
      call Show ( 'Laplacian type not recognized', CONSOLE % ERROR )
      call Show ( 'Poisson_ASCG_Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- L

  end subroutine CombineMomentsLocal


  subroutine ExchangeSolution ( P, Solution )

    class ( Poisson_ASCG_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution

    call Solution % ExchangeGhostData ( )

  end subroutine ExchangeSolution


  subroutine ApplyBoundarySolution ( P, Solution )

    class ( Poisson_ASCG_Form ), intent ( inout ) :: &
      P
    class ( FieldSetForm ), intent ( inout ) :: &
      Solution

!    real ( KDR ), dimension ( :, :, : ), pointer :: &
!      SV_3D

    call Solution % ApplyBoundaryConditions ( )

    ! call Solution % Show ( )
    ! associate ( SC  =>  Solution % FieldSet_C ( 1 ) % Element )
    ! associate ( SV  =>  SC % Storage_FSC % Storage % Value )
    ! select type ( C  =>  SC % Chart )
    !  class is ( Chart_GS_Form )

    ! call C % SetFieldPointer ( SV ( :, 1 ), SV_3D )
    ! call Show ( SV_3D, '>>> ' // trim ( SC % Field ( 1 ) ) )  

    ! end select !-- C
    ! end associate !-- SV
    ! end associate !-- SC

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


end module Poisson_ASCG__Form
