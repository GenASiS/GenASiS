module Poisson_ASC__Form

  !-- Poisson_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_ASC__Form
  use Poisson_Template

  implicit none
  private

  type, public, extends ( PoissonTemplate ) :: Poisson_ASC_Form
    class ( Atlas_SC_Form ), pointer :: &
      Atlas => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Solve
    final :: &
      Finalize
  end type Poisson_ASC_Form

    private :: &
      SolveMultipole_CSL

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


  subroutine Solve ( P, Solution, Source )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution
    class ( FieldAtlasTemplate ), intent ( in ) :: &
      Source

    class ( VariableGroupForm ), pointer :: &
      Source_VG, &
      Solution_VG

    select type ( Source )
    class is ( Storage_ASC_Form )
    Source_VG => Source % Storage ( )

    select type ( Solution )
    class is ( Storage_ASC_Form )
    Solution_VG => Solution % Storage ( )

    select case ( trim ( P % SolverType ) )
    case ( 'MULTIPOLE' )
      select type ( C => P % Atlas % Chart )
      class is ( Chart_SL_Template )
        call SolveMultipole_CSL ( P, Solution_VG, C, Source_VG )
      class default
        call Show ( 'Chart type not supported', CONSOLE % ERROR )
        call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
        call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    case default
      call Show ( 'Solver type not supported', CONSOLE % ERROR )
      call Show ( P % SolverType, 'Type', CONSOLE % ERROR )
      call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- SolverType

    end select !-- Solution
    end select !-- Source

  end subroutine Solve


  impure elemental subroutine Finalize ( P )

    type ( Poisson_ASC_Form ), intent ( inout ) :: &
      P

    call P % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SolveMultipole_CSL ( P, Solution, C, Source )
 
    type ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( VariableGroupForm ), intent ( inout ) :: &
      Solution
    class ( Chart_SL_Template ), intent ( in ) :: &
      C
    class ( VariableGroupForm ), intent ( in ) :: &
      Source

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iR, &  !-- iRadius
      iE     !-- iEquation
    real ( KDR ) :: &
      R, &  !-- Radius
      S     !-- Solution
    real ( KDR ), dimension ( : ), allocatable :: &
      Zero
    class ( GeometryFlatForm ), pointer :: &
      G

    call Show ( 'Poisson solve, multipole', P % IGNORABILITY + 2 )
    call Show ( P % Name, 'Name', P % IGNORABILITY + 2 )

    associate ( L => P % LaplacianMultipole )

    ! !-- Compute potential on Interior leaf cells

    ! call LM % SetRadialGrid ( Source )
    ! call LM % ComputeMoments ( Source )

    ! call Show ( 'Computing potential', PEC % IGNORABILITY )

    ! do iL = 1, C % nLevels
    !   !-- Assume the relevant variable is the first in the source VG.
    !   associate ( iSltn => Solution ( iL ) % Selected ( 1 ) )
    !   associate &
    !     ( S_Interior => C % Level ( iL ) % Submesh ( SUBMESH % INTERIOR ), &
    !       G => C % Level ( iL ) % Geometry, &
    !       A => C % Level ( iL ) % Adaption, &
    !       S => Source ( iL ) % Value ( :, 1 ), &
    !       Phi => Solution ( iL ) % Value ( :, iSltn ) )

    !   call Show ( iL, 'iLevel' )

    !   do iC = S_Interior % oCell + 1, &
    !           S_Interior % oCell + S_Interior % nProperCells

    !     if ( A % Value ( iC, A % INTERIOR_PARENT ) > 0.0_KDR ) cycle

    !     call LM % ComputeSolidHarmonics ( G % Value ( iC, G % CENTER ) )
    !     call ComputePotential &
    !            ( LM, G % Value ( iC, G % CENTER ), G % Value ( iC, G % WIDTH ), &
    !              S ( iC ) * G % Value ( iC, G % VOLUME ), &
    !              Phi ( iC ) )

    !   end do

    !   end associate !-- S_Interior, etc.
    !   end associate !-- iSltn
    ! end do

    ! !-- Restrict interior values and ghost exchange

    ! if ( present ( iLevelOption ) ) then
    !   !-- argument of Sort must be intent ( inout )
    !   allocate ( iLevel ( size ( iLevelOption ) ) )
    !   iLevel = iLevelOption
    !   call Sort ( iLevel )
    !   do iL = size ( iLevel ), 1, -1
    !     associate ( iLvl => iLevel ( iL ) )
    !     if ( iLvl < C % nLevels ) &
    !       call C % Restrict &
    !              ( Solution ( iLvl + 1 ), Solution ( iLvl ), &
    !                Interior = .true., Exterior = .false., iLevel = iLvl + 1 )
    !     call C % Level ( iLvl ) % StartGhostExchange &
    !            ( Solution ( iLvl : iLvl ) )
    !     call C % Level ( iLvl ) % FinishGhostExchange ( )
    !     end associate !-- iLvl
    !   end do
    ! end if

    call L % ComputeMoments ( Source )

    G => C % Geometry ( )

    allocate ( Zero ( L % nAngularMomentCells ) )
    Zero = 0.0_KDR

    do iE = 1, L % nEquations

      !$OMP parallel do private ( iC, iR, R, S )
      do iC = 1, G % nValues

        if ( .not. C % IsProperCell ( iC ) ) &
          cycle

        call L % ComputeSolidHarmonics &
               ( C % CoordinateSystem, &
                 G % Value ( iC, G % CENTER ( 1 ) : G % CENTER ( 3 ) ), &
                 C % nDimensions, R, iR )

        S = 0.0_KDR

        if ( iR > 1 .and. iR < L % nRadialCells ) then
          call SolveKernel &
                 ( S, &
                   L % MRC ( :, iR - 1, iE ), L % MIC ( :, iR, iE ), &
                   L % MRC ( :, iR, iE ), L % MIC ( :, iR + 1, iE ), &
                   L % SolidHarmonic_RC ( : ), L % SolidHarmonic_IC ( : ), &
                   L % Delta, L % RadialEdge ( iR ), &
                   L % RadialEdge ( iR + 1 ), R )
        else if ( iR == 1 ) then
          call SolveKernel &
                 ( S, &
                   Zero, L % MIC ( :, iR, iE ), &
                   L % MRC ( :, iR, iE ), L % MIC ( :, iR + 1, iE ), &
                   L % SolidHarmonic_RC ( : ), L % SolidHarmonic_IC ( : ), &
                   L % Delta, L % RadialEdge ( iR ), &
                   L % RadialEdge ( iR + 1 ), R )
        else if ( iR == L % nRadialCells ) then
          call SolveKernel &
                 ( S, &
                   L % MRC ( :, iR - 1, iE ), L % MIC ( :, iR, iE ), &
                   L % MRC ( :, iR, iE ), Zero, &
                   L % SolidHarmonic_RC ( : ), L % SolidHarmonic_IC ( : ), &
                   L % Delta, L % RadialEdge ( iR ), &
                   L % RadialEdge ( iR + 1 ), R )
        end if !-- iR

        if ( L % MaxOrder > 0 ) then
          if ( iR > 1 .and. iR < L % nRadialCells ) then
            call SolveKernel &
                   ( S, &
                     L % MRS ( :, iR - 1, iE ), L % MIS ( :, iR, iE ), &
                     L % MRS ( :, iR, iE ), L % MIS ( :, iR + 1, iE ), &
                     L % SolidHarmonic_RS ( : ), L % SolidHarmonic_IS ( : ), &
                     L % Delta, L % RadialEdge ( iR ), &
                     L % RadialEdge ( iR + 1 ), R )
          else if ( iR == 1 ) then
            call SolveKernel &
                   ( S, &
                     Zero, L % MIS ( :, iR, iE ), &
                     L % MRS ( :, iR, iE ), L % MIS ( :, iR + 1, iE ), &
                     L % SolidHarmonic_RS ( : ), L % SolidHarmonic_IS ( : ), &
                     L % Delta, L % RadialEdge ( iR ), &
                     L % RadialEdge ( iR + 1 ), R )
          else if ( iR == L % nRadialCells ) then
            call SolveKernel &
                   ( S, &
                     L % MRS ( :, iR - 1, iE ), L % MIS ( :, iR, iE ), &
                     L % MRS ( :, iR, iE ), Zero, &
                     L % SolidHarmonic_RS ( : ), L % SolidHarmonic_IS ( : ), &
                     L % Delta, L % RadialEdge ( iR ), &
                     L % RadialEdge ( iR + 1 ), R )
          end if !-- iR
        end if !-- MaxOrder

        Solution % Value ( iC, Solution % iaSelected ( iE ) )  &
          =  - S / ( 4.0_KDR  *  CONSTANT % PI )

      end do !-- iC
      !$OMP end parallel do
    end do !-- iE

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )

    case default
      call Show ( 'CoordinateSystem not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', &
                  CONSOLE % ERROR )
      call Show ( 'SolveMultipole_CSL', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    end associate !-- L

    nullify ( G )

  end subroutine SolveMultipole_CSL


  subroutine SolveKernel &
               ( S, M_R_In, M_I_In, M_R_Out, M_I_Out, SH_R, SH_I, Delta, &
                 R_In, R_Out, R )

    real ( KDR ), intent ( inout ) :: &
      S
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_R_In, M_I_In, &
      M_R_Out, M_I_Out, &
      SH_R, SH_I, &
      Delta
    real ( KDR ), intent ( in ) :: &
      R_In, R_Out, R

    integer ( KDI ) :: &
      iA  !-- iAngular
    real ( KDR ) :: &
      F, &  !-- Fraction
      M_R, M_I

    F  =  ( R - R_In ) / ( R_Out - R_in )

    do iA = 1, size ( Delta )

      M_R  =  F  *  M_R_In ( iA )   +   ( 1.0_KDR - F )  *  M_R_Out ( iA )
      M_I  =  F  *  M_I_In ( iA )   +   ( 1.0_KDR - F )  *  M_I_Out ( iA )

      S  =  S   +   Delta ( iA )  &
                    *  ( M_R  *  SH_I ( iA )  +  M_I  *  SH_R ( iA ) )

    end do

  end subroutine SolveKernel


end module Poisson_ASC__Form
