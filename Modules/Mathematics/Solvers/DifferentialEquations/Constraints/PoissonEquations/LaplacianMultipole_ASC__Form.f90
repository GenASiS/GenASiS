module LaplacianMultipole_ASC__Form

  !-- LaplacianMultipole_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Template

  implicit none
  private

  type, public, extends ( LaplacianMultipoleTemplate ) :: &
    LaplacianMultipole_ASC_Form
      class ( ChartTemplate ), pointer :: &
        Chart => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      SetParameters
    procedure, private, pass :: &
      AllocateSolidHarmonics
    final :: &
      Finalize
  end type LaplacianMultipole_ASC_Form

contains


  subroutine Initialize ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipole_ASC' 

    call L % InitializeTemplate ( A, MaxDegree, nEquations )

call Show ( '>>> Aborting during development', CONSOLE % ERROR )
call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
call PROGRAM_HEADER % Abort ( )

!    call L % SetMomentStorage ( )

  end subroutine Initialize


  subroutine SetParameters ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    integer ( KDI ) :: &
      iL, &
      iM
    class ( GeometryFlatForm ), pointer :: &
      G

    L % nEquations  =   nEquations
    L % MaxDegree   =   MaxDegree
    L % MaxOrder    =   MaxDegree

    select type ( A )
    class is ( Atlas_SC_Form )

      L % Chart  =>  A % Chart

      G => A % Geometry ( )
      L % UseDevice = G % AllocatedDevice

    end select !-- A

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
       if ( C % nDimensions  <  3 ) &
         L % MaxOrder  =  0
       if ( C % nDimensions  <  2 ) &
         L % MaxDegree  =  0
    case default
      call Show ( 'CoordinateSystem not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'SetParameters', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'SetParameters', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C
    
    associate &
      ( L   => L % MaxDegree, &
        M   => L % MaxOrder, &
        nAM => L % nAngularMomentCells )
    nAM = 0
    do iM = 0, M
      do iL = iM, L
        nAM = nAM + 1
      end do
    end do
    end associate !-- L, etc.

    call Show ( L % MaxDegree, 'MaxDegree (l)', L % IGNORABILITY )
    call Show ( L % MaxOrder, 'MaxOrder (m)', L % IGNORABILITY )
    call Show ( L % nAngularMomentCells, 'nAngularMomentCells', &
                L % IGNORABILITY )
    call Show ( L % UseDevice, 'UseDevice', L % IGNORABILITY )

    nullify ( G )

  end subroutine SetParameters


  subroutine AllocateSolidHarmonics ( L )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    integer ( KDI ) :: &
      iSH

    allocate ( L % SolidHarmonics_1D ( 4 ) )
      !-- Current, Previous_1, Previous_2, PreviousDiagonal

    do iSH = 1, size ( L % SolidHarmonics_1D )

      associate ( SH => L % SolidHarmonics_1D ( iSH ) )

      select type ( C => L % Chart )
      class is ( Chart_SL_Template )

        associate ( nC => product ( C % nCells ) )
        call SH % Initialize ( [ nC, 4 ] )
               !-- RegularCos, IrregularCos, RegularSin, IrregularSin
        if ( L % UseDevice ) &
          call SH % AllocateDevice ( )
        end associate !-- nC

      class default
        call Show ( 'Chart type not supported', CONSOLE % ERROR )
        call Show ( 'AllocateSolidHarmonics', 'subroutine', CONSOLE % ERROR )
        call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C

      end associate !-- SH

    end do !-- iSH

  end subroutine AllocateSolidHarmonics


  impure elemental subroutine Finalize ( L )

    type ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    nullify ( L % Chart )

    call L % FinalizeTemplate ( )

  end subroutine Finalize


end module LaplacianMultipole_ASC__Form
