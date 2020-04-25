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
    final :: &
      Finalize
    procedure, private, pass :: &
      SetParametersAtlas
    procedure, private, pass :: &
      AllocateSolidHarmonics
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

  end subroutine Initialize


  impure elemental subroutine Finalize ( L )

    type ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    nullify ( L % Chart )

    call L % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetParametersAtlas ( L, A )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A

    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( A )
    class is ( Atlas_SC_Form )

      if ( A % nDimensions  <  3 ) &
        L % MaxOrder  =  0
      if ( A % nDimensions  <  2 ) &
        L % MaxDegree  =  0

      G  =>  A % Geometry ( )
      L % UseDevice  =  G % AllocatedDevice

      L % Chart  =>  A % Chart

    end select !-- A

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
       L % nRadialCells  =  C % nCells ( 1 )
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
    
    nullify ( G )

  end subroutine SetParametersAtlas


  subroutine AllocateSolidHarmonics ( L )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    integer ( KDI ) :: &
      iSH

    allocate ( L % SolidHarmonics_1D ( 4 ) )
      !-- 4: Current, Previous_1, Previous_2, PreviousDiagonal

    do iSH = 1, size ( L % SolidHarmonics_1D )

      associate ( SH => L % SolidHarmonics_1D ( iSH ) )

      select type ( C => L % Chart )
      class is ( Chart_SL_Template )

        associate ( nC => product ( C % nCells ) )
        call SH % Initialize ( [ nC, 4 ] )
               !-- 4: RegularCos, IrregularCos, RegularSin, IrregularSin
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


end module LaplacianMultipole_ASC__Form
