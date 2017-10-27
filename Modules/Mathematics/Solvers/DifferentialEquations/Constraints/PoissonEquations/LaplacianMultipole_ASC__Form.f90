module LaplacianMultipole_ASC__Form

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
    final :: &
      Finalize
    procedure, private, pass :: &
      ComputeMomentsLocal
  end type LaplacianMultipole_ASC_Form

    private :: &
      SetRadialEdgeSpherical

contains


  subroutine Initialize ( LM, A, MaxDegree )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree

    if ( LM % Type == '' ) &
      LM % Type = 'a LaplacianMultipole_ASC' 

    call LM % InitializeTemplate ( A, MaxDegree )

    select case ( trim ( A % Chart % CoordinateSystem ) )
    case ( 'SPHERICAL' )
      call SetRadialEdgeSpherical ( LM )
      call LM % SetMomentStorage ( )
    case default
      call Show ( 'CoordinateSystem not supported', CONSOLE % ERROR )
      call Show ( A % Chart % CoordinateSystem, 'CoordinateSystem', &
                  CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
    end select !-- CoordinateSystem

  end subroutine Initialize


  impure elemental subroutine Finalize ( LM )

    type ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM

    nullify ( LM % Chart )

    call LM % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetParameters ( LM, A, MaxDegree )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree

    integer ( KDI ) :: &
      iL, &
      iM

    LM % MaxDegree  =   MaxDegree
    LM % MaxOrder   =   MaxDegree

    select type ( A )
    class is ( Atlas_SC_Form )
      LM % Chart  =>  A % Chart
    end select !-- A

    select type ( C => LM % Chart )
    class is ( Chart_SL_Template )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
       if ( C % nDimensions  <  3 ) &
         LM % MaxOrder  =  0
       if ( C % nDimensions  <  2 ) &
         LM % MaxDegree  =  0
    case default
      call Show ( 'CoordinateSystem not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'SetParameters', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
    end select !-- CoordinateSystem

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'SetParameters', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
    end select !-- C
    
    associate &
      ( L   => LM % MaxDegree, &
        M   => LM % MaxOrder, &
        nAM => LM % nAngularMomentCells )
    nAM = 0
    do iM = 0, M
      do iL = iM, L
        nAM = nAM + 1
      end do
    end do
    end associate !-- L, etc.

    call Show ( LM % MaxDegree, 'MaxDegree (l)', LM % IGNORABILITY )
    call Show ( LM % MaxOrder, 'MaxOrder (m)', LM % IGNORABILITY )
    call Show ( LM % nAngularMomentCells, 'nAngularMomentCells', &
                LM % IGNORABILITY )

  end subroutine SetParameters


  subroutine ComputeMomentsLocal ( LM, Source )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM
    type ( VariableGroupForm ), intent ( in ) :: &
      Source !-- array over levels    

    integer ( KDI ) :: &
      iC
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( C => LM % Chart )
    class is ( Chart_SL_Template )

      G => C % Geometry ( )

      !$OMP parallel do private ( iC )
      do iC = 1, G % nValues
        if ( .not. C % IsProperCell ( iC ) ) &
          cycle
        call LM % ComputeSolidHarmonics &
               ( C % CoordinateSystem, &
                 G % Value ( iC, G % CENTER ( 1 ) : G % CENTER ( 3 ) ), &
                 C % nDimensions )
        call LM % ComputeMomentContributions &
               ( C % CoordinateSystem, &
                 G % Value ( iC, G % CENTER ( 1 ) : G % CENTER ( 3 ) ), &
                 G % Value ( iC, G % WIDTH ( 1 ) : G % WIDTH ( 3 ) ), &
                 G % Value ( iC, G % VOLUME_JACOBIAN ), &
                 Source % Value ( iC, 1 ), &
                 C % nDimensions )
      end do
      !$OMP end parallel do

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'ComputeMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
    end select !-- C
    
    nullify ( G )

  end subroutine ComputeMomentsLocal


  subroutine SetRadialEdgeSpherical ( LM )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ), dimension ( : ), allocatable :: &
      Center, &
      Width

    select type ( C => LM % Chart )
    class is ( Chart_SL_Template )

      LM % Origin        =  0.0_KDR
      LM % nRadialCells  =  C % nCells ( 1 )

      associate ( nRC  =>  LM % nRadialCells )
      allocate ( Center ( nRC ) )
      allocate ( Width ( nRC ) )
      allocate ( LM % RadialEdge ( nRC + 1 ) )

      call C % SetGeometryCell &
             ( Width, Center, nC = nRC, nGL = 0, iD = 1, iaF = 1 )

      LM % RadialEdge ( 1 )  =  0.0_KDR
      do iC = 1, nRC
        LM % RadialEdge ( iC + 1 )  =  LM % RadialEdge ( iC )  +  Width ( iC )
      end do !-- iC
      
      end associate !-- nRC

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'SetRadialEdgeSpherical', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC_Form', 'module', CONSOLE % ERROR )    
    end select !-- C

  end subroutine SetRadialEdgeSpherical


end module LaplacianMultipole_ASC__Form
