module LaplacianMultipole_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: LaplacianMultipoleForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nRadialCells = 0, &
      nAngularMomentCells = 0, &
      MaxMoment
    real ( KDR ), dimension ( 3 ) :: &
      Origin
    real ( KDR ), dimension ( : ), allocatable :: &
      RadialEdge, &
      SolidHarmonic_RC, SolidHarmonic_RS, &  !-- Regular Cos, Sin
      SolidHarmonic_IC, SolidHarmonic_IS, &  !-- Irregular Cos, Sin
      Delta
    character ( LDF ) :: &
      Name = ''
    class ( Atlas_SC_Form ), pointer :: &
      Atlas => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type LaplacianMultipoleForm

    private :: &
      SetRadialEdge

contains


  subroutine Initialize ( LM, Atlas, MaxMoment )

    class ( LaplacianMultipoleForm ), intent ( inout ) :: &
      LM
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      Atlas
    integer ( KDI ), intent ( in ) :: &
      MaxMoment

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iL, &
      iM

    LM % IGNORABILITY = Atlas % IGNORABILITY
    LM % Name = 'Laplacian_' // trim ( Atlas % Name )

    call Show ( 'Initializing a LaplacianMultipole', LM % IGNORABILITY )
    call Show ( LM % Name, 'Name', LM % IGNORABILITY )
    call Show ( MaxMoment, 'MaxMoment', LM % IGNORABILITY )

    LM % MaxMoment = MaxMoment
    LM % Atlas => Atlas

    LM % nAngularMomentCells = ( 1 + MaxMoment ) * ( 2 + MaxMoment ) / 2 
    call Show ( LM % nAngularMomentCells, 'nAngularMomentCells', &
                LM % IGNORABILITY )

    allocate ( LM % SolidHarmonic_RC ( LM % nAngularMomentCells ) )
    allocate ( LM % SolidHarmonic_RS ( LM % nAngularMomentCells ) )
    allocate ( LM % SolidHarmonic_IC ( LM % nAngularMomentCells ) )
    allocate ( LM % SolidHarmonic_IS ( LM % nAngularMomentCells ) )
    allocate ( LM % Delta ( LM % nAngularMomentCells ) )

    associate ( L => MaxMoment )
    iV = 0
    do iM = 0, L
      do iL = iM, L
        iV = iV + 1
        if ( iM == 0 ) then
          LM % Delta ( iV ) = 1.0_KDR
        else
          LM % Delta ( iV ) = 2.0_KDR
        end if
      end do
    end do
    end associate !-- L

    call SetRadialEdge ( LM )

  end subroutine Initialize


  impure elemental subroutine Finalize ( LM )

    type ( LaplacianMultipoleForm ), intent ( inout ) :: &
      LM

    nullify ( LM % Atlas )

    if ( allocated ( LM % Delta ) ) &
      deallocate ( LM % Delta )
    if ( allocated ( LM % SolidHarmonic_IS ) ) &
      deallocate ( LM % SolidHarmonic_IS )
    if ( allocated ( LM % SolidHarmonic_IC ) ) &
      deallocate ( LM % SolidHarmonic_IC )
    if ( allocated ( LM % SolidHarmonic_RS ) ) &
      deallocate ( LM % SolidHarmonic_RS )
    if ( allocated ( LM % SolidHarmonic_RC ) ) &
      deallocate ( LM % SolidHarmonic_RC )
    if ( allocated ( LM % RadialEdge ) ) &
      deallocate ( LM % RadialEdge )

    if ( LM % Name == '' ) return

    call Show ( 'Finalizing a LaplacianMultipole', LM % IGNORABILITY )
    call Show ( LM % Name, 'Name', LM % IGNORABILITY )
    
  end subroutine Finalize


  subroutine SetRadialEdge ( LM )

    type ( LaplacianMultipoleForm ), intent ( inout ) :: &
      LM

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ), dimension ( : ), allocatable :: &
      Center, &
      Width

    select type ( C => LM % Atlas % Chart )
    class is ( Chart_SL_Template )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )

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

    case default
      call Show ( 'CoordinateSystem not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'SetRadialEdge', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_Form', 'module', CONSOLE % ERROR )
    end select !-- CoordinateSystem

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'SetRadialEdge', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_Form', 'module', CONSOLE % ERROR )    
    end select !-- C

  end subroutine SetRadialEdge


end module LaplacianMultipole_Form
