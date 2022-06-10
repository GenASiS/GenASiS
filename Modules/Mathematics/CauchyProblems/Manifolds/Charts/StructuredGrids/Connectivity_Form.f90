!-- Connectivity indexes connections (faces and edges) of segments (1D), &
!   quadrilaterals (2D), or hexahedra (3D).

module Connectivity_Form

  use Basics

  implicit none
  private

  type, public :: ConnectivityForm
    integer ( KDI ) :: &
      nDimensions, &
      nFaces, &
      nEdges, &
      nConnections
    integer ( KDI ), dimension ( 3 ) :: &
      iaInner = 0, &
      iaOuter = 0, &
      iaInnerInner = 0, &
      iaOuterInner = 0, &
      iaInnerOuter = 0, &
      iaOuterOuter = 0
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Show => Show_C
  end type ConnectivityForm

    integer ( KDI ), dimension ( 3 ), private, parameter :: &
      INNER = [ 1, 3, 5 ], &
      OUTER = [ 2, 4, 6 ], &
      INNER_INNER = [ 5, 9,  1 ], & 
      OUTER_INNER = [ 6, 10, 2 ], &  
      INNER_OUTER = [ 7, 11, 3 ], &
      OUTER_OUTER = [ 8, 12, 4 ]

    private :: &
      nBoundaryHypercubes, &
      Factorial

contains


  subroutine Initialize &
               ( C, nDimensions, IncludeFacesOption, IncludeEdgesOption )

    class ( ConnectivityForm ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    logical ( KDL ), intent ( in ), optional :: &
      IncludeFacesOption, &
      IncludeEdgesOption

    integer ( KDI ) :: &
      iD, jD, kD  !-- iDimension, etc.
    logical ( KDL ) :: &
      IncludeFaces, &
      IncludeEdges

    IncludeFaces = .true.
    if ( present ( IncludeFacesOption ) ) IncludeFaces = IncludeFacesOption

    IncludeEdges = .false.
    if ( present ( IncludeEdgesOption ) ) IncludeEdges = IncludeEdgesOption

    C % nDimensions = nDimensions
    
    associate &
      ( nD => nDimensions, &
        oF => 0, &            !-- oFace
        oE => C % nFaces )    !-- oEdge

    if ( IncludeFaces ) then
      do iD = 1, 3
        C % iaInner ( iD ) = oF + INNER ( iD )
        C % iaOuter ( iD ) = oF + OUTER ( iD )
      end do !-- iD
      C % nFaces = nBoundaryHypercubes ( nD, nD - 1 )  
    else !-- .not. IncludeFaces
      C % nFaces = 0
    end if !-- IncludeFaces

    if ( IncludeEdges ) then
      do iD = 1, 3
        jD = mod ( iD, 3 ) + 1
        kD = mod ( jD, 3 ) + 1
        if ( jD > nD .or. kD > nD ) cycle
        C % iaInnerInner ( iD ) = oE + INNER_INNER ( iD )
        C % iaOuterInner ( iD ) = oE + OUTER_INNER ( iD )
        C % iaInnerOuter ( iD ) = oE + INNER_OUTER ( iD )
        C % iaOuterOuter ( iD ) = oE + OUTER_OUTER ( iD )
      end do !-- iD
      C % nEdges = nBoundaryHypercubes ( nD, nD - 2 )
    else !-- .not. IncludeEdges
      C % nEdges = 0
    end if !-- IncludeEdges

    C % nConnections = C % nFaces + C % nEdges

    end associate !-- oF, etc.

  end subroutine Initialize


  subroutine Show_C ( C, IgnorabilityOption, nLeadingLinesOption, &
                      nTrailingLinesOption )

    class ( ConnectivityForm ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    call Show ( 'Connectivity', IgnorabilityOption, &
                nLeadingLinesOption = nLeadingLinesOption )
    call Show ( C % nDimensions, 'nDimensions', IgnorabilityOption )
    call Show ( C % nFaces, 'nFaces', IgnorabilityOption )
    call Show ( C % nEdges, 'nEdges', IgnorabilityOption )
    call Show ( C % nConnections, 'nConnections', IgnorabilityOption )
    call Show ( C % iaInner, 'iaInner', IgnorabilityOption )
    call Show ( C % iaOuter, 'iaOuter', IgnorabilityOption )
    call Show ( C % iaInnerInner, 'iaInnerInner', IgnorabilityOption )
    call Show ( C % iaOuterInner, 'iaOuterInner', IgnorabilityOption )
    call Show ( C % iaInnerOuter, 'iaInnerOuter', IgnorabilityOption )
    call Show ( C % iaOuterOuter, 'iaOuterOuter', IgnorabilityOption, &
                nTrailingLinesOption = nTrailingLinesOption )

  end subroutine Show_C


  function nBoundaryHypercubes ( N, M ) result ( nBH )

    !-- See Wikipedia "Hypercube"

    integer ( KDI ), intent ( in ) :: &
      N, M
    integer ( KDI ) :: &
      nBH

    if ( M >= 0 ) then
      nBH = ( 2 ** ( N - M ) ) &
            * Factorial ( N ) /  ( Factorial ( M ) * Factorial ( N - M ) )
    else
      nBH = 0
    end if

  end function nBoundaryHypercubes


  function Factorial ( N ) result ( F )

    integer ( KDI ), intent ( in ) :: &
      N
    integer ( KDI ) :: &
      F

    integer ( KDI ) :: &
      i

    if ( N == 0 ) then
      F = 1
    else
      F = product ( [ ( i, i = 1, N ) ] )
    end if

  end function Factorial


end module Connectivity_Form
